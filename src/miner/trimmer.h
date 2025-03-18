// Copyright 2025 The MWC Developers
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef TRIMMER_H
#define TRIMMER_H

#include "row_builder.h"
#include "bit_streamer.h"

template<uint8_t EDGE_BITS, uint8_t BUCKET_BITS>
static void bucket_hashes1(std::vector<MTL::Buffer *> *hash_data,
        std::vector<BitStreamer<EDGE_BITS-BUCKET_BITS>> * result,
        std::atomic_int32_t * tasks) {

    uint32_t BUCKETS_NUM = 1 << BUCKET_BITS;

    while (true) {
        int32_t task_idx = tasks->fetch_add(1, std::memory_order_relaxed);
        if (task_idx >= hash_data->size())
            break;

        MTL::Buffer * hashes = (*hash_data)[task_idx];
        uint32_t sz = hashes->length() / sizeof(uint32_t);
        uint32_t * hash = (uint32_t *) hashes->contents();
        uint32_t * hashLim = hash + sz;

        uint32_t RES_HASH_MASK = (1 << (EDGE_BITS - BUCKET_BITS)) - 1;

        while (hash < hashLim) {
            if (*hash==UINT_MAX)
                break;

            assert(*hash < uint64_t(1) << EDGE_BITS );
            uint32_t bucket_idx = (*hash) >> (EDGE_BITS - BUCKET_BITS);
            assert(bucket_idx < BUCKETS_NUM);

            (*result)[bucket_idx].encode_number((*hash) & RES_HASH_MASK);
            hash++;
        }

        // Data is processed, release it. Caller expected that
        hashes->release();
    }

    // BitStreamer needs to be finalized
    for (int i = 0; i < BUCKETS_NUM; ++i) {
        (*result)[i].write_flush();
    }

}

template<uint8_t EDGE_BITS, uint ELEMENT_SIZE, uint8_t BUCKET_BITS, uint8_t MASK_BITS,
            uint32_t BIT_PACKER_INDEX_INTERVAL>
static void bucket_hashes2( std::vector<std::vector<BitStreamer<EDGE_BITS-BUCKET_BITS>>> * nonces_step1,
            RowBuilder<EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> *row,
            TrimNonces * trim_nonces,
            std::atomic_int32_t * tasks)
{
    uint32_t BUCKETS_NUM = 1 << BUCKET_BITS;
    uint32_t BUCKET_SIZE = 1 << (EDGE_BITS - BUCKET_BITS*2);
    const uint64_t edge_mask = (1 << MASK_BITS) - 1;

    while (true) {
        int32_t task_idx = tasks->fetch_add(1, std::memory_order_relaxed);
        if (task_idx >= BUCKETS_NUM)
            break;

        // Let's split hashes by buckets first, then process each bucket one by one
        BitStreamer<EDGE_BITS-BUCKET_BITS*2> nonce_by_bucket[BUCKETS_NUM];
        uint32_t HASH_MASK = (1 << (EDGE_BITS - BUCKET_BITS*2)) - 1;
        for (auto it = nonces_step1->begin(); it != nonces_step1->end(); ++it) {
            BitStreamer<EDGE_BITS-BUCKET_BITS> & in = (*it)[task_idx];
            if (in.get_size()==0)
                continue;

            in.reset_read();
            uint32_t sz = in.get_size();
            for (uint32_t i = 0; i < sz; ++i) {
                uint32_t hash = in.decode_number();

                uint32_t bucket_idx = hash >> (EDGE_BITS - BUCKET_BITS * 2);
                assert(bucket_idx < BUCKETS_NUM);
                nonce_by_bucket[bucket_idx].encode_number(hash & HASH_MASK);
            }
        }

        const uint32_t base_bucket_idx = task_idx * BUCKETS_NUM;

        // Now let's process the buckets one by one
        for (uint32_t i = 0; i < BUCKETS_NUM; ++i) {
            if (nonce_by_bucket[i].get_size() == 0) {
                continue;
            }
            BitStreamer<EDGE_BITS-BUCKET_BITS*2> & nonces = nonce_by_bucket[i];
            nonces.write_flush();
            nonces.reset_read();

            const uint32_t bucket_idx = base_bucket_idx + i;
            assert(BUCKET_SIZE%64==0);
            uint64_t * mask = &(row->get_mask()[bucket_idx * BUCKET_SIZE/64*MASK_BITS]);
            Bucket<EDGE_BITS,BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> & bucket = row->get_buckets()[bucket_idx];

            uint32_t sz = nonces.get_size();
            for (uint32_t n = 0; n < sz; ++n) {
                uint32_t hash = nonces.decode_number();

                assert(hash < BUCKET_SIZE);

                // Updating mask
                const size_t bit_pos = hash * MASK_BITS;
                const size_t ll_pos = bit_pos / 64;
                const size_t offset = bit_pos % 64;
                const uint64_t edge_m = (edge_mask<<offset);
                if ( (mask[ll_pos] & edge_m) != 0 ) {
                    mask[ll_pos] -= uint64_t(1) << offset;

                    if ( (mask[ll_pos] & edge_m) == 0 ) {
                        // something is collapsed, so we can propagate trimming for ANOTHER section
                        const size_t pair_offset = hash % 2 == 0 ? offset+MASK_BITS : offset-MASK_BITS;
                        assert(pair_offset >=0 && pair_offset <= 64-MASK_BITS);
                        assert( (mask[ll_pos] & (edge_mask<<pair_offset)) != 0 );
                        mask[ll_pos] &= ~(edge_mask<<pair_offset);

                        uint nsz = bucket.get_nonces( hash ^ 1, trim_nonces->get_nonce_buf( 1<<MASK_BITS ) );
                        trim_nonces->update_written_num(nsz);
                    }
                }
            }
        }
    }

    trim_nonces->finish();
}


// u/v - both rows. u/vbufs - buffers with nonces that needs to be converted into the hashes.
template<uint8_t EDGE_BITS, uint ELEMENT_SIZE, uint8_t BUCKET_BITS, uint8_t MASK_BITS, uint32_t
        BIT_PACKER_INDEX_INTERVAL>
uint32_t trim(MetalOps &metal,
          RowBuilder<EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> *row,
          std::vector<MTL::Buffer *> *trimmed,
          uint32_t trim_nonces_per_bucket,
          const uint64_t v[4], uint nonceN,
          int threads_num) {
    assert(threads_num>0);

    metal.nonces_to_hash(*trimmed,v,nonceN);

    const uint BUCKETS_NUM = 1 << BUCKET_BITS;
    std::vector<std::vector<BitStreamer<EDGE_BITS-BUCKET_BITS>>> nonces_step1(threads_num);

    std::vector<std::thread> threads;
    std::atomic_int32_t tasks1_requested(0);
    for (int i = 0; i < threads_num; ++i) {
        nonces_step1[i].resize(BUCKETS_NUM);
        threads.emplace_back(bucket_hashes1<EDGE_BITS, BUCKET_BITS>, trimmed, &nonces_step1[i], &tasks1_requested);
    }

    // Wait for all threads to finish
    for (auto &thread: threads) {
        thread.join();
    }
    // buffers are expected to be released into the threads
    trimmed->clear();

    // Second step - process to the buckets level. All in a single task
    threads.clear();
    std::atomic_int32_t tasks2_requested(0);
    std::vector< TrimNonces > trimmed_nonces(threads_num);
    for (int i = 0; i < threads_num; ++i) {
        trimmed_nonces[i].init(&metal, trim_nonces_per_bucket);
        threads.emplace_back(bucket_hashes2<EDGE_BITS,ELEMENT_SIZE, BUCKET_BITS, MASK_BITS,BIT_PACKER_INDEX_INTERVAL>, &nonces_step1, row,
                    &trimmed_nonces[i], &tasks2_requested);
    }

    // Wait for all threads to finish
    for (auto &thread: threads) {
        thread.join();
    }
    threads.clear();

    // Let's check how many was trimmed this time
    uint32_t trimmed_total = 0;
    for (auto & tn : trimmed_nonces) {
        trimmed->insert(trimmed->end(), tn.get_trimmed_nonces().begin(), tn.get_trimmed_nonces().end());
        trimmed_total += tn.get_total_added_nonces();
        tn.get_trimmed_nonces().clear();
    }
    trimmed_nonces.clear();

    return trimmed_total;
}

#endif //TRIMMER_H
