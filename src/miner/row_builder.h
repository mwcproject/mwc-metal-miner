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

#ifndef ROW_BUILDER_H
#define ROW_BUILDER_H

#include <bitset>
#include <thread>
#include "bucket.h"
#include "inplace_bucket_sorting.h"
#include "metal.h"
#include "trim_nonces.h"

#ifdef ENABLE_TRACE

extern bool enable_data_dumping;


template <uint ELEMENT_SIZE, uint8_t HASH_BITS>
static void dump_buffers(const std::vector<MTL::Buffer*> & buffers,
            uint32_t nonce_base0,
            const std::vector<std::vector<HashNonce>> & overflow_buffers) {

    if (!enable_data_dumping)
        return;

    uint32_t nonce_base;
    for (uint i = 0; i < buffers.size(); i++) {
        MTL::Buffer* buf = buffers[i];
        int el_num = buf->length()/ELEMENT_SIZE;
        uint bucket_size = el_num / buffers.size();
        std::cout << "Buffer " << i << " size " << el_num << std::endl;
        uint64_t * data = (uint64_t *)buf->contents();
        for (uint j=0; j<el_num; j++) {
            if (j % bucket_size == 0)
                nonce_base = nonce_base0;

            HashNonce hn = extract_hash_nonce<ELEMENT_SIZE, HASH_BITS>(data, j);
            if (j%16==0) {
                std::cout << std::hex << j << "   ";
            }
            nonce_base = nonce_base + hn.nonce;
            std::cout << std::hex << hn.nonce << ";" << nonce_base << ";" << hn.hash << "  ";
            if (j%16==15)
                std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    for (uint i=0; i<overflow_buffers.size(); i++) {
        const std::vector<HashNonce> & ov_buf = overflow_buffers[i];
        std::cout << "Overflow Buffer " << i << " size " << ov_buf.size() << std::endl;
        for (uint j=0; j<ov_buf.size(); j++) {
            const HashNonce & hn = ov_buf[j];
            if (j%16==0) {
                std::cout << std::hex << j << "   ";
            }
            std::cout << std::hex << hn.nonce << ";" << hn.hash << "  ";
            if (j%16==15)
                std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}
#endif

template <uint8_t EDGE_BITS, uint ELEMENT_SIZE, uint8_t BUCKET_BITS, uint8_t MASK_BITS, uint32_t BIT_PACKER_INDEX_INTERVAL>
class RowBuilder {
private:
    // Mask will be 4 bits per edge, it is up to 15 edges can be collided (others will be abandoned)
    std::vector<uint64_t> mask;
    // 256x256 buckets data
    std::vector<Bucket<EDGE_BITS,BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL>> buckets;

public:
    // Output: resulting_trimmed_nonces. Caller should release those buffers
    RowBuilder(MetalOps & metal, uint64_t hash_v[4], uint nonceN, int threads_num,
        std::vector<MTL::Buffer *> & resulting_trimmed_nonces) {

        const uint32_t BUCKETS_NUM = 1 << BUCKET_BITS;
        std::vector<MTL::Buffer*> buffers = metal.allocate_buffers();
        assert(buffers.size() == BUCKETS_NUM);

        metal.seed_hash(buffers, hash_v, nonceN );

        // Now we can do bucket sorting...
        if (threads_num<=0) {
            threads_num = std::thread::hardware_concurrency();
        }

        assert(threads_num>0);
        std::vector<std::vector<HashNonce>> overflow_buffers(threads_num);

#ifdef ENABLE_TRACE
        if (enable_data_dumping)
        {
            // Dumping the buffers with seed
            std::cout << "Seed data: " << (uint64_t(1)<<EDGE_BITS) << std::endl;
            dump_buffers<ELEMENT_SIZE, EDGE_BITS>( buffers, -1, overflow_buffers);
        }
#endif

        std::vector<std::thread> threads;
        std::atomic_int32_t tasks1_requested(0);
        for (int i = 0; i < threads_num; ++i) {
            threads.emplace_back(sorting_data1, &overflow_buffers[i], &tasks1_requested, &buffers );
        }

        {
            // While we are waiting, we can prepare for the mask construction (step 3)
            uint32_t bucket_total = BUCKETS_NUM * BUCKETS_NUM; // Total number of buckets (2 passes total)
            // Mask will be 4 bits per edge, it is up to 15 edges can be collided (others will be abandoned)
            uint64_t elements_num = (uint64_t(1) << EDGE_BITS);
            assert(elements_num%64==0);
            mask.resize( elements_num*MASK_BITS/64, 0 );
            buckets.resize(bucket_total);
        }

        // Wait for all threads to finish
        for (auto& thread : threads) {
            thread.join();
        }

#ifdef ENABLE_TRACE
        if (enable_data_dumping)
        {
            // Dumping the buffers with seed
            std::cout << "Step1 data: " << (uint64_t(1)<<EDGE_BITS) << std::endl;
            dump_buffers<ELEMENT_SIZE, EDGE_BITS-BUCKET_BITS>(buffers, 0, overflow_buffers);
        }
#endif

        // Initiating the second round...
        threads.resize(0);
        std::atomic_int32_t tasks2_requested(0);
        for (int i = 0; i < threads_num; ++i) {
            threads.emplace_back(sorting_data2, &overflow_buffers[i], &tasks2_requested, &buffers );
        }
        // Wait for all threads to finish because we still need to process the overheads.
        for (auto& thread : threads) {
            thread.join();
        }

#ifdef ENABLE_TRACE
        // Dumping the buffers with seed
        if (enable_data_dumping)
        {
            std::cout << "Step2 data: " << (uint64_t(1)<<EDGE_BITS) << std::endl;
            dump_buffers<ELEMENT_SIZE, EDGE_BITS-BUCKET_BITS*2>(buffers, 0, overflow_buffers);
        }
#endif

        // At this point need to process overhead. Buckets without overhead can be processes earlier, BUT we want to manage memory
        // Buffers needs to be released asap. That is why we better here split overhead into the 256 buckets. So we can start other threads asap.
        // Also, BUCKETS_NUM is CPU cache friendly number.

        std::vector<std::vector<HashNonce>> overflow_by_columnt(BUCKETS_NUM);
        const uint READ_HASH_TO_BUCKET_SHIFT = EDGE_BITS-BUCKET_BITS;

        for (const auto & ov_buf : overflow_buffers) {
            for (const auto & ov : ov_buf) {
                uint32_t idx = ov.hash >> READ_HASH_TO_BUCKET_SHIFT;
                assert(idx<BUCKETS_NUM);
                overflow_by_columnt[idx].push_back(ov);
            }
        }
        // Now let's process columnts in different threads. So every thread will process BUCKETS_NUM buckets in a single task.
        // It is easier to manage the memory, we don't want overcomplicate with a memory management. This way memory overhead
        // is num_cpu/256  - it is low number for most platforms. Miner is done for CPU with 10-20 cores.

        // Building the mask from the buckets data. Bucket size is about 4096 elements for C31, so all data should fit the L1 cache well
        threads.resize(0);
        std::atomic_int32_t masks_build_requests(0);
        std::vector<TrimNonces> trim_nonces(threads_num);
        for (int i = 0; i < threads_num; ++i) {
            trim_nonces[i].init(&metal, 1<<(EDGE_BITS-MASK_BITS) );
            threads.emplace_back(build_mask, &masks_build_requests, &buffers, &overflow_by_columnt,
                &mask, &buckets, &trim_nonces[i]);
        }
        // Wait for all threads to finish because we still need to process the overheads.
        for (auto& thread : threads) {
            thread.join();
        }

        assert(resulting_trimmed_nonces.empty());
        for (auto& tn : trim_nonces) {
            resulting_trimmed_nonces.insert(resulting_trimmed_nonces.begin(),
                tn.get_trimmed_nonces().begin(), tn.get_trimmed_nonces().end());
            tn.get_trimmed_nonces().clear();
        }

#ifdef ENABLE_TRACE
        if (enable_data_dumping) {
            std::cout << "Resulting mask, size: " << mask.size() << std::endl;
            for (uint j=0; j<mask.size(); j++) {
                uint64_t m = mask[j];
                if (j%2==0) {
                    std::cout << std::hex << (j*16) << "   ";
                }
                std::cout << std::hex << m << " ";
                if (j%2==1)
                    std::cout << std::endl;
            }
            std::cout << std::endl;
        }
#endif
    }

    std::vector<uint64_t> & get_mask() { return mask; }
    std::vector<Bucket<EDGE_BITS,BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL>> & get_buckets() {return buckets;}

    bool has_mask(uint32_t mask_idx) const {
        assert( mask_idx < mask.size() * (64/MASK_BITS) );
        const size_t bit_pos = mask_idx * MASK_BITS;
        const size_t ll_pos = bit_pos / 64;
        const size_t offset = bit_pos % 64;
        const uint64_t edge_m = ( uint64_t((1 << MASK_BITS)-1) << offset );
        const uint64_t m = mask[ll_pos];
        return (m & edge_m) != 0;
    }
private:
    // Initial sortign by row
    static void sorting_data1(std::vector<HashNonce> * overflow_buffer, std::atomic_int32_t * requested_task,
            std::vector<MTL::Buffer*> * buffers) {
        const uint64_t elements_num = uint64_t(1) << EDGE_BITS;

        const uint32_t BUCKETS_NUM = 1 << BUCKET_BITS;

        // Buckets
        uint64_t * bucket_addr[BUCKETS_NUM];
        // States of the buckets
        BucketIndexLastNonce  bucket_state[BUCKETS_NUM];
        const uint32_t bucket_size = elements_num / BUCKETS_NUM / BUCKETS_NUM;

        uint8_t * backup_buffer = (uint8_t * )std::aligned_alloc(64, elements_num * ELEMENT_SIZE);

        while (true) {
            int32_t task_idx = requested_task->fetch_add(1, std::memory_order_relaxed);
            if (task_idx >= BUCKETS_NUM)
                break;

            uint32_t baseOffset = task_idx * bucket_size * ELEMENT_SIZE;

            // Note, for first it is expected unit 0 - 1 => FFFF... and it will be back to 0 on first nonce calculation
            uint32_t row_base_nonce = uint32_t(task_idx * bucket_size * BUCKETS_NUM);
            uint32_t base_nonce = row_base_nonce - 1;//BUCKETS_NUM;

            for (uint i = 0; i < BUCKETS_NUM; ++i) {
                uint8_t * baseAddr = (uint8_t *)(*buffers)[i]->contents(); // Buffers are vertically orienter
                bucket_addr[i] = (uint64_t *)( baseAddr + baseOffset );
                bucket_state[i].init( base_nonce, row_base_nonce);
                base_nonce += bucket_size; // ++;
            }

            inplace_bucket_sort<BUCKET_BITS, ELEMENT_SIZE, EDGE_BITS>( bucket_addr,
                    // States of the buckets
                    bucket_state,
                    bucket_size,
                    0,
                    overflow_buffer,
                    backup_buffer);
        }

        std::free(backup_buffer);
    }

    // Secondary sorting by columnts
    static void sorting_data2(std::vector<HashNonce> * overflow_buffer, std::atomic_int32_t * requested_task,
            std::vector<MTL::Buffer*> * buffers) {

        const uint64_t elements_num = uint64_t(1) << EDGE_BITS;

        const uint32_t BUCKETS_NUM = 1 << BUCKET_BITS;

        // Buckets
        uint64_t * bucket_addr[BUCKETS_NUM];
        // States of the buckets
        BucketIndexLastNonce  bucket_state[BUCKETS_NUM];
        const uint32_t bucket_size = elements_num / BUCKETS_NUM / BUCKETS_NUM;
        uint8_t * backup_buffer = (uint8_t * )std::aligned_alloc(64, elements_num * ELEMENT_SIZE);

        while (true) {
            int32_t task_idx = requested_task->fetch_add(1, std::memory_order_relaxed);
            if (task_idx >= BUCKETS_NUM)
                break;

            uint8_t * column_buffer_base = (uint8_t *)(*buffers)[task_idx]->contents();

            for (uint i = 0; i < BUCKETS_NUM; ++i) {
                bucket_addr[i] = (uint64_t *)( column_buffer_base + bucket_size*ELEMENT_SIZE*i );
                uint32_t base_nonce = i * bucket_size * BUCKETS_NUM;
                bucket_state[i].init( base_nonce, 0 );
            }

            uint32_t hash_base = uint32_t(task_idx) << (EDGE_BITS-BUCKET_BITS);

            inplace_bucket_sort<BUCKET_BITS, ELEMENT_SIZE, EDGE_BITS - BUCKET_BITS>( bucket_addr,
                    // States of the buckets
                    bucket_state,
                    bucket_size,
                    hash_base,
                    overflow_buffer,
                    backup_buffer);
        }
        std::free(backup_buffer);
    }

    // Sorting for buckets. Single task - process the whole column. That should take
    static void build_mask(std::atomic_int32_t * requested_task,
            std::vector<MTL::Buffer*> * buffers,
            std::vector<std::vector<HashNonce>> * overflow_by_column,
            std::vector<uint64_t> * mask,
            std::vector<Bucket<EDGE_BITS,BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL>> * buckets,
            TrimNonces * trim_nonces) {

        const uint32_t BUCKETS_NUM = 1 << BUCKET_BITS;
        const uint64_t elements_num = uint64_t(1) << EDGE_BITS;
        const uint32_t bucket_size = elements_num / BUCKETS_NUM / BUCKETS_NUM;

        const uint READ_HASH_TO_BUCKET_SHIFT = EDGE_BITS-(BUCKET_BITS*2);
        const uint BUCKET_MASK = (1 << BUCKET_BITS) - 1;

        while (true) {
            int32_t task_idx = requested_task->fetch_add(1, std::memory_order_relaxed);
            if (task_idx >= BUCKETS_NUM)
                break;

            // Let's process overhead first
            std::vector<std::vector<HashNonce>> overflow_by_bucket(BUCKETS_NUM);
            for (const auto & ov : (*overflow_by_column)[task_idx]) {
                uint32_t idx = (ov.hash >> READ_HASH_TO_BUCKET_SHIFT) & BUCKET_MASK;
                overflow_by_bucket[idx].push_back(ov);
            }

            // we can release memory for overflows
            (*overflow_by_column)[task_idx].clear();

            uint8_t * column_buffer_base = (uint8_t *)(*buffers)[task_idx]->contents();

            assert(bucket_size % 64 == 0);

            // Processing buckets one by one
            for (uint i = 0; i < BUCKETS_NUM; ++i) {

                const uint bucket_idx = task_idx * BUCKETS_NUM + i;
                assert(bucket_size%64==0);
                const uint maskIdx = bucket_idx * bucket_size/64*MASK_BITS;

                process_bucket( (uint64_t *)( column_buffer_base + bucket_size*ELEMENT_SIZE*i),
                                overflow_by_bucket[i],
                                &(*mask)[maskIdx],
#ifndef NDEBUG
                                bucket_idx * bucket_size,
#endif
                                &(*buckets)[bucket_idx],
                                trim_nonces);
            }
            // Releasing the buffer memory
            (*buffers)[task_idx]->release();
        }

        trim_nonces->finish();
    }

    // start nonce is 0, not requesting
    static void process_bucket(uint64_t * bucket_data, std::vector<HashNonce> & overhead_data,
            uint64_t * mask,
#ifndef NDEBUG
            uint mask_base_idx,
#endif
            Bucket<EDGE_BITS, BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> * res_bucket,
            TrimNonces * trim_nonces ) {
        std::vector<HashNonce> src_data;
        const uint64_t elements_num = uint64_t(1) << EDGE_BITS;
        const uint32_t BUCKETS_NUM = 1 << BUCKET_BITS;
        const uint32_t BUCKET_SIZE = elements_num / BUCKETS_NUM / BUCKETS_NUM;
        src_data.reserve(BUCKET_SIZE + BUCKET_SIZE/2);

        const uint64_t edge_mask = (1 << MASK_BITS) - 1;
        const uint64_t double_edge_mask = (1 << MASK_BITS*2) - 1;

        std::bitset<BUCKET_SIZE> banned_items;

        const uint HASH_BITS = EDGE_BITS-BUCKET_BITS*2;
        const uint HASH_LOW_MASK = (1 << HASH_BITS)-1;

#ifndef NDEBUG
        uint added_records = 0;

        // validate that this bucket mask is all zeroes
        for (uint i = 0; i < BUCKET_SIZE*MASK_BITS/64; ++i) {
            uint64_t m = mask[i];
            assert(m==0);
        }
#endif

        for (uint i = 0; i < BUCKET_SIZE; ++i) {
            HashNonce val = extract_hash_nonce<ELEMENT_SIZE,HASH_BITS>( bucket_data, i );
            assert(val.hash < BUCKET_SIZE);
            if (val.is_zero() && i>0)
                break;

#ifndef NDEBUG
   //         if (val.hash + mask_base_idx == 0x10e) {
     //           int rr=0;
       //     }
#endif

            const size_t bit_pos = val.hash * MASK_BITS;
            const size_t ll_pos = bit_pos / 64;
            const size_t offset = bit_pos % 64;
            const uint64_t edge_m = (edge_mask<<offset);
            if ( (mask[ll_pos] & edge_m) != edge_m) {
                mask[ll_pos] += uint64_t(1) << offset;
#ifndef NDEBUG
                added_records++;
#endif
            }
            else {
                banned_items.set(i);
            }
        }

        // processing separately overhead data
        for (int i = overhead_data.size() - 1; i >= 0; --i) {
            auto val = overhead_data[i];
            val.hash &= HASH_LOW_MASK;
            assert(val.hash < BUCKET_SIZE);
            const size_t bit_pos = val.hash * MASK_BITS;
            const size_t ll_pos = bit_pos / 64;
            const size_t offset = bit_pos % 64;
            const uint64_t edge_m = (edge_mask<<offset);
            if ( (mask[ll_pos] & edge_m) != edge_m) {
                mask[ll_pos] += uint64_t(1) << offset;
#ifndef NDEBUG
                added_records++;
#endif
            }
            else {
                // should be very rare case
                overhead_data.erase(overhead_data.begin()+i);
            }
        }

        // The mask is ready. Now we can calculate the total numbers, offsets and then do the final sorting by putting the data into the place
        // Problem that even buckets are small, they can be too dence, so uint16_t will not be enough
        std::vector<int16_t> offsets16(BUCKET_SIZE, 0);
        std::vector<int32_t> offsets32;
        uint32_t curOffs = 0;
        size_t offIdx = 0;
        for (uint i = 0; i < BUCKET_SIZE*MASK_BITS/64; ++i) {
            uint64_t m = mask[i];
            //we need top do trimming here as well. Checking the pairs
            for (uint j = 0; j < 64/MASK_BITS/2; ++j) {
                const uint64_t m1 = edge_mask << (j*2*MASK_BITS);
                const uint64_t m2 = edge_mask << ((j*2+1)*MASK_BITS);
                bool none1 = (m & m1) == 0;
                bool none2 = (m & m2) == 0;

                if (none1 != none2) {
                    // need to trim some that pair, resetting the mask
                    m &= ~(double_edge_mask << j*2*MASK_BITS);
                }
            }
            mask[i] = m;

            for (uint j = 0; j < 64/MASK_BITS; ++j) {
                uint16_t num_items = uint16_t(m & edge_mask);
                m >>= MASK_BITS;

                if (curOffs< 0x7FF0) {
                    offsets16[offIdx] = num_items==0 ? -1 : curOffs;
                }
                else {
                    if (offsets32.empty()) {
                        // switching to 32 bit indexes. It is a rare case, but it is possible
                        offsets32.resize(BUCKET_SIZE, 0);
                        auto it16 = offsets16.begin();
                        for ( auto it32 = offsets32.begin(); it32 != offsets32.end(); ++it32, ++it16 ) {
                            *it32 = *it16;
                        }
                        offsets16.clear();
                    }
                    offsets32[offIdx] = num_items==0 ? -1 : curOffs;
                }

                if (offIdx % BIT_PACKER_INDEX_INTERVAL == 0 && offIdx>0) {
                    assert(offIdx/BIT_PACKER_INDEX_INTERVAL>0);
                    res_bucket->interval_indexes[offIdx/BIT_PACKER_INDEX_INTERVAL-1].set(offIdx, res_bucket->intervals.get_writer_bit_pos(), curOffs);
                }
                res_bucket->intervals.encode_number(num_items);
                ++offIdx;
                curOffs += num_items;
            }
        }
        assert(curOffs<=added_records);
        assert(offIdx==BUCKET_SIZE);
        // res_bucket->intervals & interval_indexes are ready.
        res_bucket->intervals.write_flush();

        int32_t nonceR = 0;

        res_bucket->nonces.resize(curOffs);

        // Offsets are ready, now we can sort nonces by the hash masks
        for (uint i = 0; i < BUCKET_SIZE; ++i) {
            HashNonce val = extract_hash_nonce<ELEMENT_SIZE,HASH_BITS>( bucket_data, i );
            assert(val.hash < BUCKET_SIZE);
            if (val.is_zero() && i>0)
                break;

            nonceR = nonceR + val.nonce;
            val.nonce = nonceR;

            if (banned_items[i])
                continue;

            if (offsets32.empty()) {
                if (offsets16[val.hash]<0) {
                    trim_nonces->add_nonce(val.nonce);
                }
                else {
                    res_bucket->nonces[ offsets16[val.hash]++ ] = val.nonce;
                }
            } else {
                if (offsets32[val.hash]<0) {
                    trim_nonces->add_nonce(val.nonce);
                }
                else {
                    res_bucket->nonces[ offsets32[val.hash]++ ] = val.nonce;
                }
            }
        }


        // processing separately overhead data
        for (int i = overhead_data.size() - 1; i >= 0; --i) {
            auto val = overhead_data[i];
            val.hash &= HASH_LOW_MASK;
            assert(val.hash < BUCKET_SIZE);
            if (offsets32.empty()) {
                if (offsets16[val.hash]<0) {
                    trim_nonces->add_nonce(val.nonce);
                }
                else {
                    res_bucket->nonces[ offsets16[val.hash]++ ] = val.nonce;
                }
            } else {
                if (offsets32[val.hash]<0) {
                    trim_nonces->add_nonce(val.nonce);
                }
                else {
                    res_bucket->nonces[ offsets32[val.hash]++ ] = val.nonce;
                }
            }
        }
        // res_bucket->nonces are ready.
    }
};

#endif //ROW_BUILDER_H
