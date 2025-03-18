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

#ifndef INPLACE_BUKET_SORING_H
#define INPLACE_BUKET_SORING_H

#include <cstdint>
#include "assert.h"
#include "data_storage.h"

struct BucketIndexLastNonce {
    // Indexes are read and write because closer to the end, the rotation circle might be broken, that is why
    // we might read more then write
    uint32_t  indexR;
    uint32_t  indexW;
    uint32_t  backedU8Idx;
    // Because nonces are encoded as increments, read/write heads will be different
    uint32_t  nonceR; // reading
    uint32_t  nonceW; // writing

    inline void init(uint32_t _nonceR, uint32_t _nonceW) {
        indexR = indexW = backedU8Idx = 0;
        nonceR = _nonceR;
        nonceW = _nonceW; // write nonce ALLWAYS start from 0.
    }
};

template <uint32_t ELEMENT_SIZE,  uint8_t HASH_BITS>
static inline HashNonce read_hash_nonce(uint64_t * bucker_addr, BucketIndexLastNonce & bucket_state, uint8_t * backup_buffer ) {
    HashNonce hn;
    if (bucket_state.indexR>=bucket_state.indexW) {
        // reading from the original data
        hn = extract_hash_nonce<ELEMENT_SIZE,HASH_BITS>( bucker_addr , bucket_state.indexR );
    }
    else {
        // reading from the backup (write suppose to make it ready.)
        hn = extract_hash_nonce<ELEMENT_SIZE,HASH_BITS>( (uint64_t *)backup_buffer , bucket_state.indexR );
    }
    if (!hn.is_zero() || bucket_state.indexR==0) {
        hn.nonce += bucket_state.nonceR;
        bucket_state.nonceR = hn.nonce;
        bucket_state.indexR++;
    }
    //std::cout << "read from " << (bucket_state.indexR-1) << " hash: " << hn.hash << " nonce:" << hn.nonce << std::endl;
    return hn;
}

template <uint32_t ELEMENT_SIZE, uint8_t HASH_BITS>
static inline void write_hash_nonce(uint64_t * bucket_addr, BucketIndexLastNonce & bucket_state, HashNonce & hn, uint8_t * backup_buffer) {
    //std::cout << "Saving into " << (bucket_state.indexW) << " hash: " << hn.hash << " nonce:" << hn.nonce << std::endl;

    hn.hash &= (1 << HASH_BITS)-1;
    assert(hn.nonce >= bucket_state.nonceW);
    hn.nonce -= bucket_state.nonceW;
    const uint64_t maxNonce = (uint64_t(1) << (ELEMENT_SIZE*8 - HASH_BITS)) - 1;
    if (hn.nonce > maxNonce) {
        hn.nonce = maxNonce;
        assert(false); // might happens, but want to see how often and if need to debug
    }
    bucket_state.nonceW += hn.nonce;
    // Check if backup is needed
    if (bucket_state.indexW>=bucket_state.indexR && bucket_state.indexW>=bucket_state.backedU8Idx/ELEMENT_SIZE) {
        // backup all. Tmp solution, looking for bug

        // let's backup to the next 64 bytes. all data is supposed to be aligned, so 64 can be applicable to indexes
        uint idx1 = std::max(bucket_state.indexW*ELEMENT_SIZE,bucket_state.backedU8Idx);
        uint idx2 = (bucket_state.indexW + 1) * ELEMENT_SIZE;
        idx2 += 64 - idx2 % 64; // Ending at next 64 byte block (cache line)
        assert(idx2>idx1);
        assert(idx2%64==0);
        bucket_state.backedU8Idx = idx2;
        assert(bucket_state.backedU8Idx/ELEMENT_SIZE > bucket_state.indexW);
        std::memcpy( backup_buffer + idx1, ((uint8_t *)bucket_addr) + idx1, idx2-idx1 );
    }
    store_hash_nonce<ELEMENT_SIZE,HASH_BITS>(bucket_addr, bucket_state.indexW, hn);
    bucket_state.indexW++;
}

template <uint8_t BUCKET_BITS, size_t ELEMENT_SIZE,  uint8_t HASH_BITS>
void inplace_bucket_sort( uint64_t * bucket_addr[1<<BUCKET_BITS],
        // States of the buckets
        BucketIndexLastNonce  bucket_state[1<<BUCKET_BITS],
        const uint32_t bucket_size,
        const uint32_t hash_base,
        std::vector<HashNonce> * overflow_buffer,
        uint8_t * backup_buffer) {
    // we need to get some starting data. That data will bounce back and force until it will finished or
    // some data will go into the overflow_buffer. It is easy to proove.
    // because of  that we will need to mark such values. Let's use 0 counter value for that. 0 can't exist, all counters are unique

    assert(HASH_BITS > BUCKET_BITS);

    const uint BUCKETS_NUM = 1<<BUCKET_BITS;

    const uint READ_HASH_TO_BUCKET_SHIFT = HASH_BITS-BUCKET_BITS;

    // It is sequentual scanning...
    for (uint bucket = 0; bucket < BUCKETS_NUM; bucket++) {
        assert(bucket_state[bucket].indexR==0);
        uint8_t * bucket_backup = backup_buffer + bucket_size * bucket * ELEMENT_SIZE;
        uint64_t * buck_addr = bucket_addr[bucket];
        BucketIndexLastNonce & buck_state = bucket_state[bucket];
        for (uint i = 0; i < bucket_size; i++) {
            assert(buck_state.indexR==i);
            HashNonce element = read_hash_nonce<ELEMENT_SIZE,HASH_BITS>(buck_addr, buck_state, bucket_backup );
            if (element.is_zero() && i>0)
                break;

            uint32_t next_bucket = element.hash >> READ_HASH_TO_BUCKET_SHIFT;
            // check if we can write into the next bucket
            if (bucket_state[next_bucket].indexW < bucket_size) {
                write_hash_nonce<ELEMENT_SIZE, HASH_BITS - BUCKET_BITS>(
                        bucket_addr[next_bucket], bucket_state[next_bucket], element, backup_buffer + bucket_size * next_bucket * ELEMENT_SIZE);
            }
            else {
                // we can write, so we need to stash that data item.
                element.hash += hash_base;
                overflow_buffer->push_back(element);
            }
        }
    }

/* // It is loop implementation. Works great but nonces will be no in the sequence - it is really bad.
      const uint64_t READ_MAX_HASH = uint64_t(1)<<HASH_BITS;
      int cur_seed_bucket = 0;
      while (cur_seed_bucket < BUCKETS_NUM) {
        // generating the seed

        HashNonce cur_element;

        if (bucket_state[cur_seed_bucket].indexR < bucket_size ) {
            // Extract curElement
            cur_element = read_hash_nonce<ELEMENT_SIZE,HASH_BITS>(bucket_addr[cur_seed_bucket], bucket_state[cur_seed_bucket], backup_buffer + bucket_size * cur_seed_bucket * ELEMENT_SIZE );
            if (cur_element.is_zero()) {
                // No data state
                cur_seed_bucket++;
                continue;
            }
            assert(cur_element.hash < READ_MAX_HASH);

            // If belong to this bucket, no need to start the rotation chain reaction
            if (cur_element.hash >> READ_HASH_TO_BUCKET_SHIFT == cur_seed_bucket) {
                assert(bucket_state[cur_seed_bucket].indexR>=bucket_state[cur_seed_bucket].indexW);
                // write should be allways possible
                write_hash_nonce<ELEMENT_SIZE, HASH_BITS - BUCKET_BITS>(
                        bucket_addr[cur_seed_bucket], bucket_state[cur_seed_bucket], cur_element, backup_buffer + bucket_size * cur_seed_bucket * ELEMENT_SIZE);
                continue;
            }
        }
        else {
            cur_seed_bucket++;
            continue;
        }

        // Here we are ready to start the exchange circle
        while (true) {
            uint32_t bucket = cur_element.hash >> READ_HASH_TO_BUCKET_SHIFT;
            HashNonce nextElement;

            if (bucket == cur_seed_bucket) {
                // the circle is ended successfully. We got back at the point where we started from
                // write should be allways possible because we read and started form there
                write_hash_nonce<ELEMENT_SIZE, HASH_BITS - BUCKET_BITS>(
                        bucket_addr[cur_seed_bucket], bucket_state[cur_seed_bucket], cur_element, backup_buffer + bucket_size * cur_seed_bucket * ELEMENT_SIZE);
                break;
            }

            // let's check if we can write.
            bool gotNext = false;
            if (bucket_state[bucket].indexW < bucket_size) {
                // we can write, now need to read first
                if (bucket_state[bucket].indexR < bucket_size) {
                    nextElement = read_hash_nonce<ELEMENT_SIZE,HASH_BITS>(bucket_addr[bucket], bucket_state[bucket], backup_buffer + bucket_size * bucket * ELEMENT_SIZE);
                    gotNext = !nextElement.is_zero(); // hash 0 - mean empty value
                }

                write_hash_nonce<ELEMENT_SIZE, HASH_BITS - BUCKET_BITS>(
                        bucket_addr[bucket], bucket_state[bucket], cur_element, backup_buffer + bucket_size * bucket * ELEMENT_SIZE);
            }
            else {
                // we can write, so we need to stash that data item.
                overflow_buffer->push_back(cur_element);
                // the cycle is finished on stash.
                break;
            }

            if (!gotNext) {
                // circle is finished without stashing. Last item was able to write, but
                // nothing to read. The chain is broken. Just exiting the cirle loop.
                break;
            }

            // Here we are good, we can continue process the chain.
            cur_element = nextElement;
        }

    }*/

    // Now we need to finalize ALL, some data was read but never write after
    for (int i = 0; i < BUCKETS_NUM; i++) {
        HashNonce zero(0,0);
        uint64_t * bkt = bucket_addr[i];
        BucketIndexLastNonce & state = bucket_state[i];
        while(state.indexW < bucket_size) {
            store_hash_nonce<ELEMENT_SIZE,HASH_BITS-BUCKET_BITS>(bkt, state.indexW, zero);
            state.indexW++;
        }
    }

}

#endif //INPLACE_BUKET_SORING_H
