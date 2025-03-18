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

#include "../miner/inplace_bucket_sorting.h"
#include "inplace_bucket_sorting_test.h"

template<int SIZE, int ELEMENT_SIZE, int HASH_BITS, int BUCKET_BITS>
void run_inplace_sorting1(HashNonce hn_data[SIZE]) {
   // Let's sumulate the horizontal buckets layout (pass1)

    // Normal case, perfect distribution

    uint64_t data[SIZE*ELEMENT_SIZE/8];

    const int BUCKET_NUM = 1<<BUCKET_BITS;

    std::vector<HashNonce> buckets[BUCKET_NUM];

    uint32_t hash_mask = (uint64_t(1) << HASH_BITS) - 1;
    std::cout << "test_inplace_sorting1" << std::endl;

    for (int i = 0; i < SIZE; i++) {
        HashNonce hn = hn_data[i];
        hn.hash &= hash_mask;
        std::cout << "In: " << std::hex << hn.nonce << " Hash: " << hn.hash << std::endl;
        store_hash_nonce<ELEMENT_SIZE,HASH_BITS>(data, i, HashNonce(hn.hash, 1) );

        uint bucket_idx = hn.hash >> (HASH_BITS-BUCKET_BITS);
        buckets[bucket_idx].push_back( HashNonce(hn.hash & ((1 << (HASH_BITS-BUCKET_BITS)) - 1), hn.nonce) );
    }

    uint64_t * bucket_addr[BUCKET_NUM];
    BucketIndexLastNonce  bucket_state[BUCKET_NUM];
    uint cur_nonce[BUCKET_NUM];
    for (int i = 0; i < BUCKET_NUM; i++) {
        assert((SIZE*ELEMENT_SIZE * i / BUCKET_NUM) % 8 ==0);
        auto idx = (SIZE*ELEMENT_SIZE * i / BUCKET_NUM) / 8;
        bucket_addr[i] = &data[idx];
        cur_nonce[i] = 0;
        bucket_state[i].init( SIZE * i / BUCKET_NUM - 1, 0);
    }

    // Try to read and check
    uint bucket_sz = SIZE/BUCKET_NUM;
    for (int i = 0; i < SIZE; i++) {
        uint bucket_num = i / bucket_sz;
        uint bucket_idx = i % bucket_sz;

        HashNonce hn = extract_hash_nonce<ELEMENT_SIZE,HASH_BITS>( bucket_addr[bucket_num], bucket_idx );
        assert(hn.hash == hn_data[i].hash);
        assert(hn.nonce == 1);
    }

    std::vector<HashNonce> overflow_buffer;

    uint8_t * backup_buffer = (uint8_t * )std::aligned_alloc(64, SIZE/BUCKET_NUM * ELEMENT_SIZE);

    inplace_bucket_sort<BUCKET_BITS, ELEMENT_SIZE, HASH_BITS>( bucket_addr,
        bucket_state,
        SIZE/BUCKET_NUM,
        0,
        &overflow_buffer,
        backup_buffer);

    std::free(backup_buffer);

    // Now let's validate the results for all buckets
    const uint bucket_len = SIZE/4;
    for (int i = 0; i < SIZE; i++) {
        uint bucket_idx = i / bucket_len;

        HashNonce hn = extract_hash_nonce<ELEMENT_SIZE, HASH_BITS-BUCKET_BITS>(data, i);
        if (hn.is_zero() && i>0) {
            std::cout << "OUT: " << i << " NONE" << std::endl;
            continue;
        }
        std::cout << "In: " << i << "   Hash: " << hn.hash << "  Nonce: " << hn.nonce <<  std::endl;

        cur_nonce[bucket_idx] = cur_nonce[bucket_idx] + hn.nonce;
        hn.nonce = cur_nonce[bucket_idx];

        bool found = false;
        std::vector<HashNonce> & bct = buckets[bucket_idx];
        for (int j=bct.size()-1; j >= 0; j--) {
            if (bct[j].nonce == hn.nonce) {
                found = true;
                bct.erase(bct.begin()+j);
                break;
            }
        }

        assert(found);
    }

    uint expected_overflow_size = 0;
    for (int i = 0; i < 4; i++) {
        expected_overflow_size += buckets[i].size();
    }
    assert(overflow_buffer.size() == expected_overflow_size);
}



void test_inplace_sorting1() {
   // Let's sumulate the horizontal buckets layout (pass1)

    // Normal case, perfect distribution
    const int SIZE = 64;
    const int HASH_BITS = 16;
    uint32_t hash_mask = (1 << HASH_BITS) - 1;

    HashNonce hn_data[SIZE];

    for (int i = 0; i < SIZE; i++) {
        uint32_t hash = (i*27 + 13) & hash_mask;
        hn_data[i] = HashNonce(hash, i);
    }

    run_inplace_sorting1<SIZE, 4, HASH_BITS, 2>(hn_data);
    run_inplace_sorting1<SIZE, 4, 24, 2>(hn_data);
    run_inplace_sorting1<SIZE, 5, HASH_BITS, 2>(hn_data);
    run_inplace_sorting1<SIZE, 5, 31, 2>(hn_data);
    run_inplace_sorting1<SIZE, 5, 32, 2>(hn_data);
}

void test_inplace_sorting2() {
    // Sorting assymetric data, so the buffer oveflow will be used a lot

    const int SIZE = 32;
    const int HASH_BITS = 16;
    uint32_t hash_mask = (1 << HASH_BITS) - 1;

    HashNonce hn_data[SIZE];

    for (int i = 0; i < SIZE/2; i++) {
        uint32_t hash = (i*27 + 13) & hash_mask;
        hn_data[i] = HashNonce(hash, i);
    }
    for (int i = SIZE/2; i<SIZE; i++) {
        uint32_t hash = 13 & hash_mask;
        hn_data[i] = HashNonce(hash, i);
    }
    run_inplace_sorting1<SIZE, 4, HASH_BITS, 2>(hn_data);
    run_inplace_sorting1<SIZE, 4, 24, 2>(hn_data);
    run_inplace_sorting1<SIZE, 5, HASH_BITS, 2>(hn_data);
    run_inplace_sorting1<SIZE, 5, 31, 2>(hn_data);
    run_inplace_sorting1<SIZE, 5, 32, 2>(hn_data);
}

// We have some issue here, related to the backup buffer. Need special test to find the bug
void test_bucket_read_write() {
    const size_t ELEMENT_SIZE = 5;
    const uint8_t HASH_BITS = 10;

    const uint BUCKET_SIZE = 128;

    for (int threshold = 100; threshold < 990; threshold++) {
        uint8_t bucket_data[ELEMENT_SIZE*BUCKET_SIZE];
        // Let's write some data into the buffer
        for (int i = 0; i < BUCKET_SIZE; i++) {
            HashNonce hn(i*2, 1);
            store_hash_nonce<ELEMENT_SIZE,HASH_BITS>((uint64_t *)bucket_data, i, hn);
        }
        // Check by reading
        for (int i = 0; i < BUCKET_SIZE; i++) {
            HashNonce hn = extract_hash_nonce<ELEMENT_SIZE,HASH_BITS>((uint64_t *)bucket_data, i);
            assert(hn.nonce == 1);
            assert(hn.hash == i*2);
        }


        uint8_t backup_data[ELEMENT_SIZE*BUCKET_SIZE];
        uint32_t rand = 2346757234;

        const uint32_t nonceR = 10;
        const uint32_t nonceW = 20;
        BucketIndexLastNonce state;
        state.init(nonceR,nonceW);

        while (state.indexW < BUCKET_SIZE) {
            rand = rand + rand * 13 + 7;
            if (rand % 1000 < threshold && state.indexR<BUCKET_SIZE ) {
                // read
                uint i = state.indexR;
                HashNonce rhn = read_hash_nonce<ELEMENT_SIZE, HASH_BITS>((uint64_t *)bucket_data, state, backup_data );
                assert(rhn.nonce == nonceR + i + 1);
                assert(rhn.hash == i*2);
            }
            else {
                // write
                uint i = state.indexW;
                HashNonce  hn(i*2+1, state.nonceW + i );
                write_hash_nonce<ELEMENT_SIZE, HASH_BITS-2>((uint64_t *)bucket_data, state, hn, backup_data);
            }
        }

        while (state.indexR < BUCKET_SIZE) {
            uint i = state.indexR;
            HashNonce rhn = read_hash_nonce<ELEMENT_SIZE, HASH_BITS>((uint64_t *)bucket_data, state, backup_data );
            assert(rhn.nonce == nonceR + i + 1);
            assert(rhn.hash == i*2);
        }


        for (  int i = 0; i < 64; i++) {
        }

        // Check what was written
        for (int i = 0; i < BUCKET_SIZE; i++) {
            HashNonce hn = extract_hash_nonce<ELEMENT_SIZE,HASH_BITS-2>((uint64_t *)bucket_data, i);
            assert(hn.nonce == i);
            assert(hn.hash == i*2+1);
        }
    }

}

