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

#define ENABLE_TRACE

#include "../miner/row_builder.h"
#include "row_builder_test.h"

#include "../miner/metal.h"
#include "../miner/sip_hash.h"
#include <set>

static uint64_t v[4] = { 0x736f6d6570736575ULL, 0x646f72616e646f6dULL,0x6c7967656e657261ULL,0x7465646279746573ULL};

bool enable_data_dumping = false;

void switch_enable_data_dumping(bool enable) {
    enable_data_dumping = enable;
}


template <uint8_t EDGE_BITS, uint ELEMENT_SIZE, uint8_t BUCKET_BITS>
void test_row_builder_impl(bool do_dump) {
    enable_data_dumping = do_dump;

    // Simple test case with minimal number of elements, buckets num 4x4
    const uint8_t MASK_BITS = 4;
    const uint32_t BIT_PACKER_INDEX_INTERVAL = 32;
    const uint nonceN = 0;

    // Calculating the row data (mask and nonces) manually non optimal way
    std::vector<HashNonce> row_data;
    uint64_t row_edge_num = uint64_t(1) << EDGE_BITS;
    uint edge_mask = (uint)((uint64_t(1) << EDGE_BITS) - 1);
#ifdef ENABLE_TRACE
    if (enable_data_dumping) {
        std::cout << "The hash seed values: " << std::endl;
    }
#endif

    // Limit for C31/C32 test cases. I don't have enough RAM to test on my Mac
   // const uint max_hashes = /*EDGE_BITS>31 ? 0x1ffff :*/ 0x1fffffff;

    for (uint64_t i = 0; i < row_edge_num; i++) {
        uint64_t nonce = i*2 + nonceN;
        uint hash = sip_hash(v, nonce) & edge_mask;
        //if (hash<max_hashes) {
            row_data.push_back( HashNonce(hash, i) );
#ifdef ENABLE_TRACE
            if (enable_data_dumping) {
                if (i%16==0) {
                    std::cout << std::hex << i << "   ";
                }
                std::cout << std::hex << i << ";" << hash << "  ";
                if (i%16==15)
                    std::cout << std::endl;
            }
#endif
        //}
    }
#ifdef ENABLE_TRACE
    if (enable_data_dumping) {
        std::cout << std::endl;
    }
#endif

    std::sort(row_data.begin(), row_data.end(), [](const HashNonce& a, const HashNonce& b) {
        return a.hash < b.hash;
    });

    MetalOps metal(EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS);

    std::vector<MTL::Buffer *> trimmed_nonces;
    RowBuilder<EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> row(metal, v, nonceN, 1, trimmed_nonces);

    //std::vector< std::set<uint32_t> > nonces2hash( std::min(uint64_t(max_hashes), uint64_t(edge_mask)+1));
    std::vector< std::set<uint32_t> > nonces2hash( uint64_t(edge_mask)+1);
    for (const auto & hn : row_data ) {
        assert(hn.hash<=edge_mask);
        //if (hn.hash<max_hashes)
            nonces2hash[hn.hash].insert(hn.nonce);
    }

    // Validating the mask
    const std::vector<uint64_t> & mask = row.get_mask();
    assert(MASK_BITS==4);
    const uint MM = (1 << MASK_BITS) - 1;


    const uint BUCKET_MASK = (1 << (BUCKET_BITS*2)) - 1;
    const std::vector<Bucket<EDGE_BITS,BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL>> & buckets = row.get_buckets();
    assert(buckets.size() == BUCKET_MASK+1);
    assert(row_edge_num % buckets.size() == 0);
    const uint BUCKET_LEN = row_edge_num / buckets.size();

    // 16 items per uint64_t mask
    uint edge_idx = 0;
    IntervalsIndex bucket_read_index;
    uint last_bucket = buckets.size()+1;
    uint32_t nonces_buffer[1<<MASK_BITS];

    size_t nonce_num = uint64_t(1) << EDGE_BITS;
    std::vector<bool> nonces_set( nonce_num, false );
#ifdef ENABLE_TRACE
    if (enable_data_dumping) {
        std::cout << "Trimmed nonces:" << std::endl;
    }
    uint trimmed_counter = 0;
#endif

    for (auto tn : trimmed_nonces) {
        uint l = tn->length() / sizeof(uint32_t);
        uint32_t * data = (uint32_t *) tn->contents();
        for (uint i = 0; i < l; i++) {
            if (data[i] == UINT_MAX) {
                break;
            }
#ifdef ENABLE_TRACE
            if (enable_data_dumping) {
                std::cout << data[i] << "  ";
                if (++trimmed_counter % 16 == 0)
                    std::cout << std::endl;
            }
#endif
            nonces_set[data[i]] = true;
        }
        tn->release();
    }
    trimmed_nonces.clear();


    for (uint64_t m : mask) {
        for (uint i = 0; i < 16; i++) {
            /*if (edge_idx>=max_hashes) {
                break;
            }*/
            uint expected_nonce_num = m & MM;
            m >>= MASK_BITS;
            const std::set<uint32_t> & nonces1 = nonces2hash[edge_idx];
            const std::set<uint32_t> & nonces2 = nonces2hash[edge_idx ^ 1];

            uint bucket_idx = edge_idx / BUCKET_LEN;
            uint bucket_local_idx = edge_idx % BUCKET_LEN;

            if (bucket_idx!=last_bucket) {
                last_bucket = bucket_idx;
                bucket_read_index.reset();
            }

            const Bucket<EDGE_BITS,BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> & bkt = buckets[bucket_idx];

            int sz = bkt.get_nonces(bucket_local_idx, bucket_read_index, nonces_buffer);
            assert(sz==expected_nonce_num);

            if (nonces1.empty() || nonces2.empty()) {
                assert(expected_nonce_num == 0);
                // supposed to be trimmed
                for (auto n : nonces1) {
                    // checking if we have such nonces into the trimmed
                    assert(nonces_set[n]);
                    nonces_set[n] = false;
                }
            }
            else {
                assert( expected_nonce_num == nonces1.size() );
                // Let's read the nonces from the raw

                for (uint q = 0; q < expected_nonce_num; q++) {
                    assert(nonces1.count(nonces_buffer[q]) == 1);
                }
            }

            edge_idx++;
        }
    }

    // check if all trimmed nonces are gone
    for (int r=0; r<nonce_num; r++) {
        assert(nonces_set[r]==false);
    }
}

void test_row_builder_4b() {
    test_row_builder_impl<12, 4, 2>(true);
}

void test_row_builder_5b() {
    test_row_builder_impl<12, 5, 2>(false);
}

void test_row_builder_C25() {
    test_row_builder_impl<20, 5, 1>(false);
    test_row_builder_impl<25, 5, 8>(false);
}

void test_row_builder_C27() {
    test_row_builder_impl<27, 5, 8>(false);
}

static bool aa=false;

void test_row_builder_C28() {
    test_row_builder_impl<28, 5, 8>(false);

    if (aa) {
        // for compiler to check types in templates
        test_row_builder_impl<32, 5, 9>(false);
    }
}

