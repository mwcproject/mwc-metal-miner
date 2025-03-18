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

#include "../miner/trimmer.h"
#include "trimmer_test.h"
#include "../miner/metal.h"

void test_trim_step1() {
    const uint8_t EDGE_BITS = 12;
    const uint8_t BUCKET_BITS = 4;

    MetalOps metal(EDGE_BITS, 5, BUCKET_BITS);

    std::vector<uint32_t> generated_hashes;

    uint32_t nonce_cur = 27 * 31;

    std::vector<MTL::Buffer *> hash_data;
    for (int i = 0; i < 3; i++) {
        int nonces_num = (i+2) * 100;
        MTL::Buffer* buf = metal.create_buffer(nonces_num * sizeof(uint32_t));
        hash_data.push_back(buf);
        uint32_t * data = (uint32_t *) buf->contents();
        for (int q=0; q<nonces_num; q++) {
            if (q<nonces_num-10) {
                nonce_cur = nonce_cur*31 + 7;
                uint32_t hash = nonce_cur & ((1 << EDGE_BITS) - 1);
                if (hash == UINT_MAX)
                    hash = 7;
                generated_hashes.push_back(hash);
                data[q] = hash;
            }
            else {
                data[q] = UINT_MAX;
            }
        }
    }

    std::vector<BitStreamer<EDGE_BITS-BUCKET_BITS>> result(1 << BUCKET_BITS);
    std::atomic_int32_t tasks(0);

    bucket_hashes1<EDGE_BITS, BUCKET_BITS>(&hash_data,
        &result,
        &tasks);

    for (int q=0; q<(1 << BUCKET_BITS); q++) {
        result[q].reset_read();
    }

    // Validating...
    for (uint32_t hash : generated_hashes) {
        uint32_t bucket_idx = hash >> (EDGE_BITS - BUCKET_BITS);
        uint32_t bucket_hash = hash & ((1 << (EDGE_BITS - BUCKET_BITS))-1);
        uint32_t got_hash = result[bucket_idx].decode_number();
        assert(got_hash == bucket_hash);
    }

}
