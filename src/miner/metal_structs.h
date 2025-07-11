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

#ifndef METAL_STRUCTS_H
#define METAL_STRUCTS_H

#include <Metal/Metal.hpp>
#include "features.h"

// Note, those structures must match to what we already have at metal_code.cpp
// The same alignment is required, so don't use anything fancy

struct Step1Params {
    // Hash keys
    uint64_t key0;
    uint64_t key1;
    uint64_t key2;
    uint64_t key3;

    // Bucket size and there addresses
    uint32_t bucket_size;
    uint32_t bucket_blocks[128];
};

struct Step2Params {
    // Hash keys
    uint64_t key0;
    uint64_t key1;
    uint64_t key2;
    uint64_t key3;

    // Current bucket index
    uint32_t row_idx;

    // Sizes on input and output memory blocks
    uint32_t in_block_size;
    uint32_t out_block_size;

    // Addressed to source data, high 4 bits are reserved for the block size
    uint32_t in_blocks[128];

    // Addresses for the buckets data, high 4 bits are reserved for the block size
    uint32_t out_blocks[128];
};

// Trimming data
struct TrimmingParams {
    // Mask SIZE multiplier
    uint32_t mask_scale_k;
    // Mask indexes multiplier
    uint32_t mask_scale_mul;
    // Size of the block
    uint32_t block_size;
    // Type of the U/V pass.  Number of the current bucket
    uint32_t pass_bucket;

    // addresses of the blocks. Must be LAST item because of caller c++ code
    uint32_t blocks[128*128];
};

struct DiscoverParams {
    // Hash keys
    uint64_t key0;
    uint64_t key1;
    uint64_t key2;
    uint64_t key3;

    // half hash hash table length
    uint32_t hh_len;
    // half hash hash table hash functions
    uint32_t hh_k;
    // Resulting size of the nonces buffer
    uint32_t nonce_sz;
    // Resulting nonces position
    uint32_t nonce_pos;

    // Type of hashes. U or V
    bool upass;
};


#ifdef USE_METRICS
struct alignas(16) MetricsData {
    uint32_t skipped_edges;
};
#endif


#endif //METAL_STRUCTS_H
