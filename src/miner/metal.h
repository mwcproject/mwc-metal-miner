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

#ifndef METAL_H
#define METAL_H

#include <Metal/Metal.hpp>

class MetalOps {
    MTL::Device* device;
    MTL::Library* sip_hash_library;

    uint32_t EDGE_BITS;
    uint32_t ELEMENT_SIZE;
    uint32_t BUCKET_BITS;
public:
    MetalOps(uint32_t EDGE_BITS, uint32_t ELEMENT_SIZE, uint32_t BUCKET_BITS);
    ~MetalOps();

    // Build buffers. Caller will own them because we want to process them and release after one by one.
    std::vector<MTL::Buffer*> allocate_buffers();

    // Create a new not managed/tracked buffer. Caller is responsible to release it
    MTL::Buffer* create_buffer(uint buffer_size);

    // Calculate the SEED hashes for all buffers
    void seed_hash(const std::vector<MTL::Buffer*>& buffers,
            const uint64_t v[4],
            uint nonceN);

    // Convert nonces into hashes
    void nonces_to_hash(const std::vector<MTL::Buffer*>& buffers, const uint64_t v[4], uint nonceN);

    // Finds nonces for bunch of U/V pairs
    void uv_to_nonces(const uint64_t v[4], MTL::Buffer* U_hash_table, MTL::Buffer* VV_table, MTL::Buffer* NN_table);
};

#endif //METAL_H
