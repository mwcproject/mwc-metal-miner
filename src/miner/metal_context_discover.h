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

#ifndef METAL_CONTEXT_DISCOVER_H
#define METAL_CONTEXT_DISCOVER_H

#include <Metal/Metal.hpp>
#include "small_hash_table.h"

struct MetalOps;
struct NonceInfo;

// Very simple class for a single nonce discovery metod
// Buffers are lazy init, will be reused if possible.
struct MetalContextDiscover {

    MTL::Buffer* params = nullptr;
    MTL::Buffer* hh_keys = nullptr;
    MTL::Buffer* nonces = nullptr;

    ~MetalContextDiscover() {release();}

    void init(MetalOps * metal, const SmallHashTable<bool> & keys, bool upass, uint64_t v[4]);
    void read_results(std::vector<NonceInfo> & discovered_nonces );

    void release();
};


#endif //METAL_CONTEXT_DISCOVER_H
