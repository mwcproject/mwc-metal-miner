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

#ifndef METAL_CODE_H
#define METAL_CODE_H

#include <string>

// The metal code is here
std::string generate_metal_sources( uint32_t EDGE_BITS, uint32_t CYCLE_LEN, uint32_t BUCKETS_NUM, uint32_t BUCKETS_8B_NUM,
    const uint32_t EDGE_PER_BUCKET_PREV_PRIME, const uint32_t RESULTING_EDGE_PER_BUCKET,
    uint32_t COLLAPSE_STASH_SIZE, uint32_t COLLAPSE_STASH_COUNTER );

#endif //METAL_CODE_H
