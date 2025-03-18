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

#ifndef MINER_DATA_H
#define MINER_DATA_H

#include <cstdint>

struct HashNonce {
    uint32_t hash;
    uint32_t nonce; // nonce is unsigned BUT expected to have internal offset, so incremental could be negative

    HashNonce() {}
    HashNonce(uint32_t _hash, uint32_t _nonce) : hash(_hash), nonce(_nonce) {}

    inline bool is_zero() const {return hash==0 && nonce==0;}
};

#endif //MINER_DATA_H
