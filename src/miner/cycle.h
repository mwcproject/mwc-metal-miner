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

#ifndef CYCLE_H
#define CYCLE_H

#include <vector>
#include <assert.h>
#include <unordered_map>
#include <unordered_set>
#include <set>

struct EdgeExt;

/**
 * Data about single edge nonce
 */
struct NonceInfo {
    uint32_t nonce; // Nonce that represent the edge
    uint32_t uhash; // U hash (U node)
    uint32_t vhash; // V hash (V node)

    NonceInfo(uint32_t _nonce, uint32_t _uhash, uint32_t _vhash) :
                nonce(_nonce), uhash(_uhash), vhash(_vhash) {}
    NonceInfo() = default;
    NonceInfo(const NonceInfo &) = default;

    inline uint32_t get_other_nonce(uint32_t hash) const {
        const uint32_t UBit = (uint32_t(1) << 31);
        assert(hash==(uhash|UBit) || hash==vhash);
        return hash==vhash ? (uhash | UBit) : vhash;
    }
};

/**
 * Cuckatoo cycle - solution
 */
struct Cycle {
    // Hashes that was fault at pahse 2. It is what belongs to the cycle
    std::vector<uint32_t> root_hashes; // hashes that are know from the nodes
    // Recovered nonces. There are nonces that are belong to cycle and extras.
    std::vector<NonceInfo> nonces;
    // hash 2 nonce index
    std::unordered_multimap<uint32_t, uint32_t> hash_2_nonce_idx;

    // Half hashes without Ubit to process.
    // It is tasks for next discovery stage.
    std::unordered_set<uint32_t> u_hash_q;
    std::unordered_set<uint32_t> v_hash_q;

    // found solution with nonces. Might be empty if this Cycle was false positive.
    std::vector<uint32_t> cycle;

    explicit Cycle(const std::vector<EdgeExt> & nodes);
    Cycle() = default;
    Cycle(const Cycle&) = default;
    Cycle(Cycle&&) = default;
    ~Cycle() {}

    void init(const std::vector<EdgeExt> & nodes);

    std::unordered_set<uint32_t> & get_hashes_to_discover(bool upass) {
        return upass ? u_hash_q : v_hash_q;
    }

    // append a new nonce.
    void add_nonces(const std::vector<NonceInfo> & nonces, bool upass);

    // Detect cycle is possible
    bool detect_cycle(uint32_t cycle_len);

    const std::vector<uint32_t> & get_cycle() const {return cycle;}
};


#endif //CYCLE_H
