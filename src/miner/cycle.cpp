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

#include "cycle.h"
#include "cycle_finder.h"
#include <iostream>

Cycle::Cycle(const std::vector<EdgeExt> & nodes) {
    init(nodes);
}

void Cycle::init(const std::vector<EdgeExt> & nodes) {
    root_hashes.resize(0);
    nonces.resize(0);
    hash_2_nonce_idx.clear();
    u_hash_q.clear();
    v_hash_q.clear();
    cycle.clear();

    const uint32_t UBit = uint32_t(1) << 31;
    const uint32_t UMask = UBit - 1;

    for (const auto & n : nodes) {
        root_hashes.push_back(n.hash1);
        root_hashes.push_back(n.hash2);

        if (n.hash1 & UBit) {
            u_hash_q.insert( (n.hash1 & UMask) );
            u_hash_q.insert( (n.hash1 & UMask) ^ 0x01 );
        }
        else {
            v_hash_q.insert( n.hash1 );
            v_hash_q.insert( n.hash1 ^ 0x01 );
        }

        if (n.hash2 & UBit) {
            u_hash_q.insert( (n.hash2 & UMask) );
            u_hash_q.insert( (n.hash2 & UMask) ^ 0x01 );
        }
        else {
            v_hash_q.insert( n.hash2 );
            v_hash_q.insert( n.hash2 ^ 0x01 );
        }

        for (const auto & interim : n.interim_hashes) {
            if (interim.second.hash & UBit) {
                u_hash_q.insert( (interim.second.hash & UMask) );
                u_hash_q.insert( (interim.second.hash & UMask) ^ 0x01 );
            }
            else {
                v_hash_q.insert( interim.second.hash );
                v_hash_q.insert( interim.second.hash ^ 0x01 );
            }
        }
    }
}

void Cycle::add_nonces(const std::vector<NonceInfo> & ns, bool upass) {
    const uint32_t UBit = uint32_t(1) << 31;

    for (const auto & n : ns) {
        auto range = hash_2_nonce_idx.equal_range( n.uhash | UBit );
        bool found = false;
        for (auto it = range.first; it != range.second; ++it) {
            if (nonces[it->second].nonce == n.nonce)
                found = true;
        }
        if (found)
            continue;

        // Adding nonce
        // Add next hash into the request
        if (upass) {
            if (hash_2_nonce_idx.count(n.vhash)==0) {
                v_hash_q.insert( n.vhash ^ 0x01 );
            }
        }
        else {
            if (hash_2_nonce_idx.count(n.uhash | UBit)==0) {
                u_hash_q.insert( n.uhash ^ 0x01 );
            }
        }

        nonces.push_back(n);
        uint32_t pos = nonces.size()-1;

        hash_2_nonce_idx.insert( {n.uhash | UBit, pos} );
        hash_2_nonce_idx.insert( {n.vhash, pos} );
    }
}

struct NodePath {
    NodePath(uint32_t _nonce, uint32_t _hash, uint32_t _nonce_idx, int _prev_path_idx, uint32_t _len) :
        nonce(_nonce), hash(_hash), nonce_idx(_nonce_idx), prev_path_idx(_prev_path_idx), len(_len) {}

    uint32_t nonce;
    uint32_t hash;
    uint32_t nonce_idx;
    int prev_path_idx;
    uint32_t len;
};

static bool visited(uint32_t nonce, const std::vector<NodePath> & paths, int idx ) {
    while (idx>=0) {
        if (paths[idx].nonce==nonce)
            return true;
        idx = paths[idx].prev_path_idx;
    }
    return false;
}

// Detect cycle is possible
bool Cycle::detect_cycle(uint32_t cycle_len ) {
    if (!cycle.empty())
        return true;

    // May be possible
    assert(!root_hashes.empty());

    uint32_t start_hash = root_hashes[0];
    uint32_t end_hash = start_hash ^ 0x01;

    std::vector<NodePath> paths;

    auto range = hash_2_nonce_idx.equal_range( start_hash );
    for (auto it = range.first; it != range.second; ++it) {
        const NonceInfo & n = nonces[it->second];
        paths.push_back( NodePath(n.nonce, start_hash, it->second, -1, 1 ) );
    }

    for (int k=0; k<paths.size(); k++) {
        NodePath pth = paths[k];
        if (pth.len>cycle_len)
            continue;

        const NonceInfo & n = nonces[pth.nonce_idx];
        uint32_t next_hash = n.get_other_nonce(pth.hash);

        if (pth.len == cycle_len && next_hash == end_hash) {
            // solution is found, let's rewind it

            cycle.push_back( pth.nonce );
            while(pth.prev_path_idx>=0) {
                pth = paths[pth.prev_path_idx];
                cycle.push_back( pth.nonce );
            }
            assert(cycle.size() == cycle_len);
            break;
        }

        auto next = hash_2_nonce_idx.equal_range( next_hash ^ 0x01 );
        for (auto it = next.first; it != next.second; ++it) {
            const NonceInfo & nn = nonces[it->second];
            if (!visited(nn.nonce, paths, k))
                paths.push_back( NodePath(nn.nonce, next_hash ^ 0x01, it->second, k, pth.len+1 ) );
        }

    }
    return !cycle.empty();
}
