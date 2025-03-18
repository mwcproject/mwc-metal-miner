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

#ifndef TEST_LEAN_MINER_H
#define TEST_LEAN_MINER_H
#include <cstdint>
#include <vector>
#include "../miner/sip_hash.h"

template <uint8_t EDGE_BITS>
class TestLeanMiner {
private:
    uint64_t hash_v[4];
    std::vector<bool> U;
    std::vector<bool> V;
    std::vector<bool> active_noncesU;
    std::vector<bool> active_noncesV;
public:
    TestLeanMiner(uint64_t _hash_v[4]) {
        hash_v[0] = _hash_v[0];
        hash_v[1] = _hash_v[1];
        hash_v[2] = _hash_v[2];
        hash_v[3] = _hash_v[3];

        uint64_t EDGE_NUM = uint64_t(1) << EDGE_BITS;
        U.resize(EDGE_NUM, 1);
        V.resize(EDGE_NUM, 1);

        active_noncesU.resize(EDGE_NUM, true);
        active_noncesV.resize(EDGE_NUM, true);
    }

    void seedU(bool dump, bool updateV) {
        uint64_t nonce_num = uint64_t(1) << EDGE_BITS;
        uint32_t hash_mask = uint32_t((uint64_t(1) << EDGE_BITS)-1 );

        uint64_t EDGE_NUM = uint64_t(1) << EDGE_BITS;
        std::vector<bool> nU(EDGE_NUM, false);
        std::vector<bool> nV(EDGE_NUM, false);

        for (uint64_t i = 0; i < nonce_num; i++) {
            if (!active_noncesU[i])
                continue;

            uint32_t u = sip_hash(hash_v, i*2) & hash_mask;
            uint32_t v = sip_hash(hash_v, i*2 + 1) & hash_mask;
            if (U[u] && V[v]) {
                nU[u] = true;
                nV[v] = true;
            }
            else {
#ifdef ENABLE_TRACE
                if (dump) {
                    std::cout << std::hex << "TrimmedU: " << i << " (" << (i*2) << ")  H: " << u << " - " << v << std::endl;
                }
#endif
                active_noncesU[i] = false;
            }
        }
        U = nU;
        if (updateV) {
            V = nV;
        }
    }

    void seedV(bool dump, bool updateU) {
        uint64_t nonce_num = uint64_t(1) << EDGE_BITS;
        uint32_t hash_mask = uint32_t((uint64_t(1) << EDGE_BITS)-1 );

        uint64_t EDGE_NUM = uint64_t(1) << EDGE_BITS;
        std::vector<bool> nU(EDGE_NUM, false);
        std::vector<bool> nV(EDGE_NUM, false);

        for (uint64_t i = 0; i < nonce_num; i++) {
            if (!active_noncesV[i])
                continue;

            uint32_t u = sip_hash(hash_v, i*2) & hash_mask;
            uint32_t v = sip_hash(hash_v, i*2 + 1) & hash_mask;
            if (U[u] && V[v]) {
                nU[u] = true;
                nV[v] = true;
            }
            else {
#ifdef ENABLE_TRACE
                if (dump) {
                    std::cout << std::hex << "TrimmedV: " << i << " (" << (i*2) << ")  H: " << u << " - " << v << std::endl;
                }
#endif
                active_noncesV[i] = false;
            }
        }
        if (updateU) {
            U = nU;
        }
        V = nV;
    }

    void trimU(bool dump) {
        uint64_t EDGE_NUM = uint64_t(1) << EDGE_BITS;
        for (uint64_t i = 0; i < EDGE_NUM; i+=2) {
            if (U[i] != U[i+1]) {
                U[i] = false;
                U[i+1] = false;
            }
        }
        // Updating active_noncesU
        uint64_t nonce_num = uint64_t(1) << EDGE_BITS;
        uint32_t hash_mask = uint32_t((uint64_t(1) << EDGE_BITS)-1 );

        for (uint64_t i = 0; i < nonce_num; i++) {
            if (!active_noncesU[i])
                continue;

            uint32_t u = sip_hash(hash_v, i*2) & hash_mask;
            if (!U[u]) {
                if (active_noncesU[i]) {
#ifdef ENABLE_TRACE
                    if (dump) {
                        uint32_t v = sip_hash(hash_v, i*2+1) & hash_mask;
                        std::cout << std::hex << "TrimmedU: " << i << " (" << (i*2) << ")  H: " << u << " - " << v << std::endl;
                    }
#endif
                    active_noncesU[i] = false;
                }
            }
        }
    }

    void trimV(bool dump) {
        uint64_t EDGE_NUM = uint64_t(1) << EDGE_BITS;
        for (uint64_t i = 0; i < EDGE_NUM; i+=2) {
            if (V[i] != V[i+1]) {
                V[i] = false;
                V[i+1] = false;
            }
        }

        // Updating active_noncesU
        uint64_t nonce_num = uint64_t(1) << EDGE_BITS;
        uint32_t hash_mask = uint32_t((uint64_t(1) << EDGE_BITS)-1 );

        for (uint64_t i = 0; i < nonce_num; i++) {
            if (!active_noncesV[i])
                continue;

            uint32_t v = sip_hash(hash_v, i*2+1) & hash_mask;
            if (!V[v]) {
                if (active_noncesV[i]) {
#ifdef ENABLE_TRACE
                    if (dump) {
                        uint32_t u = sip_hash(hash_v, i*2) & hash_mask;
                        std::cout << std::hex << "TrimmedU: " << i << " (" << (i*2) << ")  H: " << u << " - " << v << std::endl;
                    }
#endif
                    active_noncesV[i] = false;
                }
            }
        }

    }

    const std::vector<bool> & getU() const {return U;}
    const std::vector<bool> & getV() const {return V;}

    const std::vector<bool> & get_active_noncesU() const {return active_noncesU;}
    const std::vector<bool> & get_active_noncesV() const {return active_noncesV;}


};

#endif //TEST_LEAN_MINER_H
