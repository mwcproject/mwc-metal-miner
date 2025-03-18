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
#define ENABLE_TESTS

#include "../miner/cuckatoo.h"
#include "cuckatoo_test.h"

void test_cuckatoo_12_4b() {
    uint64_t v[4] = { 0x736f6d6570736575ULL, 0x646f72616e646f6dULL,0x6c7967656e657261ULL,0x7465646279746573ULL};

    CuckatooSolver<12,2,4,6> solver;
    solver.set_threads(1);

    // let's try to find some solutions with edges, at least 2
    for (int i = 0; i < 1000; i++) {
        v[i%4]++;

        solver.set_hash(v);
        std::vector<CycleSolution> graphs;
        solver.build_graph(graphs, false, true);
        std::cout << "Found solutions " << graphs.size() << " for i=" << i << std::endl;
    }
}

void test_cuckatoo_12_5b() {
    uint64_t v[4] = { 0x736f6d6570736575ULL, 0x646f72616e646f6dULL,0x6c7967656e657261ULL,0x7465646279746573ULL};

    CuckatooSolver<12,2,5,6> solver;
    solver.set_threads(1);

    // let's try to find some solutions with edges, at least 2
    for (int i = 0; i < 1000; i++) {
        v[i%4]++;
        solver.set_hash(v);
        std::vector<CycleSolution> graphs;
        solver.build_graph(graphs,false, true);
        std::cout << "Found non empty solution " << graphs.size() << " for i=" << i << std::endl;
    }
}

void test_cuckatoo_25() {
    uint64_t v[4] = { 0x736f6d6570736575ULL, 0x646f72616e646f6dULL,0x6c7967656e657261ULL,0x7465646279746573ULL};

    for (int i = 0; i < 50; i++) {
        // 256 buckets
        v[i%4]++;
        CuckatooSolver<25,8,5,42> solver;
        // using all threads here
        solver.set_hash(v);
        std::vector<CycleSolution> graphs;
        solver.build_graph(graphs, false, false);
        if (graphs.size()>0) {
            std::cout << "Found non empty solution " << graphs.size() << " for i=" << i << std::endl;
        }
    }
}

void test_cuckatoo_29() {
    uint64_t v[4] = { 0x736f6d6570736575ULL, 0x646f72616e646f6dULL,0x6c7967656e657261ULL,0x7465646279746573ULL};

    for (int i = 0; i < 10; i++) {
        // 256 buckets
        v[i%4]++;
        CuckatooSolver<29,8,5,42> solver;
        // using all threads here
        solver.set_hash(v);
        std::vector<CycleSolution> graphs;
        solver.build_graph( graphs, false, false);
        if (graphs.size()>0) {
            std::cout << "Found non empty C29 solution " << graphs.size() << " for i=" << i << std::endl;
        }
    }
}


void test_cuckatoo_15_with_solution() {
    // it is a known solution from the mwc-node test for C15
    uint64_t v[4] = {  0x8caefe3dad6bc1f2ULL, 0xb5f021dbfa151bc5ULL, 0x626dd2f4c6f35c1fULL, 0x19c9434751e9d6e3ULL };
    std::vector<uint32_t> known_cycle = {
        // first known solutions (nonces)
        0x116, 0x3ae,0x6e6, 0x982, 0xb40, 0xc4f, 0xe06, 0xebe, 0x14d5, 0x1840, 0x1d7a, 0x1e02,
            0x2055, 0x2472, 0x269b, 0x26d2, 0x2b53, 0x2cfb, 0x2da0, 0x32f8, 0x3379, 0x374f, 0x3892, 0x402d,
            0x40ce, 0x4bc9, 0x4c88, 0x4fd7, 0x5266, 0x5be7, 0x6120, 0x668c, 0x6982, 0x6ba3, 0x6bec, 0x6d73,
            0x6f4d, 0x705b, 0x7330, 0x735d, 0x748a, 0x7943,
        // second known solution (nonces)
        0x30d, 0x3ae, 0x6e6,0xb2a,0x14ba,0x15fd,0x1d06,0x1d7a,0x2055,0x2081,0x2472,0x259f,0x269b,0x2bb8,0x2cfb,
        0x374f,0x3a2e,0x402d,0x43e5,0x47b4,0x48bc,0x4bc9,0x4d0e,0x4dc4,0x4fd7,0x5233,0x5268,0x5334,0x5d59,0x5ebd,
        0x6120,0x6131,0x6597,0x6d5d,0x6dbd,0x6e03,0x6ee6,0x70b2,0x7177,0x735d,0x743e,0x7f69
    };

    CuckatooSolver<15,2,4,42> solver;
    solver.set_hash(v);
    solver.set_threads(1);
    std::vector<CycleSolution> graphs;
    solver.build_graph(graphs, false, true, known_cycle);
    assert(graphs.size()==2); // there are two unique solutions
    std::vector<uint64_t> res_nonces;
    std::vector<uint8_t> ret_hash = solver.resolve_found_to_nonces(graphs, res_nonces);
    assert(!ret_hash.empty());
    assert(res_nonces.size()==42);
    assert(res_nonces[0]==781); // it is a best solution (smallest hash)
}

