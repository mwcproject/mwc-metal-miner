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

#ifndef CUCKATOO_H
#define CUCKATOO_H

#include <iomanip>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include "row_builder.h"
#include "trimmer.h"
#include "blake.h"

#ifdef ENABLE_TESTS
#include "../tests/test_lean_miner.h"
#endif

#ifdef ENABLE_TRACE
struct SrcData {
    uint32_t nonce;
    uint32_t hashU;
    uint32_t hashV;
};
#endif

struct CycleSolution {
    std::vector<std::pair<uint32_t,uint32_t>> uv_pairs;
};

// Dead mean that was tested and it didn't belong to the graph. No reasons to go this way
#define NODE_DEAD  0x40000000
// All Visited data mask
#define NODE_VISITED_DATA  0x3FFFFFFF
// Letf/right branch
#define NODE_LEFT 0x20000000
// graph extraction ID
#define NODE_VISITED_ID  0x00FFFFFF
// Level from the root. For 42 graph, max value is 42/2 = 21.
#define NODE_VISITED_LEVEL 0x1F000000
#define NODE_LEVEL_SHIFT 24

const uint64_t V_HASH_PLUS = uint64_t(1)<<40;
struct Node {
    uint64_t from;
    uint64_t to;

    inline void init(uint64_t _from, uint64_t _to) {
        from = _from;
        to = _to;
    }

    inline bool is_dead() const {return (state & NODE_DEAD)!=0;}

    inline uint32_t get_visited_id() const {return state & NODE_VISITED_ID;}
    inline bool is_visited(uint32_t visited_id) const {return /*(state & NODE_DEAD)!=0 ||*/ (state & NODE_VISITED_ID)==visited_id;}
    inline bool is_visited_with(uint32_t visited_id, bool left, uint8_t level) const { return get_visited_id()==visited_id && is_left()==left && get_level()==level; }
    inline bool is_left() const {return (state & NODE_LEFT)!=0;}
    inline uint32_t get_level() const {return (state & NODE_VISITED_LEVEL) >> NODE_LEVEL_SHIFT;}
    inline void set_dead() {state |= NODE_DEAD;}

    void set_visited(uint32_t visitedId, bool left, uint32_t level) {
        assert(visitedId < NODE_VISITED_ID);
        state = (state & ~NODE_VISITED_DATA) | visitedId | (left ? NODE_LEFT : 0) | (level << NODE_LEVEL_SHIFT);
    }


private:
    uint32_t state = 0;
};

uint32_t next_prime(uint32_t start);

template <uint8_t EDGE_BITS>
std::vector<uint8_t> compress_nonces(const std::vector<uint64_t>& nonces) {
    std::vector<uint8_t> compress((nonces.size() * EDGE_BITS + 7)/8, 0); // Initialize the compressed vector with zeros
    int currentBit = 0; // Track the current bit position in the compress vector
    for (uint64_t nonce : nonces) {
        // Extract the low EDGE_BIT bits from the nonce
        assert(nonce < ((1ULL << EDGE_BITS) - 1));
        assert(EDGE_BITS>8);

        int need_write = EDGE_BITS;
        while (need_write>0) {
            if (currentBit % 8 == 0) {
                // full byte can be written
                compress[currentBit / 8] = uint8_t(nonce);
                nonce >>= 8;
                currentBit += std::min(8, need_write);
                need_write -= 8;
            }
            else {
                int can_write = 8 - currentBit % 8;
                assert(can_write<8);
                compress[currentBit / 8] |= uint8_t(nonce << (8-can_write));
                currentBit += std::min(can_write, need_write);
                nonce >>= can_write;
                need_write -= can_write;
            }
        }
    }
    assert(currentBit>=compress.size()*8-7 && currentBit<=compress.size()*8);
    return compress;
}

class Solver {
 public:
    virtual ~Solver() {}
    virtual void set_hash(uint64_t _hash_v[4]) = 0;
    virtual void set_threads(int num) = 0;
    virtual void build_graph(std::vector<CycleSolution> & res_graphs,  bool dump, bool test_with_lean = true,
            const std::vector<uint32_t> & known_cycle = {}) = 0;
    virtual std::vector<uint8_t>  resolve_found_to_nonces(std::vector<CycleSolution> & solutions, std::vector<uint64_t> & res_nonces) = 0;
};

// cuckatoo solution
template <uint8_t EDGE_BITS, uint8_t BUCKET_BITS, uint ELEMENT_SIZE, uint GRAPH_SIZE>
class CuckatooSolver : public Solver {
private:
    static const uint8_t MASK_BITS = 4;
    static const uint32_t BIT_PACKER_INDEX_INTERVAL = 64;

    MetalOps metal;
    uint64_t hash_v[4];
    int threads_num;
public:
    CuckatooSolver() : metal(EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS) {
        threads_num = std::thread::hardware_concurrency();
    }

    virtual ~CuckatooSolver() {}

    virtual void set_hash(uint64_t _hash_v[4]) override {
        hash_v[0] = _hash_v[0];
        hash_v[1] = _hash_v[1];
        hash_v[2] = _hash_v[2];
        hash_v[3] = _hash_v[3];
    }

    virtual void set_threads(int num) override {
        threads_num = num;
    }

    // returns back number of
    virtual void build_graph(std::vector<CycleSolution> & res_graphs,  bool dump, bool test_with_lean = true,
            const std::vector<uint32_t> & known_cycle = {}) override {
#ifdef ENABLE_TRACE
        std::vector<uint32_t> expected_u;
        std::vector<uint32_t> expected_v;
        uint32_t hash_mask = uint32_t((uint64_t(1) << EDGE_BITS) - 1);
        if (!known_cycle.empty()) {
                if (dump) {
                    std::cout << "Known cycle data: " << std::endl;
                }
                for (auto nonce : known_cycle) {
                    uint32_t u = sip_hash(hash_v, nonce*2) & hash_mask;
                    uint32_t v = sip_hash(hash_v, nonce*2+1) & hash_mask;
                    expected_u.push_back(u);
                    expected_v.push_back(v);

                    if (dump) {
                        std::cout << std::hex << nonce << "  " << u << "  " << v << std::endl;
                    }
                }
        }
        std::sort(expected_u.begin(), expected_u.end());
        std::sort(expected_v.begin(), expected_v.end());

        if (dump)
        {
            // let's generate and pront the row data
            std::vector<SrcData> data;
            uint64_t nonce_num = uint64_t(1) << EDGE_BITS;
            for (uint64_t i = 0; i < nonce_num; i++) {
                SrcData dt;
                dt.nonce = i;
                dt.hashU = sip_hash(hash_v, i*2) & hash_mask;
                dt.hashV = sip_hash(hash_v, i*2+1) & hash_mask;
                data.push_back(dt);
            }
            std::sort(data.begin(), data.end(), [](const SrcData& a, const SrcData& b) {return a.nonce < b.nonce;});
            std::cout << "Source data by NONCES" << std::endl;
            for (const auto & d : data) {
                std::cout << std::hex << d.nonce << "  " << d.hashU << "  " << d.hashV << std::endl;
            }
            std::sort(data.begin(), data.end(), [](const SrcData& a, const SrcData& b) {return a.hashU < b.hashU;});
            std::cout << std::endl << "Source data by HASH U" << std::endl;
            for (const auto & d : data) {
                std::cout << std::hex << d.nonce << "  " << d.hashU << "  " << d.hashV << std::endl;
            }
            std::sort(data.begin(), data.end(), [](const SrcData& a, const SrcData& b) {return a.hashV < b.hashV;});
            std::cout << std::endl << "Source data by HASH V" << std::endl;
            for (const auto & d : data) {
                std::cout << std::hex << d.nonce << "  " << d.hashU << "  " << d.hashV << std::endl;
            }
            std::cout << std::endl;
        }
#endif
        std::vector<MTL::Buffer *> trimmed_U;
//        switch_enable_data_dumping(true);
        RowBuilder<EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> U(metal, hash_v, 0, threads_num, trimmed_U);
//        switch_enable_data_dumping(false);
        std::vector<MTL::Buffer *> trimmed_V;
        RowBuilder<EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> V(metal, hash_v, 1, threads_num, trimmed_V);
#ifdef ENABLE_TESTS
#ifdef ENABLE_TRACE
        if (dump) {
            dump_rowmask_and_trimmed("U RAW", U, trimmed_U);
            dump_rowmask_and_trimmed("V Raw", V, trimmed_V);
        }
#endif

        // If tests are enable, let's use in parallel the lean implementation, it is slow on CPU but easy to implement, so let's use for the testing
        TestLeanMiner<EDGE_BITS> lean_miner(hash_v);
        if (test_with_lean) {
            lean_miner.seedU(dump, true);
            lean_miner.trimU(dump);
            lean_miner.trimV(false);

            test_expected(expected_u, expected_v,U,V,lean_miner);

            // Masks are supposed to match. Checking them
            compare_masks(U, lean_miner.getU());
            compare_masks(V, lean_miner.getV());
            check_trimmed(trimmed_U, lean_miner.get_active_noncesU());
            check_trimmed(trimmed_V, lean_miner.get_active_noncesV());
        }
#endif

        // First trim round, different from the next ones because of trimmed_U & trimmed_V
        uint32_t trimmedU = trim(metal,
                  &U,
                  &trimmed_V,
                  (1 << (EDGE_BITS - MASK_BITS)) / 4,
                  hash_v, 0,
                  threads_num);
#ifdef ENABLE_TESTS
        if (dump) {
            dump_rowmask_and_trimmed("U Step1", U, trimmed_U);
        }

        if (test_with_lean) {
            lean_miner.seedU(dump, false);
            lean_miner.trimU(dump);

            test_expected(expected_u, expected_v,U,V,lean_miner);

            compare_masks(U, lean_miner.getU());
        }
#endif
        trimmed_U.insert(trimmed_U.end(), trimmed_V.begin(), trimmed_V.end());
        trimmed_V.clear();
        [[maybe_unused]] uint32_t trimmedV = trim(metal,
                  &V,
                  &trimmed_U,
                  (1 << (EDGE_BITS - MASK_BITS)) / 4,
                  hash_v, 1,
                  threads_num);
#ifdef ENABLE_TESTS
        if (dump) {
            dump_rowmask_and_trimmed("V Step1", V, trimmed_U);
        }

        if (test_with_lean) {
            lean_miner.seedV(dump, false);
            lean_miner.trimV(dump);
            test_expected(expected_u, expected_v,U,V,lean_miner);
            compare_masks(V, lean_miner.getV());
        }

        if (dump) {
            std::cout << "Step 1: Trimmed U: " << trimmedU << "  V: " << trimmedV << std::endl;
        }
#endif

        // main loop with trimming, trimmed_U is reused for both
        while (trimmedU + trimmedU > 0) {
            uint32_t trim_nonce_per_bucket = ((trimmedU + trimmedU) / (1 << BUCKET_BITS)  + 64)  & ~63ULL;
            assert( trim_nonce_per_bucket%64==0);
            assert( trim_nonce_per_bucket > 0);

            trimmedU = trim(metal,
                              &U,
                              &trimmed_U,
                              trim_nonce_per_bucket,
                              hash_v, 0,
                              threads_num);

#ifdef ENABLE_TESTS
            if (dump) {
                dump_rowmask_and_trimmed("U Step2+", U, trimmed_U);
            }


            if (test_with_lean) {
                lean_miner.seedU(dump, false);
                lean_miner.trimU(dump);
                test_expected(expected_u, expected_v,U,V,lean_miner);
                compare_masks(U, lean_miner.getU());
            }
#endif

            trimmedV = trim(metal,
                      &V,
                      &trimmed_U,
                      trim_nonce_per_bucket,
                      hash_v, 1,
                      threads_num);

#ifdef ENABLE_TESTS
            if (test_with_lean) {
                lean_miner.seedV(dump, false);
                lean_miner.trimV(dump);
                test_expected(expected_u, expected_v,U,V,lean_miner);
                compare_masks(V, lean_miner.getV());
            }

            if (dump) {
                std::cout << "Step 2+: Trimmed U: " << trimmedU << "  V: " << trimmedV << std::endl;
            }
#endif
        }

#ifdef ENABLE_TESTS
        if (test_with_lean) {
            compare_masks(U, lean_miner.getU());
            compare_masks(V, lean_miner.getV());
            test_expected(expected_u, expected_v,U,V,lean_miner);
       }
#endif

        search_for_graphs( U,V,res_graphs );

#ifndef NDEBUG
        { // Verify if U/V hashes are really belong to the U/V data
            for (const auto &sol : res_graphs) {
                for (const auto &uv : sol.uv_pairs) {
                    uint32_t u = std::get<0>(uv);
                    uint32_t v = std::get<1>(uv);
                    assert(U.has_mask(u));
                    assert(V.has_mask(v));
                }
            }
        }
#endif
    }

    // return resulting hash
    virtual std::vector<uint8_t> resolve_found_to_nonces(std::vector<CycleSolution> & solutions, std::vector<uint64_t> & res_nonces) override {
        if (solutions.empty())
            return std::vector<uint8_t>();

        // we want to find the hash so all cycle solutions can be represented by multiple hash_maps, scare factor is 3:
        const uint32_t U_size = next_prime(solutions.size() * GRAPH_SIZE/2 * 3);
        MTL::Buffer* U_hash_table = metal.create_buffer( U_size * sizeof(uint32_t) );
        // NN - nonce-nonce
        const uint32_t VVNN_size = U_size*2;
        MTL::Buffer* VV_table = metal.create_buffer( VVNN_size * sizeof(uint32_t));
        MTL::Buffer* NN_table = metal.create_buffer( VVNN_size * sizeof(uint32_t)  );

        uint32_t * U_values = (uint32_t *) U_hash_table->contents();
        uint32_t * VV_values = (uint32_t *) VV_table->contents();
        uint32_t * NN_values = (uint32_t *) NN_table->contents();

        memset( U_values, 0, U_size * sizeof(uint32_t) );
        memset( VV_values, 0, VVNN_size * sizeof(uint32_t) );
        memset( NN_values, 0, VVNN_size * sizeof(uint32_t) );

        for (auto & sol : solutions) {
            for (const auto & uv : sol.uv_pairs ) {
                auto u = std::get<0>( uv );
                auto v = std::get<1>( uv );
                auto ui = u & 0x1;
                u &= 0xFFFFFFFE;
                uint32_t table_idx = u % U_size;
                for (int i=0;i<U_size; i++, table_idx++) {
                    uint32_t idx = table_idx % U_size;
                    if (U_values[idx]==0) {
                        // found spot, fill the data
                        U_values[idx] = u;
                        VV_values[idx*2 + ui] = v;
                        break;
                    }
                    else if (U_values[idx]==u && VV_values[idx*2 + ui]==0) {
                        VV_values[idx*2 + ui] = v;
                        break;
                    }
                }
            }
        }

        // We are ready, can do metal job
        metal.uv_to_nonces(hash_v, U_hash_table, VV_table, NN_table);

        // blake hash here is 32 bytes
        std::vector<uint8_t> best_hash(32,0xff);
        std::vector<uint8_t> hash(32);
        for (int i=0;i<solutions.size();i++) {
            const auto & uv_pairs = solutions[i].uv_pairs;
            std::vector<uint64_t> nonces;
            nonces.reserve(GRAPH_SIZE);

            assert(uv_pairs.size()==GRAPH_SIZE);
            for (const auto & uv : uv_pairs ) {
                auto u = std::get<0>( uv );
                auto v = std::get<1>( uv );
                auto ui = u & 0x1;
                u &= 0xFFFFFFFE;
                uint32_t table_idx = u % U_size;
                for (int i=0;i<U_size; i++, table_idx++) {
                    uint32_t idx = table_idx % U_size;
                    if (U_values[idx]==u && VV_values[idx*2+ui]==v) {
                        uint32_t nc = NN_values[idx*2 + ui];
                        assert(nc >0); // not found?
                        nonces.push_back(nc);
                        break;
                    }
                    if (U_values[idx]==0) {
                        assert(false); // nonce is not in the table?
                        break;
                    }
                }
            }
            assert(nonces.size() == GRAPH_SIZE);
            std::sort(nonces.begin(), nonces.end());

            std::vector<uint8_t> nonces_compressed = compress_nonces<EDGE_BITS>(nonces);

            blake2b(hash.data(), nonces_compressed.data(), nonces_compressed.size());
            for (int e=0;e<32;e++) {
                if (hash[e] < best_hash[e]) {
                    best_hash = hash;
                    res_nonces = nonces;
                    break;
                }
                else if (hash[e] > best_hash[e]) {
                    break; // current solution is better
                }
            }
        }

        U_hash_table->release();
        VV_table->release();
        NN_table->release();

        return best_hash;
    }


private:
    // return true if solution is found. If there are several solutions, returns the best one (highest hash).
    void search_for_graphs( RowBuilder<EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> &U,
                            RowBuilder<EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> &V,
                            std::vector<CycleSolution> & res_graph) {
        // Processing everything in a single thread, will optimize after
        std::vector<uint64_t> & mask = U.get_mask();
        std::vector<Bucket<EDGE_BITS,BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL>> & bucketsU = U.get_buckets();

        const uint32_t BUCKET_SIZE = 1 << (EDGE_BITS - BUCKET_BITS*2);
        const uint32_t HASH_MASK = (uint64_t(1) << EDGE_BITS)-1;

        std::vector<Node> nodes;
        nodes.reserve(1000);

        uint32_t edge_idx = 0;
        uint32_t nonces_buf[1<<MASK_BITS];
        // scanning u mask first because we want to resolve u/v collisions (apparently there are many of them exist)
        for (uint64_t m : mask) {
            for (uint i = 0; i < 8; i++) {
                if (m & 0xff) {
                    assert((m & 0x0f) !=0);
                    assert((m & 0xf0) !=0);

                    const uint32_t bucketIdx = edge_idx / BUCKET_SIZE;
                    IntervalsIndex index;
                    uint32_t num = bucketsU[bucketIdx].get_nonces(edge_idx % BUCKET_SIZE, index, nonces_buf);
                    assert(num>0);
                    for (uint32_t k = 0; k < num; k++) {
                        uint32_t nc = nonces_buf[k];
                        uint32_t hash2 = sip_hash(hash_v, uint64_t(nc)*2+1) & HASH_MASK;
                        if (V.has_mask(hash2)) {
                            nodes.emplace_back().init(uint64_t(hash2) + V_HASH_PLUS, edge_idx);
                            nodes.emplace_back().init(edge_idx, uint64_t(hash2) + V_HASH_PLUS);
                        }
                    }
                    if (nodes.size()) {
                        assert(nodes.size()>0);
                        assert(nodes.back().from == edge_idx); // At least one should be found
                    }

                    num = bucketsU[bucketIdx].get_nonces(edge_idx % BUCKET_SIZE + 1, index, nonces_buf);
                    assert(num>0);
                    for (uint32_t k = 0; k < num; k++) {
                        uint32_t nc = nonces_buf[k];
                        uint32_t hash2 = sip_hash(hash_v, uint64_t(nc)*2+1) & HASH_MASK;
                        if (V.has_mask(hash2)) {
                            nodes.emplace_back().init(uint64_t(hash2) + V_HASH_PLUS, edge_idx+1);
                            nodes.emplace_back().init(edge_idx+1, uint64_t(hash2) + V_HASH_PLUS);
                        }
                    }
                    if (nodes.size()) {
                        assert(nodes.size()>0);
                        assert(nodes.back().from == edge_idx+1); // At least one should be found
                    }
                }
                m >>= 8;
                edge_idx+=2;
            }
        }

        if (nodes.size()<GRAPH_SIZE*2) {
            return;
        }

        std::sort(nodes.begin(), nodes.end(), [](const Node & n1, const Node & n2) {
            return n1.from < n2.from;
        });

        std::unordered_map<uint64_t, uint32_t> hash2index;
        uint64_t prev_hash = UINT_MAX;
        for (uint32_t i = 0; i < nodes.size(); i++) {
            if (prev_hash != nodes[i].from) {
                prev_hash = nodes[i].from;
                hash2index[prev_hash] = i;
            }
        }
#ifndef NDEBUG
        for (const auto & node : nodes ) {
            assert(node.from >= V_HASH_PLUS || node.to>=V_HASH_PLUS); // Expected one U, one V. Might fail on 0 has value
             assert(hash2index.count(node.from)==1);
             assert(hash2index.count(node.to)==1);
        }
#endif

        // now we can go forward and try to build a graph. We are looking for the 42 edges graph, so we can skip every 4 hashes and 95% chances to hit the graph
        double target_prob = 0.99; // Let's be OK with a chance to hit graph with 99%
        assert(nodes.size() >= GRAPH_SIZE);
        // X - how many choices I need to do to not miss the group of the graph
        double w = (nodes.size()/2 - (GRAPH_SIZE/4+1)) / static_cast<double>(nodes.size()/2);
        double X = std::log(1 - target_prob) / std::log(w);
        uint32_t step = std::max(1U, uint32_t(nodes.size()/X));
        assert(step>=1);

        std::vector<uint64_t> nextR1;
        std::vector<uint64_t> nextR2;
        std::vector<uint64_t> nextL1;
        std::vector<uint64_t> nextL2;
        nextL1.reserve(16);
        nextL2.reserve(16);
        nextR1.reserve(16);
        nextR2.reserve(16);
        uint32_t visitedId = 0;
        std::vector<uint64_t> * R = &nextR1;
        std::vector<uint64_t> * Rn = &nextR2;
        std::vector<uint64_t> * L = &nextL1;
        std::vector<uint64_t> * Ln = &nextL2;

        std::set<uint32_t> solutions_sums;

        uint nodes_sz = nodes.size();
        for (uint32_t i = 1; i < nodes_sz/2;) {
            // U nodes
            assert(nodes[i].from < V_HASH_PLUS );
            if (nodes[i].from%2==1) {
                i++;
                continue;
            }

            uint64_t from = nodes[i].from;

            if (from==nodes[i-1].from) {
                // skipping to the next
                for (i++; i<nodes_sz && from==nodes[i].from; i++ ) {}
                continue;
            }

            if ( nodes[hash2index.at(from & ~1ULL)].is_dead()) {
                i+=step; // it is dead starting point, no reasons to start
                continue;
            }

            visitedId++;

            R->resize(0);
            R->push_back(from);
            assert(hash2index.count(from^1)==1);
            L->resize(0);
            L->push_back(from^1);
            nodes[hash2index.at(from & ~1ULL)].set_visited(visitedId, true, 0);
            nodes[hash2index.at(from & ~1ULL)].set_dead();

            bool LDead = true;
            bool RDead = true;

            assert(GRAPH_SIZE%2 == 0);
            assert(GRAPH_SIZE/2%2 == 1); // important for how we are going in reverce order (last point MUST be V, first is U)
            uint8_t level = 1;
            for (; level < GRAPH_SIZE/2; level++) {
                // Updating R branch
                next_step( R, Rn, visitedId, level, false,
                                hash2index,nodes);
                next_step( L, Ln, visitedId, level, true,
                                hash2index,nodes);
                if (Rn->size()>1)
                    RDead = false;
                if (Ln->size()>1)
                    LDead = false;

                std::swap(R, Rn);
                std::swap(L, Ln);
                if (R->empty() || L->empty())
                    break;

                if (RDead) {
                    assert(R->size()==1);
                    nodes[hash2index.at( R->front() & ~1ULL)].set_dead();
                }
                if (LDead) {
                    assert(L->size()==1);
                    nodes[hash2index.at( L->front() & ~1ULL)].set_dead();
                }
            }
            // The last step. Here both R & L must have the common child for solution.
            if (!R->empty() && !L->empty()) {
                // Make R move first
                next_step( R, Rn, visitedId, level, false,
                                hash2index,nodes);
                // L step let's do maunally, sinse it is expected that the target is market as visited with level 'level'
                for (uint64_t cur_id : *L) {
                    uint idx = hash2index.at(cur_id);
                    bool found = false;
                    for (auto it=nodes.begin()+idx; it!=nodes.end(); it++) {
                        if (it->from != cur_id)
                            break;
                        uint32_t to_idx = hash2index.at(it->to & ~1ULL);
                        if (nodes[to_idx].get_visited_id() == visitedId && nodes[to_idx].get_level() == level) {
                            // Also both Rn and Ln Must point into different part of the pair.
                            if (std::count(Rn->begin(), Rn->end(), it->to) > 0) {
                                // Solution is found!!!!
                                assert(!nodes[to_idx].is_left()); // it has to be right branch.
                                std::vector<std::pair<uint32_t,uint32_t>> uv_pairs = extract_solution( hash2index, nodes, visitedId, it->to );

                                uint32_t u_sum = 0;
                                for (const auto & pair : uv_pairs) {
                                    u_sum += std::get<0>(pair);
                                }

                                if ( std::get<1>(solutions_sums.insert(u_sum))) {
                                    res_graph.emplace_back().uv_pairs = uv_pairs;
                                }

                                found = true;
                                break;
                            }
                        }
                    }
                    if (found)
                        break;
                }
            }
            i++;
        }
    }

    // res_hashes must be only U hashes. We will map them back to nonces after
    std::vector<std::pair<uint32_t,uint32_t>> extract_solution(const std::unordered_map<uint64_t, uint32_t> & hash2index,
                std::vector<Node> & nodes,
                uint32_t visitedId,
                uint64_t last_nodeR) {

        std::vector<std::pair<uint32_t,uint32_t>> res_uv_hashes;
        res_uv_hashes.reserve(GRAPH_SIZE);

        assert(nodes[hash2index.at(last_nodeR & ~1ULL)].get_visited_id() == visitedId);
        assert(nodes[hash2index.at(last_nodeR & ~1ULL)].get_level() == GRAPH_SIZE/2);

        uint64_t nextL = last_nodeR^1;
        uint64_t nextR = last_nodeR;

        for (uint8_t level = GRAPH_SIZE/2-1; level>0; level--) {
            uint64_t NnextL = find_prev_node(hash2index,nodes, nextL^1, visitedId, true, level);
            uint64_t NnextR = find_prev_node(hash2index,nodes, nextR^1, visitedId, false, level );
            if (nextL < V_HASH_PLUS ) {
                assert(nextR < V_HASH_PLUS);
                assert(NnextL >= V_HASH_PLUS);
                assert(NnextR >= V_HASH_PLUS);
                res_uv_hashes.push_back(std::pair<uint32_t, uint32_t>(nextL^1, uint32_t(NnextL)));
                res_uv_hashes.push_back(std::pair<uint32_t, uint32_t>(nextR^1, uint32_t(NnextR)));
            }
            else {
                assert(nextR >= V_HASH_PLUS);
                assert(NnextL < V_HASH_PLUS);
                assert(NnextR < V_HASH_PLUS);
                res_uv_hashes.push_back(std::pair<uint32_t, uint32_t>(NnextL, uint32_t(nextL^1)));
                res_uv_hashes.push_back(std::pair<uint32_t, uint32_t>(NnextR, uint32_t(nextR^1)));
            }
            nextL = NnextL;
            nextR = NnextR;
            assert(nodes[hash2index.at(nextL & ~1ULL)].get_visited_id() == visitedId);
            assert(nodes[hash2index.at(nextL & ~1ULL)].get_level() == level);
            assert(nodes[hash2index.at(nextR & ~1ULL)].get_visited_id() == visitedId);
            assert(nodes[hash2index.at(nextR & ~1ULL)].get_level() == level);
        }

        // The root node, mast have left flag and be U
        uint64_t NnextL = find_prev_node(hash2index,nodes, nextL^1, visitedId, true, 0);
#ifndef NDEBUG
        uint64_t NnextR = find_prev_node(hash2index,nodes, nextR^1, visitedId, true, 0 );
#endif
        assert( (NnextL^1) == NnextR);
        assert(nodes[hash2index.at(NnextL & ~1ULL)].get_visited_id() == visitedId);
        assert(nodes[hash2index.at(NnextL & ~1ULL)].get_level() == 0);
        assert(NnextL < V_HASH_PLUS);
        assert(nextL >= V_HASH_PLUS);
        assert(nextR >= V_HASH_PLUS);

        res_uv_hashes.push_back(std::pair<uint32_t, uint32_t>(NnextL, uint32_t(nextL^1)));
        res_uv_hashes.push_back(std::pair<uint32_t, uint32_t>(NnextL^1, uint32_t(nextR^1)));

        assert(res_uv_hashes.size() == GRAPH_SIZE);
        return res_uv_hashes;
    }

    uint64_t find_prev_node(const std::unordered_map<uint64_t, uint32_t> & hash2index,
                std::vector<Node> & nodes,
                uint64_t node,
                uint32_t visited_id,
                bool left, uint8_t level) {
        assert(hash2index.count(node)==1);
        assert(hash2index.count(node^1)==1);

        uint idx = hash2index.at(node);
        for (auto it=nodes.begin()+idx; it!=nodes.end(); it++) {
            if (it->from != node) {
                assert(false); // not found. Should never happen because path must exist
                return 0;
            }
            assert(hash2index.count(it->to));
            if (nodes[hash2index.at(it->to & ~1ULL)].is_visited_with(visited_id, left, level)) {
                return it->to;
            }
        }
        assert(false);
        return 0;
    }

    void next_step( std::vector<uint64_t> * cur, std::vector<uint64_t> * next, uint32_t visitedId, uint8_t level, bool left,
                const std::unordered_map<uint64_t, uint32_t> & hash2index,
                std::vector<Node> & nodes) {
        next->resize(0);
        for (uint64_t cur_id : *cur) {
            uint idx = hash2index.at(cur_id);
            for (auto it=nodes.begin()+idx; it!=nodes.end(); it++) {
                if (it->from != cur_id)
                    break;
                uint32_t state_idx = hash2index.at(it->to & ~1ULL);
                if (nodes[state_idx].is_visited(visitedId))
                    continue;
                // not visited, we can mark it as visited and continue...
                nodes[state_idx].set_visited(visitedId, left, level);
                next->push_back(it->to ^ 1);
            }
        }
    }


#ifdef ENABLE_TESTS
    // compare masks ant return number of 1s
    uint32_t compare_masks(RowBuilder<EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> & row,
        const std::vector<bool> & mask) {

        const std::vector<uint64_t> & row_mask = row.get_mask();
        assert(mask.size() == uint64_t(row_mask.size())*64/MASK_BITS );
        uint64_t ROW_COUNT_MASK = (1 << MASK_BITS) - 1;
        uint32_t ones = 0;
        for (uint64_t i = 0; i < mask.size(); i++) {
            uint64_t row_bit = i * MASK_BITS;
            uint64_t row_byte = row_bit / 64;
            uint64_t row_shift = row_bit % 64;
            uint64_t mask_val = row_mask[row_byte] & (ROW_COUNT_MASK << row_shift);

            if (mask_val==0) {
                assert(mask[i] == false);
            }
            else {
                ones++;
                assert(mask[i] == true);
            }
        }
        return ones;
    }

    void check_trimmed(const std::vector<MTL::Buffer *> & trimmed, const std::vector<bool> & activeNonces) {
        for (MTL::Buffer * buf : trimmed) {
            uint sz = buf->length()/4;
            assert(sz>0);
            const uint32_t * data = (uint32_t *)buf->contents();
            for (uint i = 0; i < sz; i++) {
                uint32_t nonce = data[i];
                if (nonce == UINT_MAX)
                    break;
                assert(nonce < activeNonces.size());
                assert(activeNonces[nonce]==false);
            }
        }
    }

#ifdef ENABLE_TRACE
    void dump_rowmask_and_trimmed(const char * title,
        RowBuilder<EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> &row ,
        const std::vector<MTL::Buffer *> & trimmed) {

        std::cout << title << std::endl;

        // Printing mask
        std::vector<uint64_t> & mask = row.get_mask();
        uint address = 0;
        uint printed = 0;
        std::cout << "MASK values:" << std::endl;
        for (uint64_t m : mask) {
            // Convert number to hex with leading zeros
            std::ostringstream oss;
            oss << std::hex << std::uppercase << std::setw(16) << std::setfill('0') << m;

            // Insert separators every 4 digits
            std::string hex_str = oss.str();
            std::string formatted;
            for (size_t i = 0; i < hex_str.size(); ++i) {
                if (i > 0 && (i % 4 == 0)) {
                    formatted += "-";
                }
                formatted += hex_str[i];
            }

            if (printed++ % 4 ==0) {
                std::cout << std::endl << std::hex << std::uppercase << std::setw(0) << address << "   ";
            }

            std::cout << formatted << "  ";
            address += 16;
        }
        std::cout << std::endl << std::endl << "Trimmed data:" << std::endl;
        uint buffer_idx = 0;
        for (auto buf : trimmed) {
            uint sz = buf->length()/4;
            assert(sz>0);
            const uint32_t * data = (uint32_t *)buf->contents();
            std::cout << "Buffer #" << (buffer_idx++) << "  size: " << sz << std::endl;
            for (uint i = 0; i < sz; i++) {
                uint32_t nonce = data[i];
                if (nonce==UINT_MAX)
                    break;

                uint32_t nonceR = nonce*2;
                uint32_t HASH_MASK = uint32_t((uint64_t(1) << EDGE_BITS)-1);
                uint32_t hash0 = sip_hash(hash_v, nonceR+0) & HASH_MASK;
                uint32_t hash1 = sip_hash(hash_v, nonceR+1) & HASH_MASK;
                std::cout << std::hex << std::setw(0) << nonce << "(" << nonceR << ")  H:" << hash0 << "-" << hash1 << "    ";
                if (i%16==15) {
                    std::cout << std::endl;
                }
            }
        }
    }
#endif
#endif

#ifdef ENABLE_TRACE
    void test_expected(const std::vector<uint32_t> & expected_u, const std::vector<uint32_t> & expected_v,
                const RowBuilder<EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> & U,
                const RowBuilder<EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS, MASK_BITS, BIT_PACKER_INDEX_INTERVAL> & V,
                const TestLeanMiner<EDGE_BITS> & lean_miner) {

        // let's apply known cycle
        for (uint ex_u : expected_u) {
            assert( U.has_mask(ex_u) );
            assert( lean_miner.getU().at(ex_u));
        }
        for (uint ex_v : expected_v) {
            assert( V.has_mask(ex_v) );
            assert( lean_miner.getV().at(ex_v));
        }
    }
#endif
};

#endif //CUCKATOO_H
