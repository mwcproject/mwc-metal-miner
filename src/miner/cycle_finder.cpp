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

#include "cycle_finder.h"
#include <assert.h>
#include <atomic>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include "cycle.h"
#include "features.h"
#include "metal.h"
#include "metal_context_discover.h"
#include "sip_hash.h"
#include "small_hash_table.h"
#include "utils.h"
#include <iostream>
#include "metal_context_phase2.h"

struct PathNode {
    uint32_t prev_idx;
    uint32_t hash;
    Edge * node;
    uint8_t  len;

    PathNode(uint32_t _prev_idx, uint32_t _hash, Edge * _node, uint8_t _len) :
                prev_idx(_prev_idx), hash(_hash), node(_node), len(_len) {}
    PathNode() = default;
    PathNode(const PathNode &) = default;
};

static std::vector<Edge*> rewind_solutions(const std::vector<PathNode> & paths, uint32_t path_idx, Edge * lastNode) {
    std::vector<Edge*> result;
    result.push_back(lastNode);
    while (path_idx>0) {
        result.push_back( paths[path_idx].node);
        path_idx = paths[path_idx].prev_idx;
    }
    return result;
}

static bool is_visited(const std::vector<PathNode> & paths, uint32_t path_idx, uint32_t next_hash) {
    while (path_idx>0) {
        const PathNode & pn = paths[path_idx];
        if (pn.hash/2 == next_hash/2)
            return true;
        path_idx = paths[path_idx].prev_idx;
    }
    return false;
}

void read_node_data( MetalContextPhase2 & phase2context,
                    std::vector<Edge> & resulting_data,
                    uint32_t BUCKETS_NUM,
                    uint32_t EDGE_PER_BUCKET,
                    uint32_t PHASE1_EDGE_PER_BUCKET) {

    uint32_t * bucket_positions = (uint32_t *) phase2context.get_bucket_positions();

    uint32_t total_edges = 0;
    for (uint32_t i=0; i<BUCKETS_NUM*BUCKETS_NUM; i++) {
        total_edges += bucket_positions[i];
    }
    resulting_data.reserve(total_edges);

    if (total_edges==0)
        return;

    uint64_t * data_buffer = (uint64_t *) phase2context.get_buffer0()->contents();

    const uint32_t Ubit = uint32_t(1) << 31;
    const uint32_t Umask = Ubit - 1;

    for (uint32_t ui=0; ui<BUCKETS_NUM; ui++) {
        uint32_t base_idx = ui*BUCKETS_NUM;
        for (uint32_t vi=0;vi<BUCKETS_NUM; vi++) {
            uint32_t bucket_idx = base_idx + vi;
            uint64_t * bkt_data = data_buffer + bucket_idx*PHASE1_EDGE_PER_BUCKET;
            uint32_t sz = bucket_positions[bucket_idx];
            assert(sz<PHASE1_EDGE_PER_BUCKET);
            for (uint i=0;i<sz; i++) {
                uint64_t dt = bkt_data[i];
                if (!dt)
                    continue;
                uint32_t hash1 = uint32_t(dt);
                uint32_t hash2 = uint32_t(dt >> 32);
                bool u1 = (hash1 & Ubit) != 0;
                bool u2 = (hash2 & Ubit) != 0;

                uint8_t collapsed = 0;
#ifdef TRACK_COLLAPSED
                if (u1==false && u2==false) {
                    collapsed = hash2 / EDGE_PER_BUCKET;
                    // it is V type, use V bucket to recnstruct
                    hash2 = hash2 % EDGE_PER_BUCKET + EDGE_PER_BUCKET*vi;
                }
                else {
                    collapsed = (hash1 & Umask) / EDGE_PER_BUCKET;
                    if (u1) {
                        hash1 = (((hash1 & Umask) % EDGE_PER_BUCKET) + EDGE_PER_BUCKET*ui) | Ubit;
                    }
                    else {
                        hash1 = (hash1 % EDGE_PER_BUCKET) + EDGE_PER_BUCKET*vi;
                    }
                }
#endif

                resulting_data.push_back( Edge(hash1, hash2, collapsed ) );
            }
        }
    }
}


void discover_cycles( std::vector<Edge> & nodes,
                         uint32_t CYCLE_LEN,
                         std::vector<std::vector<EdgeExt>> & resul_solutions ) {
    resul_solutions.resize(0);

    const uint32_t total_nodes = nodes.size();

    // Hash 2 node indexes
    std::unordered_multimap<uint32_t, Edge*> hash2node;
    hash2node.max_load_factor(0.75);
    hash2node.rehash( next_prime(uint32_t(std::max(1000U,total_nodes)*2 / 0.7)) );

    for (Edge & nd : nodes) {
        hash2node.emplace(nd.hash1, &nd);
        hash2node.emplace(nd.hash2, &nd);
    }

    std::vector<PathNode> paths;

    for (uint32_t nIdx=0; nIdx<total_nodes; nIdx++ ) {
        Edge & root = nodes[nIdx];
        assert(root.can_process());

        // Single node cycle (not expected, but possible)
        if (root.count==CYCLE_LEN-1) {
            if (root.hash1/2==root.hash2/2) {
                resul_solutions.push_back( {root} );
            }
            else {
                //assert(false);
            }
            root.set_processed();
            continue;
        }

        // Processing
        // Let's do DFS search. Going from hash1, should reash hash2
        uint32_t start_hash = root.hash1;
        uint32_t end_hash = root.hash2 ^ 0x01;

        paths.resize(0);
        paths.emplace_back( UINT_MAX, start_hash, &root, root.count+1 );
        uint32_t path_idx = 0;

        while(path_idx < paths.size()) {
            const PathNode & cur_node = paths[path_idx];
            uint32_t pair_hash = cur_node.hash ^ 0x01;
            auto next_nodes = hash2node.equal_range(pair_hash);
            for (auto nn = next_nodes.first; nn != next_nodes.second; ++nn) {
                Edge * nNode = nn->second;
                if (!nNode->can_process())
                    continue;
                uint8_t nextCount = cur_node.len + nNode->count + 1;
                if (nextCount>CYCLE_LEN)
                    continue;

                uint32_t next_hash = nNode->get_other_hash(pair_hash);
                if (nextCount==CYCLE_LEN) {
                    if (next_hash == end_hash) {
                        // Solution is found!!!
                        std::vector<Edge*> solution = rewind_solutions(paths, path_idx, nNode);
                        std::vector<EdgeExt> & res = resul_solutions.emplace_back();
                        for (const auto s : solution )
                            res.push_back(EdgeExt(*s));
                        break;
                    }
                    continue;
                }

                // Can going forward.
                if (next_hash==end_hash) {
                    continue; // dead end, too short cycle
                }

                if (is_visited(paths, path_idx, next_hash)) {
                    continue;
                }

                // Not visited, short path. Can add to next review
                paths.emplace_back( path_idx, next_hash, nNode, nextCount );
            }
            path_idx++;
        }
        root.set_processed();
    }
}

const uint32_t STASH_SHIFT = 32-5;
const uint32_t STASH_MASK = (uint32_t(1) << STASH_SHIFT) - 1;


static void generate_stash_items(std::vector<std::vector<EdgeExt>> & solutions, uint32_t COLLAPSE_STASH_COUNTER,
            std::unordered_set<uint32_t> & procesed,
            std::vector<uint32_t> & resulted_keys) {
    resulted_keys.resize(0);
    for (const auto & sol : solutions) {
        for (const auto & n : sol) {
            if (n.count>COLLAPSE_STASH_COUNTER) {
                if (procesed.insert(n.hash1 & STASH_MASK).second)
                    resulted_keys.push_back( n.hash1 & STASH_MASK );
                if (procesed.insert(n.hash2 & STASH_MASK).second)
                    resulted_keys.push_back( n.hash2 & STASH_MASK );

                for ( auto interim : n.interim_hashes ) {
                    if (procesed.insert(interim.second.hash & STASH_MASK).second)
                        resulted_keys.push_back( interim.second.hash & STASH_MASK );

                    if (procesed.insert((interim.second.hash ^ 0x01) & STASH_MASK).second)
                        resulted_keys.push_back( (interim.second.hash ^ 0x01) & STASH_MASK );
                }
            }
        }
    }
}

inline static void append_hashes(const std::vector<std::pair<uint32_t, uint32_t>> & found_hashes, const InterimNodes & interim,
                    std::unordered_map<uint32_t, InterimNodes> & interim_hashes) {
    for ( const std::pair<uint32_t, uint32_t> & counter_hash : found_hashes ) {
        uint32_t counter = std::get<0>(counter_hash);
        if (counter < interim.count) {
            uint32_t hash = std::get<1>(counter_hash);
            interim_hashes[hash] = InterimNodes(hash, counter);
        }
    }

}

// return true if something was found
static void enrich_node(const SmallHashTable<std::vector<std::pair<uint32_t, uint32_t>> > & found_paired,
            EdgeExt & node) {

    std::vector<uint32_t> keys;
    keys.reserve(node.interim_hashes.size());
    for (const auto & ni : node.interim_hashes)
        keys.push_back(ni.first);

    for ( auto k : keys ) {
        InterimNodes interim = node.interim_hashes.at(k);
        if (found_paired.contains(interim.hash & STASH_MASK)) {
            append_hashes( found_paired.get(interim.hash & STASH_MASK), interim, node.interim_hashes );
        }
        if (found_paired.contains((interim.hash ^ 0x01) & STASH_MASK)) {
            append_hashes( found_paired.get((interim.hash ^ 0x01) & STASH_MASK), interim, node.interim_hashes );
        }
    }

    if ( found_paired.contains(node.hash1 & STASH_MASK) ) {
        append_hashes( found_paired.get(node.hash1 & STASH_MASK), InterimNodes(node.hash1, node.count), node.interim_hashes );
    }
    if ( found_paired.contains(node.hash2 & STASH_MASK) ) {
        append_hashes( found_paired.get(node.hash2 & STASH_MASK), InterimNodes(node.hash2, node.count), node.interim_hashes );
    }
}

// Enriching Nodes with large conters from collapse_stash data.
void enrich_cycles(std::vector<std::vector<EdgeExt>> & solutions,
            uint32_t COLLAPSE_STASH_COUNTER,
            uint32_t * collapse_stash, uint32_t COLLAPSE_STASH_SIZE) {

    std::vector<uint32_t> stash_items_to_find;
    std::unordered_set<uint32_t> procesed;
    SmallHashTable<std::vector<std::pair<uint32_t, uint32_t>> >  found_paired;

    const std::vector<std::pair<uint32_t, uint32_t>> empty;
    std::vector<uint32_t> to_process;
    while(true) {
        generate_stash_items(solutions, COLLAPSE_STASH_COUNTER, procesed, stash_items_to_find);
        to_process.insert(to_process.end(), stash_items_to_find.begin(), stash_items_to_find.end());
        if (to_process.empty())
            break;

        to_process = found_paired.reset( to_process, stash_items_to_find.size()*2, 0xFFFFFFFF, empty );

        uint32_t collapsed_size = std::min(COLLAPSE_STASH_SIZE, collapse_stash[0]);
        for (uint32_t k=0; k<collapsed_size; k++) {
            uint32_t key = collapse_stash[1+k*2] & STASH_MASK;
            if (found_paired.contains(key)) {
                found_paired.get(key).push_back( std::pair<uint32_t, uint32_t>( (collapse_stash[1+k*2]>>STASH_SHIFT) + COLLAPSE_STASH_COUNTER, collapse_stash[1+k*2+1] ) );
            }
        }

        // Checking if something is found, so we can expand our search
        for (auto & sol : solutions) {
            for (auto & n : sol) {
                enrich_node(found_paired, n);
            }
        }
    }

}

void find_nonces(std::vector<std::vector<EdgeExt>> & solutions, MetalOps & metal, MetalContextDiscover & context_discover,
            uint64_t v[4], std::vector<Cycle> & result_cycles) {
    result_cycles.reserve(solutions.size());

    for (const auto & sol : solutions) {
        result_cycles.emplace_back(sol);
    }

    assert(!result_cycles.empty());

    std::vector<uint32_t> rejected_u;
    std::vector<uint32_t> rejected_v;

    // Let's start from U. In general it doesn't matter start form U or V.
    bool upass = false;
    for (uint u=0; u<metal.COLLAPSE_STASH_COUNTER*2+1+2; u++, upass = !upass) { //
        std::vector<uint32_t> & rejected = upass ? rejected_u : rejected_v;

        std::unordered_set<uint32_t> hashes(rejected.begin(), rejected.end());

        for (auto & c : result_cycles) {
            auto & hh = c.get_hashes_to_discover(upass);
            hashes.insert( hh.begin(), hh.end() );
            hh.clear();
        }

        if (hashes.empty())
            continue;

        SmallHashTable<bool> hash_table;
        rejected = hash_table.reset( std::vector<uint32_t>(hashes.begin(), hashes.end()), hashes.size()*3, 0xFFFFFFFF, false );

        context_discover.init(&metal, hash_table, upass, v);
        metal.discover_nonces( context_discover );
        std::vector<NonceInfo> discovered_nonces;
        context_discover.read_results(discovered_nonces);

#ifndef NDEBUG
        uint32_t hash_mask = (uint32_t(1) << metal.EDGE_BITS)-1;
        for (const NonceInfo & n : discovered_nonces) {
            uint32_t uh = sip_hash( v, uint64_t(n.nonce)*2 ) & hash_mask;
            uint32_t vh = sip_hash( v, uint64_t(n.nonce)*2+1 ) & hash_mask;
            assert( n.uhash == uh );
            assert( n.vhash == vh );
        }
#endif

        bool done = true;
        for (auto & c : result_cycles) {
            c.add_nonces(discovered_nonces, upass);
            if (!c.detect_cycle(metal.CYCLE_LEN)) {
                done = false;
            }
        }

        if (done)
            break;
    }
    //
}
