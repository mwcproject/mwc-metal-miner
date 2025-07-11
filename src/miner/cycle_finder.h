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

#ifndef CYCLE_FINDER_H
#define CYCLE_FINDER_H

#include <vector>
#include <assert.h>
#include <unordered_map>


#define STATE_EXCLUDED 0x01

struct Cycle;
struct MetalOps;
struct MetalContextDiscover;
class MetalContextPhase2;

/**
 * Edge how it was represented at metal
 */
struct Edge {
    uint32_t hash1; // Hashes, both can be U or V. Suported UV/UU/VV
    uint32_t hash2;
    uint8_t  count; // Collapsed counter
    uint8_t  state = 0; // Internal state

    inline Edge(uint32_t _hash1, uint32_t _hash2, uint8_t _count ) : hash1(_hash1), hash2(_hash2), count(_count) {}
    Edge() = default;
    Edge(const Edge &) = default;
    Edge & operator=(const Edge &) = default;

    inline bool can_process() const {return (state & STATE_EXCLUDED) == 0; }
    inline void set_processed() { state |= STATE_EXCLUDED; }

    inline uint32_t get_other_hash(uint32_t hash) {
        assert(hash == hash1 || hash==hash2);
        if (hash==hash1)
            return hash2;
        else
            return hash1;
    }
};

/**
 * Collapse earlier hashes with counters
 */
struct InterimNodes {
    uint32_t hash;
    uint8_t count;

    InterimNodes(uint32_t _hash, uint8_t _count) : hash(_hash), count(_count) {}
    InterimNodes() = default;
    InterimNodes(const InterimNodes &) = default;
    InterimNodes & operator = (const InterimNodes &) = default;
};

/**
 * Edge with discovered hased
 */
struct EdgeExt {
    // two hashes and a count.
    uint32_t hash1;
    uint32_t hash2;
    uint8_t  count;
    // Discovered collapsed hashes that was dropped before.
    std::unordered_map<uint32_t, InterimNodes> interim_hashes;

    EdgeExt(Edge & n) : hash1(n.hash1), hash2(n.hash2), count(n.count) {}
    EdgeExt() = default;
    EdgeExt(const EdgeExt &) = default;
    EdgeExt(EdgeExt &&) = default;
};

/**
 * Read the Edges form the metal data
 * @param phase2context - metal phase2 context with buffers that has data
 * @param resulting_data - Resulting data that was red.
 * @param BUCKETS_NUM  - Buckets Num
 * @param EDGE_PER_BUCKET - Edge per buckets as overall
 * @param PHASE1_EDGE_PER_BUCKET - Edge per buckets at phase1 result.
 */
void read_node_data( MetalContextPhase2 & phase2context,
                    std::vector<Edge> & resulting_data,
                    uint32_t BUCKETS_NUM,
                    uint32_t EDGE_PER_BUCKET,
                    uint32_t PHASE1_EDGE_PER_BUCKET);

/**
 * Discover cycles into the trimming results. Expected that trim result is relatevely small, do DFS search in single
 * thread
 * @param nodes - phase2 resulting nodes
 * @param CYCLE_LEN - Cycle length that is needed for solution (42).
 * @param result_solutions - resulting cycles
 */
void discover_cycles( std::vector<Edge> & nodes,
                         uint32_t CYCLE_LEN,
                         std::vector<std::vector<EdgeExt>> & result_solutions );

/**
 * Enriching Nodes with large counters from collapse_stash data.
 * @param solutions - solutions
 * @param COLLAPSE_STASH_COUNTER - Counter when we start Stashing the collapsed edges.
 * @param collapse_stash - Stashing data from the metal
 * @param COLLAPSE_STASH_SIZE - Size of the stashed data
 */
void enrich_cycles(std::vector<std::vector<EdgeExt>> & solutions,
            uint32_t COLLAPSE_STASH_COUNTER,
            uint32_t * collapse_stash, uint32_t COLLAPSE_STASH_SIZE);

/**
 * Find nonces from the known hashes.
 * @param solutions - Solutions to process
 * @param metal - metal engine to run discovery cycles
 * @param context_discover - Metal contect for that job. One context per a single thread
 * @param v  - hash that was assigned to that task
 * @param result_cycles - Resulting data
 */
void find_nonces(std::vector<std::vector<EdgeExt>> & solutions, MetalOps & metal, MetalContextDiscover & context_discover,
            uint64_t v[4],
            std::vector<Cycle> & result_cycles);

#endif //CYCLE_FINDER_H
