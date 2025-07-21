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

#include "blake.h"
#include "events_tracker.h"
#include "thread_pool.h"
#include "metal.h"
#include "cycle_finder.h"
#include "cycle.h"

#ifdef ENABLE_TESTS
#include "../tests/test_sorting.h"
#endif

#include "metal_context.h"
#include "metal_context_phase1.h"
#include "metal_context_phase2.h"
#include "metal_context_discover.h"
#include "sip_hash.h"
#include <iostream>

#ifdef ENABLE_TRACE
struct SrcData {
    uint32_t nonce;
    uint32_t hashU;
    uint32_t hashV;
};
#endif

// Solution is expected to be proceesed throgn multiple pipeline steps. Here is how we can control that.
struct CuckatooSolution {
    int height;
    int jobId;
    uint64_t nonce;
    uint64_t hash_v[4];

    int current_stage = -1;

    FindCycleStartResult phase1_result;
    std::vector<std::vector<EdgeExt>> solutions;
    std::vector<Cycle> cycles;

    CuckatooSolution(int _height, int _jobId, uint64_t _nonce, uint64_t _hash_v[4], int stage ) :
            height(_height), jobId(_jobId), nonce(_nonce), current_stage(stage) {
        hash_v[0] = _hash_v[0];
        hash_v[1] = _hash_v[1];
        hash_v[2] = _hash_v[2];
        hash_v[3] = _hash_v[3];
    }
};

// How cuckatoo solution can be found.
// EDGE_BITS max value is 31 because high bit us used to track U/V type of the node.
// It i sneeded because of collapsing feature.
template <uint8_t EDGE_BITS, uint GRAPH_SIZE>
class CuckatooSolver {
private:
    // Metal related functionality
    MetalOps metal;

    // Metal context for a stage1 - the primary longest stage
    MetalContextPhase1 context_phase1;
    // Metal context for stage 2 - multiple short steps with some waiting to finish calculations
    MetalContextPhase2 context_phase2;
    // Metal context for hash to nonces mapping. It is used when cycle is found, once in about 42 calculation cycles.
    MetalContextDiscover context_discover;

public:
    // Thread group 256 should be optimal for M4
    CuckatooSolver() : metal(4, 1024UL * 1024UL * 1024UL * 11, // 11 Gb of mem is minimum.
        EDGE_BITS, GRAPH_SIZE,  64, 16, // 64 buckets, 16 of 8b is maximum that fit
        16777213U, // Primary number for 4Mb mask
        13312, // must fit 4k page. item size: 8 byte
        6*1024*1024, 5) // 6M here means 6*8=48Mb of stash data
    {
        // Init metal and contextes. Note, it is expected that phase1 can be processed in one thread, phase2 in another one.
        metal.init();
        context_phase1.init_phase1(metal);
        context_phase2.init_phase2(metal);
        // context_discover has lazy init
    }

    virtual ~CuckatooSolver() {
        // Release all the buffers and metal objects
        context_phase1.release();
        context_phase2.release();
        context_discover.release();

        metal.release();
    }

    /**
     * Stage1 processing routine. Returns back the buffer with edges fo further processing.
     * Note, every next call is single threaded because the same context is used!
     * Caller must release resulting FindCycleStartResult data
     * @param hash_v  - hashes to process.
     * @return Data for the phase2 processing. Caller must release the data form FindCycleStartResult
     */
    FindCycleStartResult phase1(uint64_t hash_v[4]) {
        context_phase1.reset_phase1(metal);
        FindCycleStartResult res = metal.find_cycles_start(context_phase1, hash_v);
        return res;
    }

    /**
     * Stage 2 processing.
     * @param task_context - context of this processing solution
     * @return true if any solution candidate is found (can be false alarm). False - nothing is there.
     */
    bool phase2(CuckatooSolution & task_context ) {
        context_phase2.reset_phase2( task_context.phase1_result);
        metal.find_cycles_phase2(context_phase2, task_context.phase1_result.lastPassU);

        std::vector<Edge> nodes;
        read_node_data( context_phase2,
                            nodes,
                            metal.BUCKETS_NUM,
                            metal.EDGE_PER_BUCKET,
                            metal.PHASE1_EDGE_PER_BUCKET);

        // Discover cycles in resulting data
        if (!nodes.empty()) {
            discover_cycles( nodes,
                             metal.CYCLE_LEN,
                             task_context.solutions );
        }

        // Enrich found cycles with related collapsed data.
        // Need it to minimize the number of steps for enrich_cycles
        if (!task_context.solutions.empty()) {
            ::enrich_cycles(task_context.solutions,
                 metal.COLLAPSE_STASH_COUNTER,
                context_phase2.get_collapse_stash(), metal.COLLAPSE_STASH_SIZE);
        }
        return !task_context.solutions.empty();
    }

    /**
     * Solution has hashes that we need convert the nonce.
     * Note, this step is GPU time consuming, BUT it supposed to be called once in 42 phase1/phase2 calls.
     * @param task_context - context of this processing solution
     */
    void enrich_cycles(CuckatooSolution & task_context) {
        find_nonces(task_context.solutions, metal, context_discover, task_context.hash_v, task_context.cycles);
    }

    /**
     * Testing method for performance build. This method run sequentually all 3 phases and prints resulting metrics
     * @param hash_v - hash that defines solution
     */
    void test_build_all(uint64_t hash_v[4]) {

#ifdef NDEBUG
        // Need few iteraiton because for the first one metrics is incorrect. Takes time for Metal to init buffers
        for (int i=0;i<10;i++)
#endif
        {
            REPORT_RESET;
            REPORT_EVENT("Start","", "main", "");

            CuckatooSolution solution(100, 1, 0, hash_v, -1);

            // Calling all 3 stages
            solution.phase1_result = phase1(hash_v);

            phase2(solution);
            if (!solution.solutions.empty()) {
                enrich_cycles(solution);
                assert(!solution.cycles.empty());
            }
            REPORT_EVENT("Find_Cycles","Start", "main", "");
            REPORT_GENERATE
        }
    }

};

#endif //CUCKATOO_H
