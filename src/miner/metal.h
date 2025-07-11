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

#ifndef METAL_H
#define METAL_H

#include <Metal/Metal.hpp>
#include "mem_pool.h"
#include "features.h"
#include "metal_structs.h"
#include "small_hash_table.h"


#ifdef NDEBUG
#define RES_TYPE  MTL::ResourceStorageModePrivate
#else
#define RES_TYPE  MTL::ResourceStorageModeShared
#endif

#define RELEASE(METAL_OBJ) if (METAL_OBJ) { METAL_OBJ->release(); METAL_OBJ = nullptr;}

/**
 * Metal Finction and pipeline state.
 * Sinse they are allways in pairs, let's keep them here
 */
struct MetalCalcPair {
    MTL::Function* function = nullptr;
    MTL::ComputePipelineState* pipeline = nullptr;

    MetalCalcPair() = default;

    ~MetalCalcPair();
    void init(MTL::Device* device, MTL::Library* library, const char * function_name);
};

/**
 * Phase1 resulting data
 * Note, FindCycleStartResult owns buffers data, that is why only move allowed, no copy
 */
struct FindCycleStartResult {
    MTL::Buffer* edges; // Data with resulting edges
    MTL::Buffer* collapse_stash_buffer; // Stash data of collapsed edges.
    bool lastPassU; // Last U/V trim process.

    FindCycleStartResult() : edges(nullptr), collapse_stash_buffer(nullptr), lastPassU(false) {}
    FindCycleStartResult(MTL::Buffer* _edges, MTL::Buffer* _collapse_stash_buffer, bool _lastPassU) :
            edges(_edges), collapse_stash_buffer(_collapse_stash_buffer), lastPassU(_lastPassU) {
        assert(edges);
        assert(collapse_stash_buffer);
    }
    FindCycleStartResult(const FindCycleStartResult&) = delete;
    FindCycleStartResult& operator=(const FindCycleStartResult&) = delete;

    // Move constructor
    FindCycleStartResult(FindCycleStartResult&& other) noexcept {
        edges = other.edges;
        collapse_stash_buffer = other.collapse_stash_buffer;
        other.edges = nullptr;
        other.collapse_stash_buffer = nullptr;
        lastPassU = other.lastPassU;
    }

    // Move assignment
    FindCycleStartResult& operator=(FindCycleStartResult&& other) noexcept {
        if (this != &other) {
            release(); // release current resources
            edges = other.edges;
            collapse_stash_buffer = other.collapse_stash_buffer;
            other.edges = nullptr;
            other.collapse_stash_buffer = nullptr;
            lastPassU = other.lastPassU;
        }
        return *this;
    }

    ~FindCycleStartResult() {
        release();
    }
private:
    void release() {
        if (edges) {
            edges->release();
        }
        if (collapse_stash_buffer) {
            collapse_stash_buffer->release();
        }
    }
};

class MetalContext;
struct CollapsedData;
struct NonceInfo;
struct MetalContextDiscover;

/**
 * Metal processing routine.
 */
struct MetalOps {
public:
    // Metal device
    MTL::Device* device = nullptr;
    // Metal compiled code
    MTL::Library* sip_hash_library = nullptr;

    // Total number of trimming steps in Phase1. Expected that the data amount is significant.
    uint32_t PRIMARY_TRIM_STEPS;
    // Memory size to use for main buffers
    uint64_t MEM_SIZE;
    // Edge bits for our solution. Should be 31 for C31
    uint32_t EDGE_BITS;
    // Expected length of the cycle
    uint32_t CYCLE_LEN;
    // Number of buckets for trimming. Can't be less than CYCLE_LEN
    uint32_t BUCKETS_NUM;
    // Number of full 8b bukets at stage 1. More is better but it require memory.
    // Other buckets are jusy hash 4b buckets.
    uint32_t BUCKETS_8B_NUM;
    // Largest prime number smoller than  EDGE_PER_BUCKET
    uint32_t EDGE_PER_BUCKET_PREV_PRIME;
    // Number of edges per bucket
    uint32_t EDGE_PER_BUCKET;
    // Number of edges per bucket expected at the resulting stage
    uint32_t PHASE1_EDGE_PER_BUCKET;
    // Size of the stash buffer for collapsed data
    uint32_t COLLAPSE_STASH_SIZE;
    // Minimal counetr value to keep edge in the stash
    uint32_t COLLAPSE_STASH_COUNTER;

    // Metal defined functions that we are using to process the data
    MetalCalcPair f_st1_build_buckets;
    MetalCalcPair f_nonce_to_8b;
    MetalCalcPair f_clean_data;
    MetalCalcPair f_clean_data_FF;
    MetalCalcPair f_build_mask_simple;
    MetalCalcPair f_build_mask;
    MetalCalcPair f_trim_mask;
    MetalCalcPair f_trimmed_to_next_buckets_st2;
    MetalCalcPair f_trimmed_to_next_buckets_st_trim;
    MetalCalcPair f_apply_collapsed_data;
    MetalCalcPair f_compact_zeroes;
    MetalCalcPair f_copy_resulting_data;
    MetalCalcPair f_discover_nonces;

    // Twi queues to process that Phase1/phase2. Expected that secondary_queue will run lightwait tasks, so it will be faster.
    MTL::CommandQueue* primary_queue = nullptr;
    MTL::CommandQueue* secondary_queue = nullptr;


public:
    // Threads number can't be more then 1024 !!!!
    MetalOps(uint32_t PRIMARY_TRIM_STEPS, uint64_t MEM_SIZE, uint32_t EDGE_BITS, uint32_t CYCLE_LEN, uint32_t BUCKETS_NUM, uint32_t BUCKETS_8B_NUM,
        uint32_t EDGE_PER_BUCKET_PREV_PRIME, uint32_t PHASE1_EDGE_PER_BUCKET,
        uint32_t COLLAPSE_STASH_SIZE, uint32_t COLLAPSE_STASH_COUNTER);
    ~MetalOps();

    // Init buffers
    void init();

    // First step of finding the graph. All CPU is expecte dto be loaded
    // Result will be provided as a b8 buffer with a data, about 60k pairs. Next it can be processed with a not a fully loaded GPU
    FindCycleStartResult find_cycles_start(MetalContext & context,
        const uint64_t hash_keys[4]);

    // Second phase processing
    void find_cycles_phase2(MetalContext & context, bool lastPassU);

    // release all the buffers
    void release();

    // Hind nonces for known hashes
    void discover_nonces(MetalContextDiscover & context );

private:
    // Primary Execution steps

    // Step1, the first command. Initial distribution by Buckets
    void call_st1_build_buckets(MetalContext & context, std::vector<MTL::CommandBuffer*> & running_commands,
#ifdef STAGE1_TESTS
                    const std::vector<MemRange> & st1_8B_buckets, const std::vector<MemRange> & st1_1B_buckets,
#endif
                    int emit_event,
                    const uint64_t hash_keys[4]);

    // Trimming and secondary distribution by buckets.
    void call_st2_trim_bucketing(MetalContext &context,
                                           int wait_event, int emit_event,
                                           const uint64_t hash_keys[4],
                                           std::vector<MTL::CommandBuffer *> &running_commands,
                                           std::vector<std::vector<MemRange>> & buckets,
                                           int bucket_idx,
                                           const std::vector<MemRange> & st1_8B_buckets, const std::vector<MemRange> & st1_1B_buckets,
                                           std::string & prev_event);

    // Process trimming step, no moving between buckets, data will be reused
    MTL::CommandBuffer* execute_trimming_steps(MetalContext &context,
                                           MTL::CommandBuffer* command,
                                           bool run_tests,
                                           const std::vector<std::vector<MemRange>> & buckets,
                                           bool passU,
                                           int bucket_idx,
                                           int trim_idx,
                                           std::string & prev_event);

    // Copy phase1 resulting data into the buffer.
    void copy_resulting_data(MetalContext &context, std::vector<MTL::CommandBuffer *> &running_commands,
                MTL::Buffer* resulting_buffer,
                const std::vector<std::vector<MemRange>> & buckets,
                std::string & prev_event);
private:
    // Tier-2 execution steps related to a single kernel command
    MTL::CommandBuffer * call_nonce_to_8b(MTL::CommandBuffer *command, MetalContext &context,
                    const MemRange &in_range, const MemRange & b8_range,
                    std::string & prev_event,
                     const uint64_t hash_keys[4] );

    MTL::CommandBuffer * call_build_mask_simple(MTL::CommandBuffer *command, MetalContext &context,
                uint bucket_idx,
                std::string & prev_event);

    MTL::CommandBuffer * call_trimmed_to_next_buckets_st2(MTL::CommandBuffer *command, MetalContext &context,
                uint bucket_idx,
                std::string & prev_event);

    MTL::CommandBuffer * call_compact_zeroes(MTL::CommandBuffer *command, MetalContext &context,
                uint bucket_idx,
                std::string & prev_event,
                uint32_t trim_param_offset);

    MTL::CommandBuffer * call_apply_collapsed_data(MTL::CommandBuffer *command, MetalContext &context,
                 bool isStep2,
                 int mask_scale_k,
                 bool passU, uint bucket_idx,
                std::string & prev_event);

    MTL::CommandBuffer * call_build_mask(MTL::CommandBuffer *command, MetalContext &context,
                    bool passU, uint32_t bucket_idx,
                    std::string & prev_event);

    MTL::CommandBuffer * call_trimmed_to_next_buckets_st_trim(MTL::CommandBuffer *command, MetalContext &context,
                    bool passU, uint32_t bucket_idx,
                    std::string & prev_event);


private:
#ifndef NDEBUG
    // Tier-3 testing related routine

    // Consystency for all
    void test_all_buckets_state(MetalContext &context, const std::vector<std::vector<MemRange>> & buckets);

    void test_initial_buckets_state(MetalContext &context, const std::vector<MemRange> & in,
            bool stage2, bool passU, int buckets_idx,
            int trim_idx, const std::vector<std::vector<MemRange>> & buckets);

    void test_st2_results(MetalContext &context, const std::vector<MemRange> & in,
            const std::vector<std::vector<MemRange>> & buckets,
            int bucket_idx,
            const uint64_t hash_keys[4],
            std::unordered_map<uint64_t, CollapsedData> & out_collapsed_data,
            std::vector<uint32_t> & out_bucket_positions);

    void test_collapsed_data_results(MetalContext &context, int bucket_idx,
        const std::unordered_map<uint64_t, CollapsedData> & collapsed_data,
        const std::vector<uint32_t> & bucket_positions_pre_collapse,
        const std::vector<std::vector<MemRange>> & buckets);

    void test_trim_results(MetalContext &context,
            const std::vector<std::vector<MemRange>> & buckets,
            bool passU,
            int bucket_idx, int trim_idx,
            const std::vector< std::vector<uint64_t>> & in_data,
            const std::vector<uint32_t> bucket_positions_pre,
            uint32_t mask_scale_k,
            uint32_t mask_mul,
            std::unordered_map<uint64_t, CollapsedData> & out_collapsed_data,
            std::vector<uint32_t> & out_bucket_positions,
            uint32_t collapse_stash_prev_size);

    void test_trim_collapdes_results(MetalContext &context,
            const std::vector<std::vector<MemRange>> & buckets,
            int trim_idx,
            const std::unordered_map<uint64_t, CollapsedData> & collapsed_data,
            const std::vector<uint32_t> & bucket_positions_pre_collapse);
#endif
private:
    void init_st1_build_buckets_params(MetalContext &context,
                const uint64_t hash_keys[4],
                Step1Params* params,
                std::vector<MemRange> & st1_8B_buckets,
                std::vector<MemRange> & st1_1B_buckets);

    std::vector<MemRange> init_build_mask_and_buckets_st2(MetalContext &context,
            const uint64_t hash_keys[4],
            Step2Params* kernel_data,
            uint32_t bucket_idx,
            const MemRange & in,
            std::vector<MemRange> & res_buckets);

    std::vector<MemRange> init_build_mask_and_buckets_st_trim(volatile TrimmingParams* kernel_data,
                bool passU,
                uint32_t bucketIdx,
                const std::vector<std::vector<MemRange>> & buckets,
                uint32_t mask_scale_k,
                uint32_t mask_scale_mul);

    std::vector<MemRange> construct_in_data(const std::vector<std::vector<MemRange>> & buckets,
                            bool passU,
                            int bucket_idx) const;
};



#endif //METAL_H
