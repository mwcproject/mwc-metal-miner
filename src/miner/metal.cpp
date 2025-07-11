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

#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION

#pragma clang diagnostic ignored "-Wc++20-extensions"

#include "metal.h"
#include "mem_pool.h"
#include <iostream>
#include "events_tracker.h"
#include <cmath>
#include <chrono>
#include <thread>
#include "metal_code.h"
#include "metal_context.h"
#include "metal_context_discover.h"
#include "features.h"

#include "../tests/test_st1_build_buckets.h"
#include "../tests/test_nonce_to_8b.h"
#include "../tests/test_build_mask.h"

// 1.02 - have few at 64 buckets.  We can accept such loss
#define BUFFERS_OVERHEAD 1.08

#define TRIMMING_RATIO_ST2 0.6

const uint32_t MASK_MUL_ROW[] = {7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
const uint32_t MASK_MUL_ROW_SIZE = sizeof(MASK_MUL_ROW) / sizeof(MASK_MUL_ROW[0]);

MetalCalcPair::~MetalCalcPair() {
    if (function) {
        function->release();
    }
    if (pipeline) {
        pipeline->release();
    }
}

void MetalCalcPair::init(MTL::Device* device, MTL::Library* library, const char * function_name) {
    assert(function==nullptr);
    assert(pipeline==nullptr);

    function = library->newFunction(NS::String::string(function_name, NS::UTF8StringEncoding));
    if (!function) {
        std::cerr << "Failed to find kernel function " << function_name << std::endl;
        exit(-1);
    }

    // Create Compute Pipeline
    NS::Error* error = nullptr;
    pipeline = device->newComputePipelineState(function, &error);
    if (!pipeline) {
        std::cerr << "Failed to create pipeline for " << function_name << ", Err: " << error->localizedDescription()->utf8String() << std::endl;
        exit(-1);
    }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Init sip_hash calculations
static MTL::Library* init_sip_hash_library(MTL::Device* device,uint32_t EDGE_BITS, uint32_t CYCLE_LEN, uint32_t BUCKETS_NUM,
    uint32_t BUCKETS_8B_NUM, const uint32_t EDGE_PER_BUCKET_PREV_PRIME, const uint32_t PHASE1_EDGE_PER_BUCKET,
    uint32_t COLLAPSE_STASH_SIZE, uint32_t COLLAPSE_STASH_COUNTER) {

	NS::Error* error = nullptr;
	MTL::Library* library = device->newLibrary(
	        NS::String::string(generate_metal_sources(EDGE_BITS, CYCLE_LEN, BUCKETS_NUM, BUCKETS_8B_NUM, EDGE_PER_BUCKET_PREV_PRIME,
	                PHASE1_EDGE_PER_BUCKET, COLLAPSE_STASH_SIZE, COLLAPSE_STASH_COUNTER).c_str(),
			NS::UTF8StringEncoding), nullptr, &error);

	if (!library) {
		std::cerr << "Failed to compile Metal sip_hash: " << error->localizedDescription()->utf8String() << std::endl;
		exit(-1);
	}

	return library;
}


MetalOps::MetalOps(uint32_t _PRIMARY_TRIM_STEPS, uint64_t _MEM_SIZE, uint32_t _EDGE_BITS, uint32_t _CYCLE_LEN, uint32_t _BUCKETS_NUM, uint32_t _BUCKETS_8B_NUM,
    uint32_t _EDGE_PER_BUCKET_PREV_PRIME, uint32_t _PHASE1_EDGE_PER_BUCKET,
    uint32_t _COLLAPSE_STASH_SIZE, uint32_t _COLLAPSE_STASH_COUNTER ) :
    PRIMARY_TRIM_STEPS(_PRIMARY_TRIM_STEPS), MEM_SIZE(_MEM_SIZE), EDGE_BITS(_EDGE_BITS), CYCLE_LEN(_CYCLE_LEN), BUCKETS_NUM(_BUCKETS_NUM), BUCKETS_8B_NUM(_BUCKETS_8B_NUM),
    EDGE_PER_BUCKET_PREV_PRIME(_EDGE_PER_BUCKET_PREV_PRIME), PHASE1_EDGE_PER_BUCKET(_PHASE1_EDGE_PER_BUCKET),
    COLLAPSE_STASH_SIZE(_COLLAPSE_STASH_SIZE), COLLAPSE_STASH_COUNTER(_COLLAPSE_STASH_COUNTER)

{
    assert( BUCKETS_NUM<=128);
    assert((uint32_t(1) << EDGE_BITS) % BUCKETS_NUM == 0);
    EDGE_PER_BUCKET = (uint32_t(1) << EDGE_BITS) / BUCKETS_NUM;
}

MetalOps::~MetalOps() {
    release();
}

void MetalOps::init() {
    assert(primary_queue == nullptr);
    assert(secondary_queue == nullptr);

    device = MTL::CreateSystemDefaultDevice();
    if (!device) {
        std::cerr << "Metal is not supported on this device!" << std::endl;
        exit(-1);
    }

	sip_hash_library = init_sip_hash_library(device, EDGE_BITS, CYCLE_LEN, BUCKETS_NUM, BUCKETS_8B_NUM, EDGE_PER_BUCKET_PREV_PRIME,
	        PHASE1_EDGE_PER_BUCKET, COLLAPSE_STASH_SIZE, COLLAPSE_STASH_COUNTER);

    f_st1_build_buckets.init(device, sip_hash_library, "st1_build_buckets");
    f_nonce_to_8b.init(device, sip_hash_library, "nonce_to_8b");
    f_clean_data.init(device, sip_hash_library, "clean_data");
    f_clean_data_FF.init(device, sip_hash_library, "clean_data_FF");
    f_build_mask_simple.init(device, sip_hash_library, "build_mask_simple");
    f_build_mask.init(device, sip_hash_library, "build_mask");
    f_trim_mask.init(device, sip_hash_library, "trim_mask");
    f_trimmed_to_next_buckets_st2.init(device, sip_hash_library, "trimmed_to_next_buckets_st2");
    f_trimmed_to_next_buckets_st_trim.init(device, sip_hash_library, "trimmed_to_next_buckets_st_trim");
    f_apply_collapsed_data.init(device, sip_hash_library, "apply_collapsed_data");
    f_compact_zeroes.init(device, sip_hash_library, "compact_zeroes");
    f_copy_resulting_data.init(device, sip_hash_library, "copy_resulting_data");
    f_discover_nonces.init(device, sip_hash_library, "discover_nonces");

    // Create Command Queue
    primary_queue = device->newCommandQueue();
    secondary_queue = device->newCommandQueue();
}

FindCycleStartResult MetalOps::find_cycles_start(MetalContext & context,
    const uint64_t hash_keys[4])
{
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();

    assert(MEM_SIZE % 2==0);
    assert((MEM_SIZE/2) % 4096 ==0);

    REPORT_EVENT("Start","", "main", "");

    // Allocate memory
    std::vector<MemRange> st1_build_buckets_ranges;

    assert(primary_queue); // call init_seed_hash first!

    // Starting st1_build_buckets
    std::vector<MemRange> st1_8B_buckets;
    std::vector<MemRange> st1_1B_buckets;

    init_st1_build_buckets_params(context, hash_keys, context.get_st1_nonce28_params(), st1_8B_buckets, st1_1B_buckets);

    std::vector<MTL::CommandBuffer*> running_commands;

    int st1_event = context.generate_next_event();

    // Step1 - initial buckets distribution
    call_st1_build_buckets(context, running_commands,
#ifdef STAGE1_TESTS
                    st1_8B_buckets, st1_1B_buckets,
#endif
                    st1_event,
                    hash_keys);

#ifndef NDEBUG
#ifdef SHOW_TRACING
    std::cout << "Stage1 is done, memory usage: " << context.get_mem_pool()->getMemoryStatus() << std::endl;
#endif
#endif

#ifdef METALL_CALLBACKS_PERF
    std::string perf_prev_event = "st1_build_buckets ends";
#else
    std::string perf_prev_event = "st1_build_buckets";
#endif

    // Expected to be a matrix BUCKETS_NUM x BUCKETS_NUM
    std::vector<std::vector<MemRange>> buckets;
    buckets.resize(BUCKETS_NUM);

    // Creating 1b/8b processing pipelines. Result must be 8b data for the V raw, in 64 buckets
    // 8b processing first because it will allow to release some memory thta will be used for 1b

    int prev_event = st1_event;

    for (int bidx = 0; bidx < BUCKETS_NUM; bidx++) {
#ifndef NDEBUG
#ifdef SHOW_TRACING
            std::cout << "First pass at " << bidx << ", memory usage: " << context.get_mem_pool()->getMemoryStatus() << std::endl;
#endif
#endif
            int next_event = context.generate_next_event();
            call_st2_trim_bucketing(context,
                                    prev_event, next_event,
                                    hash_keys,
                                    running_commands,
                                    buckets, bidx,
                                    st1_8B_buckets, st1_1B_buckets,
                                    perf_prev_event);
            prev_event = next_event;
    }

    REPORT_EVENT("stage2 ends", "stage2 start", "main", "");

//    int mask_scale_k = 1;
    int mask_mul_idx = 0;
    bool passU = true; // trimming started from V side and this valie will be inversed to false

//    uint32_t prev_edge_number = uint32_t((uint64_t(1) << EDGE_BITS)-1);

#ifndef NDEBUG
    assert(buckets.size() == BUCKETS_NUM);
    for (uint32_t i = 0; i < BUCKETS_NUM; i++) {
        assert(buckets[i].size() == BUCKETS_NUM);
    }
#endif

#if !defined(NDEBUG) && defined(STAGE_TRIM_TESTS) && defined(ALL_BUCKETS_DATA_VALIDATION)
    // Let's test all the buckets data
    test_all_buckets_state(context,buckets);
#endif

    // let's init the params first
    for (uint k=0;k<2;k++) {
        for (int bidx = 0; bidx < BUCKETS_NUM; bidx++) {
            mask_mul_idx = (mask_mul_idx+1) % MASK_MUL_ROW_SIZE;
            uint32_t mask_mul = MASK_MUL_ROW[mask_mul_idx];
            bool passU = (k==0);
            init_build_mask_and_buckets_st_trim( context.get_trimming_params(passU, bidx), passU,
                    bidx, buckets,
                    1, mask_mul ); // no scale now

        }
    }

#ifndef NDEBUG
    uint32_t last_buckets_sum = UINT_MAX;
#endif

    // Here we want to have trimming loop....
    for (uint trim_idx = 0; trim_idx < PRIMARY_TRIM_STEPS; trim_idx++) { // Ideally we should count the number of edges and stop when that number will be relatevly small.
        passU = !passU;

#ifndef NDEBUG
        { // show progress in debug model
            uint32_t * positions = (uint32_t *) context.get_bucket_positions();
            // getting sum from diagonal. We need estimate only. The real number is not interesting.
            // Note, we want have here sum for a single bucket, below is it compared with EDGE_PER_BUCKET value
            uint32_t buckets_sum = 0;
            for (uint k=0;k<BUCKETS_NUM;k++) {
                buckets_sum += positions[ k*BUCKETS_NUM +k ];
            }

            // Progress is not going down, we can exit
            if (buckets_sum >= uint(last_buckets_sum * 0.99))
                break;
#ifdef SHOW_TRACING
            std::cout << "Current edges_number: " << (buckets_sum*BUCKETS_NUM) <<
                                "  Progress: " << (double(last_buckets_sum-buckets_sum)/last_buckets_sum*100.0)  << std::endl;
#endif
            last_buckets_sum = buckets_sum;
        }
#endif

        std::string start_name = "trim_start_" + std::to_string(trim_idx);
        REPORT_EVENT( start_name.c_str() , "", "main", "");
        perf_prev_event = start_name;

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Start trimming process
        MTL::CommandBuffer* command = primary_queue->commandBuffer();
#ifndef METAL_DO_WAITS
        command->encodeWait(context.get_event(), prev_event);
#endif
        for (int bidx = 0; bidx < BUCKETS_NUM; bidx++) {
#ifndef NDEBUG
#ifdef SHOW_TRACING
            std::cout << "Trimming " << trim_idx << " at " << bidx << std::endl;
#endif
#endif

            command = execute_trimming_steps(context, command,
#ifdef STAGE_TRIM_TESTS
            true,
#else
            false,
#endif
                                           buckets,
                                           passU,
                                           bidx,
                                           trim_idx,
                                           perf_prev_event);
        }

#ifndef METAL_DO_WAITS
        command->commit();
        command->waitUntilCompleted();
#endif
        command->retain();
        command->release();

        std::string end_name = "trim_end_" + std::to_string(trim_idx);
        REPORT_EVENT( end_name.c_str(), start_name.c_str() , "main", "");
        perf_prev_event = end_name;
    }

    // TODO add final step that copy data into the resulting 256Mb buffer for a secondary step

    // result is edges plus positions
    MTL::Buffer* resulting = device->newBuffer( PHASE1_EDGE_PER_BUCKET * BUCKETS_NUM*BUCKETS_NUM * 8 + BUCKETS_NUM*BUCKETS_NUM*4, MTL::ResourceStorageModeShared );
    {
        copy_resulting_data(context, running_commands,
                resulting,
                buckets,
                perf_prev_event);
    }

#ifndef METAL_DO_WAITS
    // waiting for all commands to finish
    for (MTL::CommandBuffer* cmd : running_commands) {
        cmd->waitUntilCompleted();
        assert(cmd->status() == MTL::CommandBufferStatusCompleted);
        // released because earlier it is expected that retain was already called.
        cmd->release();
    }
#else
    assert(running_commands.empty());
#endif

    REPORT_EVENT("PrimaryFinished","Start", "main", "");

    pool->release();

    return FindCycleStartResult(resulting, context.take_collapse_stash_buffer(), passU);
}

void MetalOps::find_cycles_phase2(MetalContext & context, bool lastPassU) {
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();

    REPORT_EVENT("StartPhase2","", "main", "");
    std::string perf_prev_event = "StartPhase2";

    // Allocate memory
    assert(secondary_queue); // call init_seed_hash first!

    std::vector<std::vector<MemRange>> buckets;
    {
        buckets.resize(BUCKETS_NUM);
        assert((PHASE1_EDGE_PER_BUCKET*8) % MEM_POOL_UNITS==0);
        const uint32_t bkt_sz_page = (PHASE1_EDGE_PER_BUCKET*8) / MEM_POOL_UNITS;
        uint32_t idx = 0;
        for (uint ui=0; ui<BUCKETS_NUM; ui++) {
            std::vector<MemRange> & row = buckets[ui];
            row.resize(BUCKETS_NUM);
            for (uint vi=0; vi<BUCKETS_NUM; vi++) {
                // everyhting in the same buffer 0
                row[vi] = MemRange(0, bkt_sz_page * idx, bkt_sz_page * (idx+1), PHASE1_EDGE_PER_BUCKET*8);
                idx++;
            }
        }
    }
#ifndef NDEBUG
    assert(buckets.size() == BUCKETS_NUM);
    for (uint32_t i = 0; i < BUCKETS_NUM; i++) {
        assert(buckets[i].size() == BUCKETS_NUM);
    }
#endif

    int mask_scale_k = 1;
    int mask_mul_idx = 0;
    bool passU = lastPassU;

#if !defined(NDEBUG) && defined(PHASE2_TESTS) && defined(ALL_BUCKETS_DATA_VALIDATION)
    // Let's test all the buckets data
    test_all_buckets_state(context,buckets);
#endif

    uint32_t last_buckets_sum = UINT_MAX;

    // Here we want to have trimming loop....
    uint32_t phase2_trim_idx = 0;
    for (int trim_idx=10; trim_idx<100; trim_idx++) {
        // first let's check how many edges are already in the buckets
        uint32_t * positions = (uint32_t *) context.get_bucket_positions();
        // getting sum from diagonal. We need estimate only. The real number is not interesting.
        // Note, we want have here sum for a single bucket, below is it compared with EDGE_PER_BUCKET value
        uint32_t buckets_sum = 0;
        for (uint k=0;k<BUCKETS_NUM;k++) {
            buckets_sum += positions[ k*BUCKETS_NUM +k ];
        }

        // Progress is not going down, we can exit
        if (buckets_sum >= uint(last_buckets_sum * 0.95))
            break;

        passU = !passU;

        mask_mul_idx = (mask_mul_idx+1) % MASK_MUL_ROW_SIZE;
        uint32_t mask_mul = MASK_MUL_ROW[mask_mul_idx];

#ifndef NDEBUG
#ifdef SHOW_TRACING
        std::cout << "Current edges_number: " << (buckets_sum*BUCKETS_NUM) << "  mask_mul: " << mask_mul <<
                            "  Progress: " << (double(last_buckets_sum-buckets_sum)/last_buckets_sum*100.0)  << std::endl;
                last_buckets_sum = buckets_sum;
#endif
#endif

        const int MASK_SPARE_VAL = 120;
        if ( buckets_sum>1000 && buckets_sum * mask_scale_k < EDGE_PER_BUCKET/MASK_SPARE_VAL ) {
                    mask_scale_k = int(EDGE_PER_BUCKET/MASK_SPARE_VAL / (buckets_sum) ) + 1;
                    mask_scale_k = std::min(mask_scale_k, 256); // Limit mask by 32kb
#ifndef NDEBUG
#ifdef SHOW_TRACING
                    std::cout << "Updating mask_scale_k. New value: " << mask_scale_k << std::endl;
#endif
#endif
        }

        std::string start_name = "phase2_trim_start_" + std::to_string(phase2_trim_idx);
        REPORT_EVENT( start_name.c_str() , "", "main", "");
        perf_prev_event = start_name;

        // init prams
        for (int bidx = 0; bidx < BUCKETS_NUM; bidx++) {
                init_build_mask_and_buckets_st_trim( context.get_trimming_params(passU, bidx), passU,
                        bidx, buckets,
                        mask_scale_k, mask_mul ); // no scale now
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Start trimming process
        MTL::CommandBuffer* command = secondary_queue->commandBuffer();
        MTL::CommandBuffer* cmd_copy = command;

        for (int bidx = 0; bidx < BUCKETS_NUM; bidx++) {
#ifndef NDEBUG
#ifdef SHOW_TRACING
            std::cout << "Trimming " << trim_idx << " at " << bidx << std::endl;
#endif
#endif
            command = execute_trimming_steps(context, command,
#ifdef PHASE2_TESTS
            true,
#else
            false,
#endif
                                           buckets,
                                           passU,
                                           bidx,
                                           trim_idx,
                                           perf_prev_event);
        }

        // Waiting for all. Here we can afford it, no rush with loading the Q
#ifndef METAL_DO_WAITS
        assert(cmd_copy == command);
        // waiting for all commands to finish
        command->commit();
        command->waitUntilCompleted();
        assert(command->status() == MTL::CommandBufferStatusCompleted);
#endif
        command->retain();
        command->release();
        (void)cmd_copy;

        std::string end_name = "trim_end_" + std::to_string(phase2_trim_idx);
        REPORT_EVENT( end_name.c_str(), start_name.c_str() , "main", "");
        perf_prev_event = end_name;
    }

    REPORT_EVENT("SecondaryFinished","Start", "main", "");

    pool->release();
}

// half hash keys
void MetalOps::discover_nonces(MetalContextDiscover & context ) {
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();

    MTL::CommandBuffer* command = secondary_queue->commandBuffer();
    MTL::ComputeCommandEncoder* encoder = command->computeCommandEncoder();
    encoder->setComputePipelineState(f_discover_nonces.pipeline);

    int buf_idx = 0;
    encoder->setBuffer(context.hh_keys, 0, buf_idx++);
    encoder->setBuffer(context.nonces, 0, buf_idx++);
    encoder->setBuffer(context.params, 0, buf_idx++);

    // TODO MOVE ME INTO PARAMS
    // Note, must be SYNCED UP with metal_code.cpp
    const uint DISCOVER_NONCES_PER_JOB = 1024;

    const uint32_t edge_num = uint32_t(1) << EDGE_BITS;
    assert(edge_num % DISCOVER_NONCES_PER_JOB == 0);

    // Dispatch threads for each buffer (assuming 1 result per buffer)
    encoder->dispatchThreads(MTL::Size(edge_num / DISCOVER_NONCES_PER_JOB, 1, 1), MTL::Size(1, 1, 1));
    encoder->endEncoding();
    encoder->retain();
    encoder->release();

    // Submit task, no waiting
    command->commit();
    command->waitUntilCompleted();
    assert(command->status() == MTL::CommandBufferStatusCompleted);
    command->retain();
    command->release();

    pool->release();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
///

// Step1
void MetalOps::call_st1_build_buckets(MetalContext & context, std::vector<MTL::CommandBuffer*> & running_commands,
#ifdef STAGE1_TESTS
                    const std::vector<MemRange> & st1_8B_buckets, const std::vector<MemRange> & st1_1B_buckets,
#endif
                    int emit_event,
                    const uint64_t hash_keys[4])
{
        MTL::CommandBuffer* command = primary_queue->commandBuffer();
        MTL::ComputeCommandEncoder* encoder = command->computeCommandEncoder();
        encoder->setComputePipelineState(f_st1_build_buckets.pipeline);

        int buf_idx = 0;
        encoder->setBuffer(context.get_buffer0(), 0, buf_idx++);
        encoder->setBuffer(context.get_buffer1(), 0, buf_idx++);
        encoder->setBuffer(context.get_st1_nonce28_params_buffer(),  context.get_st1_nonce28_params_offset(), buf_idx++);
#ifdef USE_METRICS
        encoder->setBuffer(context.get_metrics_buffer(), context.get_metrics_offset(), buf_idx++);
#endif

        // TODO MOVE ME INTO PARAMS
        const uint STEP1_THREADS = 512;

        // Dispatch threads for each buffer (assuming 1 result per buffer)
        encoder->dispatchThreads(MTL::Size(STEP1_THREADS*BUCKETS_NUM, 1, 1), MTL::Size(STEP1_THREADS, 1, 1));
        encoder->endEncoding();
        encoder->retain();
        encoder->release();

#ifdef METALL_CALLBACKS_PERF
        command->addCompletedHandler([=](MTL::CommandBuffer* buffer) {
            if (buffer->status() == MTL::CommandBufferStatusCompleted) {
                REPORT_EVENT("st1_build_buckets", "Start", "main", "");
            }
        });
#endif
        command->encodeSignalEvent(context.get_event(), emit_event);

        // Submit task, no waiting
        command->commit();
        command->retain();
#ifdef METAL_DO_WAITS
        command->waitUntilCompleted();
        assert(command->status() == MTL::CommandBufferStatusCompleted);
        command->release();
#else
        running_commands.push_back(command);
#endif

        // Validation...
        {
            // let's run some tests for the resulting data.
#ifndef NDEBUG
            //uint32_t * buf0 = (uint32_t *) context.buffer0->contents();
            //uint32_t * buf1 = (uint32_t *) context.buffer1->contents();

#ifdef USE_METRICS
            // Can't read by pointer directly, CPU cache might be not in sync.
            MetricsData* metr = context.get_metrics();
            MetricsData mrt_cpy;
            memcpy(&mrt_cpy, metr, sizeof(MetricsData));
            assert(mrt_cpy.skipped_edges==0);
#endif
#ifdef STAGE1_TESTS

            uint32_t found_edges = 0;
            for (int u=0;u<st1_8B_buckets.size(); u++) {
                std::vector<uint64_t> data;
                st1_8B_buckets[u].copy_data_into(context, data);
                found_edges += validate_bucket_data(false, EDGE_BITS, BUCKETS_NUM,
                        u, -1, 0, data );
            }
            assert(found_edges < st1_8B_buckets.size()*EDGE_PER_BUCKET * 1.01);
            assert(found_edges > st1_8B_buckets.size()*EDGE_PER_BUCKET * 0.99);

            validate_empty_collapse_buffer( context.get_collapse_edges(), context.get_collapse_edges_size(),
#ifdef TRACK_COLLAPSED
                context.get_collapse_edges_counts(), context.get_collapse_edges_counts_size(),
#endif
                EDGE_PER_BUCKET);

            // here can be any of 256 resulting buffers. Checking only 1 because of the memory limits
            test_st1_build_buckets(hash_keys,
                        context,
                        context.get_st1_nonce28_params()->bucket_size,
                        st1_8B_buckets,
                        st1_1B_buckets,
                        EDGE_BITS);
#endif
#endif
        }
}

void MetalOps::call_st2_trim_bucketing(MetalContext &context,
                                       int wait_event, int emit_event,
                                       const uint64_t hash_keys[4],
                                       std::vector<MTL::CommandBuffer *> &running_commands,
                                       std::vector<std::vector<MemRange>> & buckets,
                                       int bucket_idx,
                                       const std::vector<MemRange> & st1_8B_buckets, const std::vector<MemRange> & st1_1B_buckets,
                                       std::string & prev_event) {

    MemRange b8_range;

    MTL::CommandBuffer *command = primary_queue->commandBuffer();
    MTL::CommandBuffer *cmd_copy = command;

    // nonce to 8b edges if needed
    if (bucket_idx >= BUCKETS_8B_NUM) {
        // we are at 1b buckets. Let's decode them into the bigger buffer
        const MemRange &in_range = st1_1B_buckets[bucket_idx - BUCKETS_8B_NUM];
        b8_range = context.get_mem_pool()->allocate(in_range.get_length_bytes() * 2);

        command = call_nonce_to_8b(command, context,
                    in_range, b8_range,
                    prev_event,
                    hash_keys);

#ifdef METAL_DO_WAITS
        context.get_mem_pool()->release(st1_1B_buckets[bucket_idx - BUCKETS_8B_NUM], context, nullptr); // Original data can be released, it is a relatevely big chunk
#else
        context.get_mem_pool()->release(st1_1B_buckets[bucket_idx - BUCKETS_8B_NUM], context, command); // Original data can be released, it is a relatevely big chunk
#endif

    }
    else {
        b8_range = st1_8B_buckets[bucket_idx];
    }

    std::vector<MemRange> res_buckets;
    // Let's use half of the buckets. I think we should be good enough at this point from collisions point of view.
    std::vector<MemRange> in =
            init_build_mask_and_buckets_st2(context, hash_keys, context.get_st2_params(bucket_idx), bucket_idx, b8_range, res_buckets); // no scale now

    assert(res_buckets.size() == BUCKETS_NUM);
    assert(buckets[bucket_idx].size() == 0);
    buckets[bucket_idx] = res_buckets;

    init_build_mask_and_buckets_st_trim(context.get_st2_trimming_params(bucket_idx), true, bucket_idx, buckets, 1, 1);

#if !defined(NDEBUG) && defined(STAGE2_TESTS)
    test_initial_buckets_state(context, in,
            true, true, bucket_idx, 0, buckets);
#endif

#ifndef METAL_DO_WAITS
    command->encodeWait(context.get_event(), wait_event);
#endif

    command = call_build_mask_simple(command, context, bucket_idx, prev_event);
    command = call_trimmed_to_next_buckets_st2(command, context, bucket_idx, prev_event);
    command = call_compact_zeroes(command, context, bucket_idx, prev_event, context.get_st2_trimming_params_offset(bucket_idx));

#ifndef NDEBUG
    std::unordered_map<uint64_t, CollapsedData> collapsed_data;
    std::vector<uint32_t> bucket_positions_pre_collapse;
    test_st2_results(context, in,
            buckets,
            bucket_idx,
            hash_keys,
            collapsed_data,
            bucket_positions_pre_collapse);
#endif

    command = call_apply_collapsed_data(command, context,
                true, 1, true,
                bucket_idx, prev_event);

    command->encodeSignalEvent(context.get_event(), emit_event);

#ifndef NDEBUG
    test_collapsed_data_results(context, bucket_idx,
        collapsed_data,
        bucket_positions_pre_collapse,
        buckets);
#endif

    command->retain();
#ifdef METAL_DO_WAITS
    command->release();
    command = nullptr;
#else
    assert(cmd_copy == command);
    command->commit();
    running_commands.push_back(command);
#endif
    (void)cmd_copy;

    // Original data can be released, it is a relatevely big chunk
    context.get_mem_pool()->release(b8_range, context, command);

}

void MetalOps::copy_resulting_data(MetalContext &context, std::vector<MTL::CommandBuffer *> &running_commands,
                MTL::Buffer* resulting_buffer,
                const std::vector<std::vector<MemRange>> & buckets,
                std::string & prev_event) {
    MTL::CommandBuffer *command = primary_queue->commandBuffer();
    MTL::ComputeCommandEncoder* encoder = command->computeCommandEncoder();
    encoder->setComputePipelineState(f_copy_resulting_data.pipeline);

    int buf_idx = 0;
    encoder->setBuffer(context.get_buffer0(), 0, buf_idx++);
    encoder->setBuffer(context.get_buffer1(), 0, buf_idx++);
    // any trimming TrimmingParams does work
    encoder->setBuffer(context.get_trimming_params_buffer(),  context.get_trimming_params_offset(false,0), buf_idx++);
    encoder->setBuffer(context.get_bucket_positions_buffer(),  context.get_bucket_positions_offset(), buf_idx++);
    encoder->setBuffer(resulting_buffer,  0, buf_idx++);
    encoder->setBuffer(resulting_buffer,  BUCKETS_NUM*BUCKETS_NUM*PHASE1_EDGE_PER_BUCKET * 8, buf_idx++);

    // TODO MOVE ME INTO PARAMS
    const uint COPY_THREADS = 256;

    // Dispatch threads for each buffer (assuming 1 result per buffer)
    encoder->dispatchThreads(MTL::Size(COPY_THREADS*BUCKETS_NUM*BUCKETS_NUM, 1, 1),
                        MTL::Size(COPY_THREADS, 1, 1));
    encoder->endEncoding();
    encoder->retain();
    encoder->release();

#if defined(METAL_DO_WAITS) && defined(METALL_CALLBACKS_PERF)
    command->addCompletedHandler([=, &prev_event](MTL::CommandBuffer* buffer) {
        if (buffer->status() == MTL::CommandBufferStatusCompleted) {
            REPORT_EVENT("copy_data", prev_event.c_str(), "main", "");
            prev_event = "copy_data";
        }
    });
#endif

    // Submit task, no waiting
    command->commit();
    command->retain();
#ifdef METAL_DO_WAITS
    command->waitUntilCompleted();
    assert(command->status() == MTL::CommandBufferStatusCompleted);
    command->release();
#else
    running_commands.push_back(command);
#endif

#ifdef METAL_DO_WAITS
    // Validation...
    {
        // let's run some tests for the resulting data.
#ifndef NDEBUG
        uint32_t * bucket_pos = (uint32_t *) context.get_bucket_positions();
        uint32_t * target_bucket_pos = (uint32_t *) (((uint64_t *) resulting_buffer->contents()) + BUCKETS_NUM*BUCKETS_NUM*PHASE1_EDGE_PER_BUCKET);

        for (int k=0; k<BUCKETS_NUM*BUCKETS_NUM; k++) {
            assert(bucket_pos[k] == target_bucket_pos[k]);
            assert(bucket_pos[k]<PHASE1_EDGE_PER_BUCKET);
        }

        assert(buckets.size() == BUCKETS_NUM);
        for (int ui=0; ui<BUCKETS_NUM; ui++) {
            for (int vi=0; vi<BUCKETS_NUM; vi++) {
                uint64_t * data = (uint64_t *) buckets[ui][vi].get_data_ptr(context );
                uint64_t * res_data = ((uint64_t *) resulting_buffer->contents()) + (BUCKETS_NUM*ui + vi) * PHASE1_EDGE_PER_BUCKET;
                for (int w=0; w<PHASE1_EDGE_PER_BUCKET; w++) {
                    assert(data[w] == res_data[w]);
                }
            }
        }
#endif
    }
#endif
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<MemRange> MetalOps::construct_in_data(const std::vector<std::vector<MemRange>> & buckets,
                            bool passU,
                            int bucket_idx) const {
    std::vector<MemRange> in;
    assert(buckets.size()==BUCKETS_NUM);
    assert(buckets[0].size()==BUCKETS_NUM);
    if (passU) {
        in = buckets[bucket_idx];
    }
    else {
        for (int k=0;k<BUCKETS_NUM;k++) {
            assert(buckets[k].size() == BUCKETS_NUM);
            in.push_back(buckets[k][bucket_idx]);
        }
    }
    assert(in.size() == BUCKETS_NUM);

    return in;
}

MTL::CommandBuffer* MetalOps::execute_trimming_steps(MetalContext &context,
                                       MTL::CommandBuffer* command,
                                       bool run_tests,
                                       const std::vector<std::vector<MemRange>> & buckets,
                                       bool passU,
                                       int bucket_idx,
                                       int trim_idx,
                                       std::string & prev_event) {

#if !defined(NDEBUG) && (defined(STAGE_TRIM_TESTS) || defined(PHASE2_TESTS))
    std::vector<uint32_t> bucket_positions_pre;
    std::vector< std::vector<uint64_t>> in_data;
    if (run_tests && trim_idx>=STAGE_TRIM_TESTS_START_STEP) {
        std::vector<MemRange> in = construct_in_data(buckets, passU, bucket_idx);

        test_initial_buckets_state(context, in,
                    false, passU, bucket_idx,
                    trim_idx, buckets);

        bucket_positions_pre.resize(BUCKETS_NUM*BUCKETS_NUM);
        assert(context.get_bucket_positions_size() == BUCKETS_NUM*BUCKETS_NUM*4);
        memcpy( bucket_positions_pre.data(), context.get_bucket_positions(), context.get_bucket_positions_size() );

        in_data.resize(in.size());
        assert(in_data.size() == BUCKETS_NUM);

        for (int s = 0; s < in.size(); s++) {
            in[s].copy_data_into(context, in_data[s]);
        }
    }
    uint32_t collapse_stash_prev_size = context.get_collapse_stash()[0];
#endif

#if !defined(NDEBUG) && (defined(STAGE_TRIM_TESTS) || defined(PHASE2_TESTS)) && defined(ALL_BUCKETS_DATA_VALIDATION)
    // Let's test all the buckets data
    if (run_tests && trim_idx >= STAGE_TRIM_TESTS_START_STEP) {
        test_all_buckets_state(context,buckets);
    }
#endif

    command = call_build_mask(command, context,
                        passU, bucket_idx,
                        prev_event);
    command = call_trimmed_to_next_buckets_st_trim(command, context,
                        passU, bucket_idx,
                       prev_event);
    command = call_compact_zeroes(command, context,
                            bucket_idx,
                            prev_event, context.get_trimming_params_offset(passU, bucket_idx));

#if !defined(NDEBUG) && defined(METAL_DO_WAITS) && defined(TRACK_COLLAPSED)
    uint32_t* stash_data = context.get_collapse_stash();
#ifdef SHOW_TRACING
    if (bucket_idx==BUCKETS_NUM-1) {
        std::cout << "collapse stash size: " << stash_data[0] << std::endl;
    }
#endif
    assert(stash_data[0] < COLLAPSE_STASH_SIZE);
#endif


#if !defined(NDEBUG) && (defined(STAGE_TRIM_TESTS) || defined(PHASE2_TESTS))
    uint32_t mask_scale_k = context.get_trimming_params(passU, bucket_idx)->mask_scale_k;
    uint32_t mask_mul = context.get_trimming_params(passU, bucket_idx)->mask_scale_mul;
    std::unordered_map<uint64_t, CollapsedData> collapsed_data;
    std::vector<uint32_t> bucket_threads_positions_trim;
    if (run_tests && trim_idx>=STAGE_TRIM_TESTS_START_STEP) {
        test_trim_results(context, buckets,
                passU,
                bucket_idx, trim_idx,
                in_data,
                bucket_positions_pre,
                mask_scale_k, mask_mul,
                collapsed_data,
                bucket_threads_positions_trim,
                collapse_stash_prev_size);
    }
#endif

#if !defined(NDEBUG) && (defined(STAGE_TRIM_TESTS) || defined(PHASE2_TESTS)) && defined(ALL_BUCKETS_DATA_VALIDATION)
    if (run_tests && trim_idx >= STAGE_TRIM_TESTS_START_STEP) {
        test_all_buckets_state(context,buckets);
    }
#endif

    command = call_apply_collapsed_data(command, context,
                false,
                1,
                passU,
                bucket_idx, prev_event);

#if !defined(NDEBUG) && defined(ALL_BUCKETS_DATA_VALIDATION)
    if (run_tests && trim_idx >= STAGE_TRIM_TESTS_START_STEP) {
        test_all_buckets_state(context,buckets);
    }
#endif

#if !defined(NDEBUG) && (defined(STAGE_TRIM_TESTS) || defined(PHASE2_TESTS))
    if (run_tests && trim_idx>=STAGE_TRIM_TESTS_START_STEP) {
        test_trim_collapdes_results(context,buckets,
                trim_idx,
                collapsed_data,
                bucket_threads_positions_trim);
    }
#endif

    return command;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MTL::CommandBuffer * MetalOps::call_nonce_to_8b(MTL::CommandBuffer *command, MetalContext &context,
                    const MemRange &in_range, const MemRange & b8_range,
                    std::string & prev_event,
                     const uint64_t hash_keys[4] ) {
        MTL::ComputeCommandEncoder *encoder = command->computeCommandEncoder();
        encoder->setComputePipelineState(f_nonce_to_8b.pipeline);

        int buf_idx = 0;
        encoder->setBuffer(in_range.get_buffer(context), uint64_t(in_range.start) * MEM_POOL_UNITS, buf_idx++);
        encoder->setBuffer(b8_range.get_buffer(context), uint64_t(b8_range.start) * MEM_POOL_UNITS, buf_idx++);
        encoder->setBuffer(context.get_st1_nonce28_params_buffer(), context.get_st1_nonce28_params_offset(), 0, buf_idx++);

        // Dispatch threads for each buffer (assuming 1 result per buffer)
        const uint NONCE_TO_8B_THREADS = 512;

        // Dispatch threads for each buffer (assuming 1 result per buffer)
        encoder->dispatchThreads(MTL::Size(NONCE_TO_8B_THREADS * BUCKETS_NUM, 1, 1),
                                 MTL::Size(NONCE_TO_8B_THREADS, 1, 1));
        encoder->endEncoding();
        encoder->retain();
        encoder->release();

#if defined(METAL_DO_WAITS) && defined(METALL_CALLBACKS_PERF)
        command->addCompletedHandler([=, &prev_event](MTL::CommandBuffer *buffer) {
            if (buffer->status() == MTL::CommandBufferStatusCompleted) {
                // let's run some tests for the resulting data.
                REPORT_EVENT("nonce_to_8b", prev_event.c_str(), "main", "");
                prev_event = std::string("nonce_to_8b");
            }
        });
#endif

#ifdef METAL_DO_WAITS
        command->commit();
        command->waitUntilCompleted();
        assert(command->status() == MTL::CommandBufferStatusCompleted);
        command->retain();
        command->release();

        command = primary_queue->commandBuffer();
#endif

        {
#ifndef NDEBUG
#ifdef STAGE2_TESTS
            test_nonce_to_8b_check(hash_keys, context, in_range, b8_range, EDGE_BITS);
#endif
#endif
        }

        return command;
}

MTL::CommandBuffer * MetalOps::call_build_mask_simple(MTL::CommandBuffer *command, MetalContext &context,
            uint bucket_idx,
            std::string & prev_event)
{
    MTL::ComputeCommandEncoder* encoder = command->computeCommandEncoder();
    encoder->setComputePipelineState(f_build_mask_simple.pipeline);

    int buf_idx = 0;
    encoder->setBuffer(context.get_buffer0(), 0, buf_idx++);
    encoder->setBuffer(context.get_buffer1(), 0, buf_idx++);
    encoder->setBuffer(context.get_mask_buffer(), context.get_mask_offset(), buf_idx++);
    encoder->setBuffer(context.get_st2_params_buffer(), context.get_st2_params_offset(bucket_idx), buf_idx++);

    // TODO MOVE ME INTO THE PARAMS
    const int BUILD_MASK_THREADS_NUM = 256;

    // Dispatch threads for each buffer (assuming 1 result per buffer)
    encoder->dispatchThreads(MTL::Size(BUCKETS_NUM * BUILD_MASK_THREADS_NUM, 1, 1),
                                 MTL::Size(BUILD_MASK_THREADS_NUM, 1, 1));
    encoder->endEncoding();
    encoder->retain();
    encoder->release();

    // trim_mask
    encoder = command->computeCommandEncoder();
    encoder->setComputePipelineState(f_trim_mask.pipeline);

    encoder->setBuffer(context.get_mask_buffer(), context.get_mask_offset(), 0);
    assert(context.get_mask_size()%3 == 0);
    const uint half_mask_bytes = context.get_mask_size()/3;
    encoder->setBuffer(context.get_mask_buffer(), context.get_mask_offset() + half_mask_bytes, 1);

    assert(half_mask_bytes % (4 * 4) == 0);
    int tasks = half_mask_bytes / 4 / 4;

    // Dispatch threads for each buffer (assuming 1 result per buffer)
     encoder->dispatchThreads(MTL::Size(tasks, 1, 1), MTL::Size(1, 1, 1)); // 8192 tasks
     encoder->endEncoding();
     encoder->retain();
     encoder->release();

#if defined(METAL_DO_WAITS) && defined(METALL_CALLBACKS_PERF)
        command->addCompletedHandler([=, &prev_event](MTL::CommandBuffer *buffer) {
            if (buffer->status() == MTL::CommandBufferStatusCompleted) {
                // let's run some tests for the resulting data.
                REPORT_EVENT("build_mask_simple", prev_event.c_str(), "main", "");
                prev_event = std::string("build_mask_simple");
            }
        });
#endif

#ifdef METAL_DO_WAITS
     command->commit();
     command->waitUntilCompleted();
     assert(command->status() == MTL::CommandBufferStatusCompleted);
     command->retain();
     command->release();

     command = primary_queue->commandBuffer();
#endif

    return command;
}

MTL::CommandBuffer * MetalOps::call_trimmed_to_next_buckets_st2(MTL::CommandBuffer *command, MetalContext &context,
            uint bucket_idx,
            std::string & prev_event)
{
    MTL::ComputeCommandEncoder* encoder = command->computeCommandEncoder();
    encoder->setComputePipelineState(f_trimmed_to_next_buckets_st2.pipeline);

    uint32_t buf_idx = 0;
    encoder->setBuffer(context.get_buffer0(), 0, buf_idx++);
    encoder->setBuffer(context.get_buffer1(), 0, buf_idx++);
    encoder->setBuffer(context.get_mask_buffer(), context.get_mask_offset(), buf_idx++);
    encoder->setBuffer(context.get_st2_params_buffer(), context.get_st2_params_offset(bucket_idx), buf_idx++);
    encoder->setBuffer(context.get_collapse_edges_buffer(), context.get_collapse_edges_offset(), buf_idx++);
#ifdef TRACK_COLLAPSED
    encoder->setBuffer(context.get_collapse_edges_counts_buffer(), context.get_collapse_edges_counts_offset(), buf_idx++);
#else
    encoder->setBuffer(nullptr, 0,  buf_idx++);
#endif
    encoder->setBuffer(context.get_bucket_positions_buffer(), context.get_bucket_positions_offset(), buf_idx++);
#ifdef USE_METRICS
    encoder->setBuffer(context.get_metrics_buffer(), context.get_metrics_offset(), buf_idx++);
#endif

    // TODO MOVE ME INTO CONFIGS
    const int TRIM_ST2_TREADS_NUM = 512;

    // Dispatch threads for each buffer (assuming 1 result per buffer)
    encoder->dispatchThreads(MTL::Size(BUCKETS_NUM * TRIM_ST2_TREADS_NUM, 1, 1),
                                 MTL::Size(TRIM_ST2_TREADS_NUM, 1, 1)); // 8192 tasks
    encoder->endEncoding();
    encoder->retain();
    encoder->release();

/*#ifndef NDEBUG
    MTL::BlitCommandEncoder* blit = command->blitCommandEncoder();
    blit->synchronizeResource(context.primary_data_private);
    blit->endEncoding();
#endif
*/

#if defined(METAL_DO_WAITS) && defined(METALL_CALLBACKS_PERF)
        command->addCompletedHandler([=, &prev_event](MTL::CommandBuffer *buffer) {
            if (buffer->status() == MTL::CommandBufferStatusCompleted) {
                // let's run some tests for the resulting data.
                REPORT_EVENT("trimmed_to_next_buckets_st2", prev_event.c_str(), "main", "");
                prev_event = std::string("trimmed_to_next_buckets_st2");
            }
        });
#endif

#ifdef METAL_DO_WAITS
     command->commit();
     command->waitUntilCompleted();
     assert(command->status() == MTL::CommandBufferStatusCompleted);
     command->retain();
     command->release();

     command = primary_queue->commandBuffer();
#endif

    return command;
}

MTL::CommandBuffer * MetalOps::call_compact_zeroes(MTL::CommandBuffer *command, MetalContext &context,
                uint bucket_idx,
                std::string & prev_event,
                uint32_t trim_param_offset) {

    MTL::ComputeCommandEncoder* encoder = command->computeCommandEncoder();
    encoder->setComputePipelineState(f_compact_zeroes.pipeline);

    uint32_t buf_idx = 0;
    encoder->setBuffer(context.get_buffer0(), 0, buf_idx++);
    encoder->setBuffer(context.get_buffer1(), 0, buf_idx++);
    encoder->setBuffer(context.get_trimming_params_buffer(), trim_param_offset, buf_idx++);
    encoder->setBuffer(context.get_bucket_positions_buffer(), context.get_bucket_positions_offset(), buf_idx++);

    // TODO MOVE ME INTO CONFIGS
    const int COMPACT_TREADS_NUM = 1024;

    // Dispatch threads for each buffer (assuming 1 result per buffer)
    encoder->dispatchThreads(MTL::Size(BUCKETS_NUM * COMPACT_TREADS_NUM, 1, 1),
                                 MTL::Size(COMPACT_TREADS_NUM, 1, 1)); // 8192 tasks
    encoder->endEncoding();
    encoder->retain();
    encoder->release();

#if defined(METAL_DO_WAITS) && defined(METALL_CALLBACKS_PERF)
        command->addCompletedHandler([=, &prev_event](MTL::CommandBuffer *buffer) {
            if (buffer->status() == MTL::CommandBufferStatusCompleted) {
                // let's run some tests for the resulting data.
                REPORT_EVENT("compact_zeroes", prev_event.c_str(), "main", "");
                prev_event = std::string("compact_zeroes");
            }
        });
#endif

#ifdef METAL_DO_WAITS
     command->commit();
     command->waitUntilCompleted();
     assert(command->status() == MTL::CommandBufferStatusCompleted);
     command->retain();
     command->release();

     command = primary_queue->commandBuffer();
#endif

    return command;
}

MTL::CommandBuffer * MetalOps::call_apply_collapsed_data(MTL::CommandBuffer *command, MetalContext &context,
                bool isStep2,
                int mask_scale_k,
                bool passU, uint bucket_idx,
                std::string & prev_event) {
    MTL::ComputeCommandEncoder* encoder = command->computeCommandEncoder();
    encoder->setComputePipelineState(f_apply_collapsed_data.pipeline);

    uint32_t buf_idx = 0;
    encoder->setBuffer(context.get_buffer0(), 0, buf_idx++);
    encoder->setBuffer(context.get_buffer1(), 0, buf_idx++);
    encoder->setBuffer(context.get_mask_buffer(), context.get_mask_offset(), buf_idx++);
    if (isStep2)
        encoder->setBuffer(context.get_st2_params_buffer(), context.get_st2_trimming_params_offset(bucket_idx), buf_idx++);
    else
        encoder->setBuffer(context.get_trimming_params_buffer(), context.get_trimming_params_offset(passU, bucket_idx), buf_idx++);

    encoder->setBuffer(context.get_bucket_positions_buffer(), context.get_bucket_positions_offset(), 0, buf_idx++);
    encoder->setBuffer(context.get_collapse_edges_buffer(), context.get_collapse_edges_offset(), buf_idx++);
#ifdef TRACK_COLLAPSED
    encoder->setBuffer(context.get_collapse_edges_counts_buffer(), context.get_collapse_edges_counts_offset(), 0, buf_idx++);
#else
    encoder->setBuffer(nullptr, 0, buf_idx++);
#endif
#ifdef USE_METRICS
    encoder->setBuffer(context.get_metrics_buffer(), context.get_metrics_offset(), buf_idx++);
#endif

    int collapse_tasks_num = (context.get_mask_size()/3/8 + mask_scale_k-1) / mask_scale_k;
    encoder->dispatchThreads(MTL::Size(collapse_tasks_num, 1, 1), MTL::Size(1, 1, 1));
    encoder->endEncoding();
    encoder->retain();
    encoder->release();

#if defined(METAL_DO_WAITS) && defined(METALL_CALLBACKS_PERF)
    command->addCompletedHandler([=, &prev_event](MTL::CommandBuffer *buffer) {
        if (buffer->status() == MTL::CommandBufferStatusCompleted) {
            // let's run some tests for the resulting data.
            REPORT_EVENT("apply_collapsed_data", prev_event.c_str(), "main", "");
            prev_event = std::string("apply_collapsed_data");
        }
    });
#endif

#ifdef METAL_DO_WAITS
    command->commit();
    command->waitUntilCompleted();
    assert(command->status() == MTL::CommandBufferStatusCompleted);
    command->retain();
    command->release();

    command = primary_queue->commandBuffer();
#endif

    return command;
}

MTL::CommandBuffer * MetalOps::call_build_mask(MTL::CommandBuffer *command, MetalContext &context,
                    bool passU, uint32_t bucket_idx,
                    std::string & prev_event) {
    // next building the mask
    MTL::ComputeCommandEncoder* encoder = encoder = command->computeCommandEncoder();
    encoder->setComputePipelineState(f_build_mask.pipeline);

    int buf_idx = 0;
    encoder->setBuffer(context.get_buffer0(), 0, buf_idx++);
    encoder->setBuffer(context.get_buffer1(), 0, buf_idx++);
    encoder->setBuffer(context.get_mask_buffer(), context.get_mask_offset(), 0, buf_idx++);
    encoder->setBuffer(context.get_trimming_params_buffer(), context.get_trimming_params_offset(passU, bucket_idx), buf_idx++);
    encoder->setBuffer(context.get_bucket_positions_buffer(), context.get_bucket_positions_offset(), buf_idx++);

    // TODO MOVE ME INTO THE PARAMS
    const int BUILD_MASK_THREADS_NUM = 256;

    // Dispatch threads for each buffer (assuming 1 result per buffer)
    encoder->dispatchThreads(MTL::Size(BUCKETS_NUM * BUILD_MASK_THREADS_NUM, 1, 1), MTL::Size(BUILD_MASK_THREADS_NUM, 1, 1));
    encoder->endEncoding();
    encoder->retain();
    encoder->release();

    // trim_mask
    encoder = command->computeCommandEncoder();
    encoder->setComputePipelineState(f_trim_mask.pipeline);

    assert(context.get_mask_size()%3 == 0);
    const uint half_mask_bytes = context.get_mask_size()/3;
    encoder->setBuffer(context.get_mask_buffer(), context.get_mask_offset(), 0);
    encoder->setBuffer(context.get_mask_buffer(), context.get_mask_offset() + half_mask_bytes, 1);

    assert(half_mask_bytes % (4 * 4) == 0);
    int tasks = half_mask_bytes / 4 / 4;

    // Dispatch threads for each buffer (assuming 1 result per buffer)
    encoder->dispatchThreads(MTL::Size(tasks, 1, 1), MTL::Size(1, 1, 1)); // 8192 tasks
    encoder->endEncoding();
    encoder->retain();
    encoder->release();


#if defined(METAL_DO_WAITS) && defined(METALL_CALLBACKS_PERF)
    command->addCompletedHandler([=, &prev_event](MTL::CommandBuffer *buffer) {
        if (buffer->status() == MTL::CommandBufferStatusCompleted) {
            // let's run some tests for the resulting data.
            REPORT_EVENT("build_mask", prev_event.c_str(), "main", "");
            prev_event = std::string("build_mask");
        }
    });
#endif

#ifdef METAL_DO_WAITS
    command->commit();
    command->waitUntilCompleted();
    assert(command->status() == MTL::CommandBufferStatusCompleted);
    command->retain();
    command->release();

    command = primary_queue->commandBuffer();
#endif

    return command;
}

MTL::CommandBuffer * MetalOps::call_trimmed_to_next_buckets_st_trim(
            MTL::CommandBuffer *command,
            MetalContext &context,
            bool passU, uint32_t bucket_idx,
            std::string & prev_event)
{
    // And now we can process to the next bucket
    MTL::ComputeCommandEncoder* encoder = command->computeCommandEncoder();
    encoder->setComputePipelineState( f_trimmed_to_next_buckets_st_trim.pipeline);

    uint32_t buf_idx = 0;
    encoder->setBuffer(context.get_buffer0(), 0, buf_idx++);
    encoder->setBuffer(context.get_buffer1(), 0, buf_idx++);
    encoder->setBuffer(context.get_mask_buffer(), context.get_mask_offset(), buf_idx++);
    encoder->setBuffer(context.get_trimming_params_buffer(), context.get_trimming_params_offset(passU, bucket_idx), buf_idx++);
    encoder->setBuffer(context.get_bucket_positions_buffer(), context.get_bucket_positions_offset(), buf_idx++);
    encoder->setBuffer(context.get_collapse_edges_buffer(), context.get_collapse_edges_offset(), buf_idx++);
#ifdef TRACK_COLLAPSED
    encoder->setBuffer(context.get_collapse_edges_counts_buffer(), context.get_collapse_edges_counts_offset(), buf_idx++);
    encoder->setBuffer(context.get_collapse_stash_buffer(), 0, buf_idx++);
#else
    encoder->setBuffer(nullptr, 0, buf_idx++);
    encoder->setBuffer(nullptr, 0, buf_idx++);
#endif
#ifdef USE_METRICS
    encoder->setBuffer(context.get_metrics_buffer(), context.get_metrics_offset(), buf_idx++);
#endif

    // TODO MOVE ME INTO THE PARAMS
    const int TRIM_THREADS_NUM = 512;

    // Dispatch threads for each buffer (assuming 1 result per buffer)
    encoder->dispatchThreads(MTL::Size( BUCKETS_NUM * TRIM_THREADS_NUM, 1, 1), MTL::Size(TRIM_THREADS_NUM, 1, 1));
    encoder->endEncoding();
    encoder->retain();
    encoder->release();


#if defined(METAL_DO_WAITS) && defined(METALL_CALLBACKS_PERF)
    command->addCompletedHandler([=, &prev_event](MTL::CommandBuffer *buffer) {
        if (buffer->status() == MTL::CommandBufferStatusCompleted) {
            // let's run some tests for the resulting data.
            REPORT_EVENT("trimmed_to_next_buckets_st_trim", prev_event.c_str(), "main", "");
            prev_event = std::string("trimmed_to_next_buckets_st_trim");
        }
    });
#endif

#ifdef METAL_DO_WAITS
    command->commit();
    command->waitUntilCompleted();
    assert(command->status() == MTL::CommandBufferStatusCompleted);
    command->retain();
    command->release();

    command = primary_queue->commandBuffer();
#endif

    return command;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///     TESTS
#ifndef NDEBUG
void MetalOps::test_initial_buckets_state(MetalContext &context, const std::vector<MemRange> & in,
            bool stage2, bool passU, int buckets_idx,
            int trim_idx, const std::vector<std::vector<MemRange>> & buckets) {

    validate_empty_mask_buffers(context.get_mask(), context.get_mask_size());

    std::vector<std::vector<uint64_t>> in_data;
    in_data.resize(in.size());

    std::vector<uint32_t> bucket_positions;
    bucket_positions.resize(BUCKETS_NUM*BUCKETS_NUM);
    assert(context.get_bucket_positions_size() == BUCKETS_NUM*BUCKETS_NUM*4);

    memcpy( bucket_positions.data(), context.get_bucket_positions(), context.get_bucket_positions_size() );

    for (int s = 0; s < in.size(); s++) {
        in[s].copy_data_into(context, in_data[s]);

        if (stage2) {
            validate_bucket_data(false, EDGE_BITS, BUCKETS_NUM, buckets_idx, -1, 0, in_data[s]);
        }
        else {
            if (passU)
                        validate_bucket_data(true, EDGE_BITS, BUCKETS_NUM,
                               buckets_idx,s,  bucket_positions[buckets_idx*BUCKETS_NUM + s],
                               in_data[s] );
                    else
                        validate_bucket_data(true, EDGE_BITS, BUCKETS_NUM,
                               s, buckets_idx, bucket_positions[s*BUCKETS_NUM + buckets_idx],
                               in_data[s] );
        }
    }

    validate_empty_collapse_buffer(context.get_collapse_edges(), context.get_collapse_edges_size(),
#ifdef TRACK_COLLAPSED
                                       context.get_collapse_edges_counts(), context.get_collapse_edges_counts_size(),
#endif
                                       EDGE_PER_BUCKET);

    if (!stage2) {
        if (trim_idx>5)
           check_data_duplication(context, buckets);
    }
}

void MetalOps::test_st2_results(MetalContext &context, const std::vector<MemRange> & in,
            const std::vector<std::vector<MemRange>> & buckets,
            int bucket_idx,
            const uint64_t hash_keys[4],
            std::unordered_map<uint64_t, CollapsedData> & out_collapsed_data,
            std::vector<uint32_t> & out_bucket_positions)
{
            // uint32_t *buf0 = (uint32_t *) context.get_buffer0()->contents();

            //uint32_t *maskOdd = (uint32_t *) context.get_mask();
            //uint32_t *maskEven = maskOdd + context.mask_buffer_size/3/4;
            //uint32_t *maskSec = maskEven + context.mask_buffer_size/3/4;

            // Can't read by poiner directly, CPU cache might be not in sync.
#ifdef USE_METRICS
            MetricsData *metr = context.get_metrics();
            MetricsData mrt_cpy;
            memcpy(&mrt_cpy, metr, sizeof(MetricsData));
            assert(mrt_cpy.skipped_edges == 0);
            // uint64_t * data = (uint64_t *) in[56].get_data_ptr(buffers);

#endif
#ifdef STAGE2_TESTS
            std::vector<std::vector<uint64_t>> in_data;
            assert(in.size() == BUCKETS_NUM);
            in_data.resize(in.size());
            for (int e=0; e<BUCKETS_NUM; e++) {
                in[e].copy_data_into(context, in_data[e]);
            }

            uint32_t *thr_pos = (uint32_t *) context.get_bucket_positions();
            assert(context.get_bucket_positions_size() == BUCKETS_NUM * BUCKETS_NUM * 4);
            out_bucket_positions.resize(context.get_bucket_positions_size() / 4);
            memcpy(out_bucket_positions.data(), thr_pos, context.get_bucket_positions_size());

            std::vector<std::vector<uint64_t>> out_data;
            std::vector<MemRange> res_buckets = buckets[bucket_idx];
            assert(res_buckets.size() == BUCKETS_NUM);

            out_data.resize(res_buckets.size());
            for (int s = 0; s < res_buckets.size(); s++) {
                res_buckets[s].copy_data_into(context, out_data[s]);

                validate_bucket_data(true, EDGE_BITS, BUCKETS_NUM, bucket_idx, s,
                                     out_bucket_positions[bucket_idx * BUCKETS_NUM + s], out_data[s]);
            }


            validate_filled_collapse_buffer(context.get_collapse_edges(), context.get_collapse_edges_size(),
#ifdef TRACK_COLLAPSED
                                            context.get_collapse_edges_counts(), context.get_collapse_edges_counts_size(),
#endif
                                            EDGE_PER_BUCKET, EDGE_BITS, CYCLE_LEN, BUCKETS_NUM,
                                            context.get_mask(), context.get_mask_size(),
                                            out_collapsed_data);

            test_build_buckets_st2(context.get_mask(), context.get_mask_size(),
                                bucket_idx, out_bucket_positions, in_data, out_data, EDGE_BITS,
                                BUCKETS_NUM, hash_keys, out_collapsed_data);
#endif
}

void MetalOps::test_collapsed_data_results(MetalContext &context, int bucket_idx,
        const std::unordered_map<uint64_t, CollapsedData> & collapsed_data,
        const std::vector<uint32_t> & bucket_positions_pre_collapse,
        const std::vector<std::vector<MemRange>> & buckets) {

#ifndef NDEBUG
            //uint32_t *buf0 = (uint32_t *) context.buffer0->contents();

#ifdef USE_METRICS
            MetricsData *metr = static_cast<MetricsData *>(context.get_metrics());
            MetricsData mrt_cpy;
            memcpy(&mrt_cpy, metr, sizeof(MetricsData));
            assert(mrt_cpy.skipped_edges == 0);
#endif
#endif
#if !defined(NDEBUG) && defined(STAGE2_TESTS)

            validate_empty_collapse_buffer( context.get_collapse_edges(), context.get_collapse_edges_size(),
#ifdef TRACK_COLLAPSED
                                           context.get_collapse_edges_counts(), context.get_collapse_edges_counts_size(),
#endif
                                           EDGE_PER_BUCKET);

            std::vector<uint32_t> bucket_threads_positions_now;
            uint32_t *thr_pos = (uint32_t *) context.get_bucket_positions();
            assert(context.get_bucket_positions_size() == BUCKETS_NUM * BUCKETS_NUM * 4);
            bucket_threads_positions_now.resize(context.get_bucket_positions_size() / 4);
            memcpy(bucket_threads_positions_now.data(), thr_pos, context.get_bucket_positions_size());

            validate_collapsed_data(bucket_positions_pre_collapse, bucket_threads_positions_now, BUCKETS_NUM,
                                    EDGE_PER_BUCKET, collapsed_data, buckets, context);
#endif
}

void MetalOps::test_all_buckets_state(MetalContext &context,
                const std::vector<std::vector<MemRange>> & buckets)
{
    assert(buckets.size()==BUCKETS_NUM);

    std::vector<uint32_t> bucket_positions(BUCKETS_NUM*BUCKETS_NUM, 0);
    assert(context.get_bucket_positions_size() == bucket_positions.size()*4);
    memcpy(bucket_positions.data(), context.get_bucket_positions(), bucket_positions.size()*4);


                for ( int ui=0; ui<BUCKETS_NUM; ui++) {
                    const std::vector<MemRange> & urow = buckets[ui];
                    assert(urow.size() == BUCKETS_NUM);
                    for (int vi=0; vi<BUCKETS_NUM; vi++) {
                        std::vector<uint64_t> data;
                        urow[vi].copy_data_into(context, data);

                        validate_bucket_data(true, EDGE_BITS, BUCKETS_NUM,
                               ui,vi,
                               bucket_positions[ui * BUCKETS_NUM + vi],
                               data );
                    }
                }

}

void MetalOps::test_trim_results(MetalContext &context,
            const std::vector<std::vector<MemRange>> & buckets,
            bool passU,
            int bucket_idx, int trim_idx,
            const std::vector< std::vector<uint64_t>> & in_data,
            const std::vector<uint32_t> bucket_positions_pre,
            uint32_t mask_scale_k,
            uint32_t mask_mul,
            std::unordered_map<uint64_t, CollapsedData> & out_collapsed_data,
            std::vector<uint32_t> & out_bucket_positions,
            uint32_t collapse_stash_prev_size) {
#ifdef USE_METRICS
    MetricsData* metr = context.get_metrics();
    MetricsData mrt_cpy;
    memcpy(&mrt_cpy, metr, sizeof(MetricsData));
    assert(mrt_cpy.skipped_edges==0);
#endif

    if (bucket_idx==0) {
        std::cout << "Bucket 0 size: " << ((uint32_t *) context.get_bucket_positions())[0] << std::endl;
    }

    uint32_t * thr_pos = (uint32_t *) context.get_bucket_positions();
    assert(context.get_bucket_positions_size() == BUCKETS_NUM*BUCKETS_NUM * 4);
    out_bucket_positions.resize(context.get_bucket_positions_size()/4);
    memcpy( out_bucket_positions.data(), thr_pos, context.get_bucket_positions_size() );

    if (trim_idx >= STAGE_TRIM_TESTS_START_STEP) {
        std::vector<MemRange> in = construct_in_data(buckets, passU, bucket_idx);

        std::vector< std::vector<uint64_t>> res_data;
        {
            res_data.resize( in.size() );
            for (int s = 0; s < in.size(); s++) {
                in[s].copy_data_into(context, res_data[s]);

                if (passU)
                    validate_bucket_data(true, EDGE_BITS, BUCKETS_NUM,
                           bucket_idx,s, out_bucket_positions[bucket_idx * BUCKETS_NUM + s],
                           res_data[s] );
                else
                    validate_bucket_data(true, EDGE_BITS, BUCKETS_NUM,
                           s, bucket_idx, out_bucket_positions[s * BUCKETS_NUM + bucket_idx],
                           res_data[s] );
            }
        }

        validate_filled_collapse_buffer( context.get_collapse_edges(), context.get_collapse_edges_size(),
#ifdef TRACK_COLLAPSED
                                           context.get_collapse_edges_counts(), context.get_collapse_edges_counts_size(),
#endif
                                          EDGE_PER_BUCKET,
                                          EDGE_BITS, CYCLE_LEN, BUCKETS_NUM,
                                          context.get_mask(), context.get_mask_size(),
                                          out_collapsed_data );

        std::unordered_set<uint64_t> migrated_data;
        extract_edge_migrated_data(bucket_positions_pre,
                                out_bucket_positions,
                                EDGE_BITS,
                                BUCKETS_NUM,
                                passU,
                                bucket_idx,
                                buckets,
                                context,
                                migrated_data);

        if (trim_idx>5)
            check_data_duplication(context, buckets);

        std::unordered_map<uint64_t, CollapsedData> collapsed_res;
        test_build_mask_buckets(
                context.get_mask(), context.get_mask_size(),
                bucket_idx,
                out_bucket_positions,
                in_data,
                mask_scale_k, mask_mul,
                res_data,
                EDGE_BITS, BUCKETS_NUM, EDGE_PER_BUCKET_PREV_PRIME,
                passU,
                CYCLE_LEN,
                out_collapsed_data,
                migrated_data,
                collapsed_res,
                COLLAPSE_STASH_COUNTER,
                collapse_stash_prev_size,
                context.get_collapse_stash());

        validate_expected_collapsed(out_collapsed_data, collapsed_res);
    }
}

void MetalOps::test_trim_collapdes_results(MetalContext &context,
            const std::vector<std::vector<MemRange>> & buckets,
            int trim_idx,
            const std::unordered_map<uint64_t, CollapsedData> & collapsed_data,
            const std::vector<uint32_t> & bucket_positions_pre_collapse) {
#if !defined(NDEBUG) && defined(STAGE_TRIM_TESTS)
    // Let's test all the buckets data
    std::vector<uint32_t> bucket_threads_positions_now;
    if (trim_idx >= STAGE_TRIM_TESTS_START_STEP) {
        assert(buckets.size()==BUCKETS_NUM);

        uint32_t * thr_pos = (uint32_t *) context.get_bucket_positions();
        assert(context.get_bucket_positions_size() == BUCKETS_NUM*BUCKETS_NUM * 4);
        bucket_threads_positions_now.resize(context.get_bucket_positions_size()/4);
        memcpy( bucket_threads_positions_now.data(), thr_pos, context.get_bucket_positions_size() );

    }
#endif
    {
#ifndef NDEBUG
#ifdef USE_METRICS
        MetricsData* metr = context.get_metrics();
        MetricsData mrt_cpy;
        memcpy(&mrt_cpy, metr, sizeof(MetricsData));
        assert(mrt_cpy.skipped_edges==0);
#endif
#endif

#if !defined(NDEBUG) && defined(STAGE_TRIM_TESTS)
        if (trim_idx >= STAGE_TRIM_TESTS_START_STEP) {
            validate_empty_collapse_buffer(context.get_collapse_edges(), context.get_collapse_edges_size(),
#ifdef TRACK_COLLAPSED
                            context.get_collapse_edges_counts(), context.get_collapse_edges_counts_size(),
#endif
                            EDGE_PER_BUCKET);

            validate_collapsed_data( bucket_positions_pre_collapse, bucket_threads_positions_now,
                         BUCKETS_NUM, EDGE_PER_BUCKET,
                         collapsed_data,
                         buckets,
                         context);
        }
#endif
    }
}
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void MetalOps::init_st1_build_buckets_params(MetalContext &context,
                    const uint64_t hash_keys[4],
                    Step1Params* params,
                    std::vector<MemRange> & st1_8B_buckets,
                    std::vector<MemRange> & st1_1B_buckets) {
    params->key0 = hash_keys[0];
    params->key1 = hash_keys[1];
    params->key2 = hash_keys[2];
    params->key3 = hash_keys[3];

    // need to set bucket address & block_size.  block_size in bytes for a single thread
    uint edge_per_bucket = (uint32_t(1) << EDGE_BITS) / BUCKETS_NUM;
    uint edge_alingment = MEM_POOL_UNITS/4;
    uint bucket_alloc_edges = uint((edge_per_bucket * BUFFERS_OVERHEAD) + edge_alingment - 1)/edge_alingment*edge_alingment;
    while ((bucket_alloc_edges/edge_alingment) % BUCKETS_NUM != 0) {
        bucket_alloc_edges += edge_alingment;
    }
    assert( bucket_alloc_edges*8 % MEM_POOL_UNITS == 0 );
    assert( (bucket_alloc_edges*8 / MEM_POOL_UNITS) % 64 == 0 );

    assert(bucket_alloc_edges % BUCKETS_NUM == 0);
    assert(bucket_alloc_edges > edge_per_bucket);
    assert(double(bucket_alloc_edges)/double(edge_per_bucket) >= BUFFERS_OVERHEAD);
    assert(double(bucket_alloc_edges)/double(edge_per_bucket) < BUFFERS_OVERHEAD + 0.01);

    params->bucket_size = bucket_alloc_edges;
    for (uint k=0; k<BUCKETS_8B_NUM; k++) {
        st1_8B_buckets.push_back( context.get_mem_pool()->allocate(bucket_alloc_edges*8) );
        const MemRange & mr = st1_8B_buckets.back();
        params->bucket_blocks[k] = mr.to_allocated_block();
    }
    for (uint k=BUCKETS_8B_NUM; k<BUCKETS_NUM; k++) {
        st1_1B_buckets.push_back( context.get_mem_pool()->allocate(bucket_alloc_edges*4) );
        const MemRange & mr = st1_1B_buckets.back();
        params->bucket_blocks[k] = mr.to_allocated_block();
    }
}

static void split_into_chunks(const MemRange & in_bucket,
            uint32_t in_chunks,
            std::vector<MemRange> & res_inputs)
{
    assert(in_chunks>0);
    assert(in_bucket.length() > in_chunks);
    assert(in_bucket.get_length_bytes() % MEM_POOL_UNITS == 0);

    assert( in_bucket.length() % in_chunks == 0);

    int bucket_chunk_num = in_bucket.length() / in_chunks;

    res_inputs.resize(in_chunks);

    for (uint k=0; k<in_chunks; k++) {
        res_inputs[k] = MemRange(in_bucket.buffer, in_bucket.start + bucket_chunk_num*k, in_bucket.start + bucket_chunk_num*(k+1), (bucket_chunk_num)*MEM_POOL_UNITS);
    }
}

std::vector<MemRange> MetalOps::init_build_mask_and_buckets_st2(MetalContext &context,
            const uint64_t hash_keys[4],
            Step2Params* params,
            uint32_t bucket_idx,
            const MemRange & in,
            std::vector<MemRange> & res_buckets) {

    params->key0 = hash_keys[0];
    params->key1 = hash_keys[1];
    params->key2 = hash_keys[2];
    params->key3 = hash_keys[3];

    std::vector<MemRange> inputs;
    split_into_chunks(in, BUCKETS_NUM, inputs);
    uint32_t in_sz = inputs.size();
    assert(in_sz==BUCKETS_NUM);
    assert(in_sz<=128);
    // in bucket size needs to be accumulate, let's do in b8 num.
    for (uint k=0; k<in_sz; k++) {
        params->in_blocks[k] = inputs[k].to_allocated_block();
    }
    params->in_block_size = inputs[0].get_length_bytes()/8;

    uint32_t bucket_size = (uint32_t( (uint32_t(1) << EDGE_BITS) / BUCKETS_NUM / BUCKETS_NUM * TRIMMING_RATIO_ST2 * BUFFERS_OVERHEAD * 8.0) / MEM_POOL_UNITS + 1) * MEM_POOL_UNITS;

    assert(BUCKETS_NUM<=256);
    assert(bucket_idx < BUCKETS_NUM);
    // current bucket:
    // 1 if for UPass
    params->row_idx = bucket_idx;

    assert(bucket_size % MEM_POOL_UNITS == 0);

    res_buckets.resize(BUCKETS_NUM);
    for (uint k=0; k<BUCKETS_NUM; k++) {
        res_buckets[k] = context.get_mem_pool()->allocate(bucket_size);
        params->out_blocks[k] = res_buckets[k].to_allocated_block();
    }
    params->out_block_size = bucket_size / 8;

    return inputs;
}

std::vector<MemRange> MetalOps::init_build_mask_and_buckets_st_trim(volatile TrimmingParams* kernel_data,
            bool passU,
            uint32_t bucket_idx,
            const std::vector<std::vector<MemRange>> & buckets,
            uint32_t mask_scale_k,
            uint32_t mask_scale_mul) {

    std::vector<MemRange> ins;
    assert(buckets.size()==BUCKETS_NUM);
    assert(buckets[0].size()==BUCKETS_NUM);

#ifndef NDEBUG
    { // Let's check that bucket size is the same for all
        uint64_t bkt_len = buckets[0][0].get_length_bytes();
        for ( const auto & row : buckets) {
            for (const auto & b : row) {
                assert( b.get_length_bytes() == bkt_len);
            }
        }
    }
#endif

    int bi = 0;
    for ( const auto & row : buckets) {
        for (const auto & b : row) {
            kernel_data->blocks[bi++] = b.to_allocated_block();
        }
    }

    if (passU) {
        ins = buckets[bucket_idx];
    }
    else {
        ins.resize(BUCKETS_NUM);
        for (uint k=0; k<BUCKETS_NUM; k++) {
            ins[k] = buckets[k][bucket_idx];
        }
    }

    assert(ins.size()==BUCKETS_NUM);

    kernel_data->mask_scale_k = mask_scale_k;
    kernel_data->mask_scale_mul = mask_scale_mul;
    kernel_data->pass_bucket = passU ? 1 : 0;

    assert(BUCKETS_NUM<=256);
    assert(bucket_idx < BUCKETS_NUM);
    // current bucket:
    kernel_data->pass_bucket |= uint32_t(bucket_idx) << 4;

    kernel_data->block_size = buckets[0][0].get_length_bytes()/8;

    return ins;
}


void MetalOps::release() {
    RELEASE(primary_queue);
    RELEASE(secondary_queue);

 	RELEASE(sip_hash_library);
	RELEASE(device);
}



