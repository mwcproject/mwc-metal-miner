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

#ifndef METAL_CONTEXT_PHASE2_H
#define METAL_CONTEXT_PHASE2_H

#include "metal_context.h"
#include "features.h"

struct FindCycleStartResult;

// Phase 2 context
class MetalContextPhase2 : public MetalContext {
private:
    // Phase 1 data. Can be represented as 64 buckets
    MTL::Buffer* phase1_results = nullptr;
    // Callapse data.
    MTL::Buffer* collapse_stash = nullptr;

    // device private data for a primary stage.
    MTL::Buffer* primary_data_private = nullptr;
    uint32_t  mask_buffer_size = 0;
    uint32_t  mask_buffer_offset = 0;

    uint32_t  collapse_edges_size = 0;
    uint32_t  collapse_edges_offset = 0;
#ifdef TRACK_COLLAPSED
    uint32_t  collapse_edges_counts_size = 0;
    uint32_t  collapse_edges_counts_offset = 0;
#endif

    // device shared data for a primary stage.
    MTL::Buffer* primary_data_shared = nullptr;

    uint32_t  metrics_size = 0;
    uint32_t  metrics_offset = 0;
    uint32_t trim_params_size = 0; // size per one item. Total BUCKETS_NUM * 2
    uint32_t trim_params_offset = 0;

    uint32_t BUCKETS_NUM = 0;
    uint32_t PHASE1_EDGE_PER_BUCKET;
public:

public:
    virtual ~MetalContextPhase2() {
        release();
    }

    void init_phase2(const MetalOps & metal);
    void reset_phase2(FindCycleStartResult & phase1_res);
    void release();

    virtual MemPool * get_mem_pool() override;
    virtual MTL::Buffer* get_buffer0() override;
    virtual MTL::Buffer* get_buffer1() override;

    virtual Step1Params * get_st1_nonce28_params() override;
    virtual MTL::Buffer* get_st1_nonce28_params_buffer() override;
    virtual uint32_t get_st1_nonce28_params_offset() override;

    virtual Step2Params * get_st2_params(uint32_t bucket_idx) override;
    virtual MTL::Buffer* get_st2_params_buffer() override;
    virtual uint32_t get_st2_params_offset(uint32_t bucket_idx) const override;
#ifdef USE_METRICS
    virtual MetricsData * get_metrics() override;
    virtual MTL::Buffer*  get_metrics_buffer() override;
    virtual uint32_t  get_metrics_offset() override;
#endif
    virtual MTL::Buffer*  get_trimming_params_buffer() override;
    virtual TrimmingParams * get_trimming_params(bool passU, uint32_t bucket_idx) override;
    virtual uint32_t get_trimming_params_offset(bool passU, uint32_t bucket_idx) override;

    virtual TrimmingParams * get_st2_trimming_params(uint32_t bucket_idx) override;
    virtual uint32_t get_st2_trimming_params_offset(uint32_t bucket_idx) const override;

    virtual void * get_mask() override;
    virtual uint32_t get_mask_size() override;
    virtual MTL::Buffer* get_mask_buffer() override;
    virtual uint32_t get_mask_offset() override;

    virtual void * get_collapse_edges() override;
    virtual uint32_t get_collapse_edges_size() override;
    virtual MTL::Buffer* get_collapse_edges_buffer() override;
    virtual uint32_t get_collapse_edges_offset() override;
#ifdef TRACK_COLLAPSED
    virtual void * get_collapse_edges_counts() override;
    virtual uint32_t get_collapse_edges_counts_size() override;
    virtual MTL::Buffer* get_collapse_edges_counts_buffer() override;
    virtual uint32_t get_collapse_edges_counts_offset() override;
#endif

    virtual void * get_bucket_positions() override;
    virtual uint32_t get_bucket_positions_size() override;
    virtual MTL::Buffer* get_bucket_positions_buffer() override;
    virtual uint32_t get_bucket_positions_offset() override;

    virtual MTL::Buffer* get_collapse_stash_buffer() override;
    virtual MTL::Buffer* take_collapse_stash_buffer() override;
    virtual uint32_t* get_collapse_stash() override;
};


#endif //METAL_CONTEXT_PHASE2_H
