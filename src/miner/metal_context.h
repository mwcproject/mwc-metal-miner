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

#ifndef METAL_PRIMARY_CONTEXT_H
#define METAL_PRIMARY_CONTEXT_H

#include <Metal/Metal.hpp>
#include "features.h"
#include "metal_structs.h"

struct MetalOps;
class MemPool;

// Context for the metal operations.
// Needed because phase1 and phase2 operates on same routine but very diffrent data
class MetalContext {
public:
    virtual ~MetalContext() {}

    virtual MemPool * get_mem_pool() = 0;
    virtual MTL::Buffer* get_buffer0() = 0;
    virtual MTL::Buffer* get_buffer1() = 0;

    virtual Step1Params * get_st1_nonce28_params() = 0;
    virtual MTL::Buffer* get_st1_nonce28_params_buffer() = 0;
    virtual uint32_t get_st1_nonce28_params_offset() = 0;

    virtual Step2Params * get_st2_params(uint32_t bucket_idx) = 0;
    virtual MTL::Buffer* get_st2_params_buffer() = 0;
    virtual uint32_t get_st2_params_offset(uint32_t bucket_idx) const = 0;
#ifdef USE_METRICS
    virtual MetricsData * get_metrics() = 0;
    virtual MTL::Buffer*  get_metrics_buffer() = 0;
    virtual uint32_t  get_metrics_offset() = 0;
#endif
    virtual MTL::Buffer*  get_trimming_params_buffer() = 0;
    virtual TrimmingParams * get_trimming_params(bool passU, uint32_t bucket_idx) = 0;
    virtual uint32_t get_trimming_params_offset(bool passU, uint32_t bucket_idx) = 0;

    virtual TrimmingParams * get_st2_trimming_params(uint32_t bucket_idx) = 0;
    virtual uint32_t get_st2_trimming_params_offset(uint32_t bucket_idx) const = 0;

    virtual void * get_mask() = 0;
    virtual uint32_t get_mask_size() = 0;
    virtual MTL::Buffer* get_mask_buffer() = 0;
    virtual uint32_t get_mask_offset() = 0;

    virtual void * get_collapse_edges() = 0;
    virtual uint32_t get_collapse_edges_size() = 0;
    virtual MTL::Buffer* get_collapse_edges_buffer() = 0;
    virtual uint32_t get_collapse_edges_offset() = 0;
#ifdef TRACK_COLLAPSED
    virtual void * get_collapse_edges_counts() = 0;
    virtual uint32_t get_collapse_edges_counts_size() = 0;
    virtual MTL::Buffer* get_collapse_edges_counts_buffer() = 0;
    virtual uint32_t get_collapse_edges_counts_offset() = 0;
#endif

    virtual void * get_bucket_positions() = 0;
    virtual uint32_t get_bucket_positions_size() = 0;
    virtual MTL::Buffer* get_bucket_positions_buffer() = 0;
    virtual uint32_t get_bucket_positions_offset() = 0;

    ////////////////////////////////////////////
    // Collapbe stack buffer expected to be same for all processing series, CPU/GPU accessible
    virtual MTL::Buffer* get_collapse_stash_buffer() = 0;
    virtual MTL::Buffer* take_collapse_stash_buffer() = 0;
    virtual uint32_t* get_collapse_stash() = 0;
};

uint32_t calc_aligned_offset(uint32_t offset);

void prepare_buffers(const MetalOps & metal, MTL::Buffer* primary_data_private,
            uint32_t mask_offset, uint32_t mask_size,
#ifdef TRACK_COLLAPSED
            uint32_t collapse_edges_counts_offset, uint32_t collapse_edges_counts_size,
#endif
            uint32_t bucket_positions_offset, uint32_t bucket_positions_size,
            uint32_t collapse_edges_offset, uint32_t collapse_edges_size);

#endif //METAL_PRIMARY_CONTEXT_H
