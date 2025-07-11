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

#include "metal_context_phase2.h"

#include "metal.h"

void MetalContextPhase2::init_phase2(const MetalOps & metal) {
    BUCKETS_NUM = metal.BUCKETS_NUM;
    PHASE1_EDGE_PER_BUCKET = metal.PHASE1_EDGE_PER_BUCKET;

    assert(phase1_results==nullptr);

    assert(primary_data_private==nullptr);

    assert(metal.EDGE_PER_BUCKET % 8 == 0);
    uint32_t half_mask_bytes = metal.EDGE_PER_BUCKET / 8 / 2; //
    assert(half_mask_bytes <= 4 * 1024*1024);

    mask_buffer_size = half_mask_bytes * 3;
    mask_buffer_offset = 0;

    collapse_edges_offset = calc_aligned_offset( mask_buffer_size + mask_buffer_offset );
    collapse_edges_size = metal.EDGE_PER_BUCKET*4;

    uint32_t buf_sz = calc_aligned_offset( collapse_edges_offset + collapse_edges_size);
#ifdef TRACK_COLLAPSED
    collapse_edges_counts_offset = calc_aligned_offset( collapse_edges_offset + collapse_edges_size);
    collapse_edges_counts_size = metal.EDGE_PER_BUCKET/2; // *4/8; - 1 byte per pair
    buf_sz = calc_aligned_offset( collapse_edges_counts_offset + collapse_edges_counts_size);
#endif

    primary_data_private = metal.device->newBuffer(buf_sz, RES_TYPE);
#ifndef NDEBUG
    memset(primary_data_private->contents(), 0x69, buf_sz);
#endif

    prepare_buffers(metal, primary_data_private,
            mask_buffer_offset, mask_buffer_size,
#ifdef TRACK_COLLAPSED
            collapse_edges_counts_offset, collapse_edges_counts_size,
#endif
            0, 0,
            collapse_edges_offset, collapse_edges_size);

    // Params data that will be initialized from CPU
    metrics_offset = 0;
#ifdef USE_METRICS
    metrics_size = sizeof(MetricsData);
#else
    metrics_size = 0;
#endif

    trim_params_offset = calc_aligned_offset( metrics_offset + metrics_size);
    trim_params_size = calc_aligned_offset(sizeof(TrimmingParams)); // size per one item. Total BUCKETS_NUM + BUCKETS_NUM*PRIMARY_TRIM_STEPS

    uint32_t buf2_sz = trim_params_offset + trim_params_size * BUCKETS_NUM; // u/v for all combinations
    primary_data_shared = metal.device->newBuffer(buf2_sz, MTL::ResourceStorageModeShared);

#ifndef NDEBUG
    memset(primary_data_shared->contents(), 0x69, buf2_sz);
#endif

#ifdef USE_METRICS
    memset( static_cast<char*>(primary_data_shared->contents()) + metrics_offset, 0, sizeof(MetricsData) );
#endif
}

void MetalContextPhase2::reset_phase2(FindCycleStartResult & phase1_res) {
    if (phase1_results) {
        phase1_results->release();
    }
    phase1_results = phase1_res.edges;
    phase1_res.edges = nullptr;

    if (collapse_stash) {
        collapse_stash->release();
    }
    assert(phase1_res.collapse_stash_buffer);
    collapse_stash = phase1_res.collapse_stash_buffer;
    phase1_res.collapse_stash_buffer = nullptr;

    assert(collapse_stash);
}

void MetalContextPhase2::release() {
    RELEASE(phase1_results);
    RELEASE(primary_data_private);
    RELEASE(primary_data_shared);
    RELEASE(collapse_stash);
}

MemPool * MetalContextPhase2::get_mem_pool() {
    assert(false);
    return nullptr;
}

MTL::Buffer* MetalContextPhase2::get_buffer0() {
    assert(phase1_results);
    return phase1_results;
}

MTL::Buffer* MetalContextPhase2::get_buffer1() {
    return nullptr;
}

Step1Params * MetalContextPhase2::get_st1_nonce28_params() {
    assert(false);
    return nullptr;
}
MTL::Buffer* MetalContextPhase2::get_st1_nonce28_params_buffer() {
    assert(false);
    return nullptr;
}
uint32_t MetalContextPhase2::get_st1_nonce28_params_offset() {
    assert(false);
    return 0;
}

Step2Params * MetalContextPhase2::get_st2_params(uint32_t bucket_idx) {
    assert(false);
    return nullptr;
}
MTL::Buffer* MetalContextPhase2::get_st2_params_buffer() {
    assert(false);
    return nullptr;
}
uint32_t MetalContextPhase2::get_st2_params_offset(uint32_t bucket_idx) const {
    assert(false);
    return 0;
}
#ifdef USE_METRICS
MetricsData * MetalContextPhase2::get_metrics() {
    return reinterpret_cast<MetricsData*>( static_cast<char*>(primary_data_shared->contents()) + get_metrics_offset() );
}
MTL::Buffer*  MetalContextPhase2::get_metrics_buffer() {return primary_data_shared;}
uint32_t  MetalContextPhase2::get_metrics_offset() {return metrics_offset;}
#endif

MTL::Buffer*  MetalContextPhase2::get_trimming_params_buffer() {return primary_data_shared;}
TrimmingParams * MetalContextPhase2::get_trimming_params(bool passU, uint32_t bucket_idx) {
    return reinterpret_cast<TrimmingParams*>( static_cast<char*>(primary_data_shared->contents()) + get_trimming_params_offset(passU, bucket_idx) );
}
uint32_t MetalContextPhase2::get_trimming_params_offset(bool passU, uint32_t bucket_idx) {
    assert(bucket_idx<BUCKETS_NUM);
    return trim_params_offset + trim_params_size * bucket_idx;
}

TrimmingParams * MetalContextPhase2::get_st2_trimming_params(uint32_t bucket_idx) {
    assert(false);
    return nullptr;
}
uint32_t MetalContextPhase2::get_st2_trimming_params_offset(uint32_t bucket_idx) const {
    assert(false);
    return 0;
}

void * MetalContextPhase2::get_mask() {
    return static_cast<char*>(primary_data_private->contents()) + mask_buffer_offset;
}
uint32_t MetalContextPhase2::get_mask_size() {return mask_buffer_size;}
MTL::Buffer* MetalContextPhase2::get_mask_buffer() {return primary_data_private;}
uint32_t MetalContextPhase2::get_mask_offset() {return mask_buffer_offset;}

void * MetalContextPhase2::get_collapse_edges() {
    return static_cast<char*>(primary_data_private->contents()) + collapse_edges_offset;
}
uint32_t MetalContextPhase2::get_collapse_edges_size() {return collapse_edges_size;}
MTL::Buffer* MetalContextPhase2::get_collapse_edges_buffer() {return primary_data_private;}
uint32_t MetalContextPhase2::get_collapse_edges_offset() {return collapse_edges_offset;}
#ifdef TRACK_COLLAPSED
void * MetalContextPhase2::get_collapse_edges_counts() {
    return static_cast<char*>(primary_data_private->contents()) + collapse_edges_counts_offset;
}
uint32_t MetalContextPhase2::get_collapse_edges_counts_size() {return collapse_edges_counts_size;}
MTL::Buffer* MetalContextPhase2::get_collapse_edges_counts_buffer() {return primary_data_private;}
uint32_t MetalContextPhase2::get_collapse_edges_counts_offset() {return collapse_edges_counts_offset;}
#endif

void * MetalContextPhase2::get_bucket_positions()  {
    assert(phase1_results);
    return static_cast<char*>(phase1_results->contents()) + get_bucket_positions_offset();
}
uint32_t MetalContextPhase2::get_bucket_positions_size() {return BUCKETS_NUM*BUCKETS_NUM*4;}
MTL::Buffer* MetalContextPhase2::get_bucket_positions_buffer() {return phase1_results;}
uint32_t MetalContextPhase2::get_bucket_positions_offset() {return BUCKETS_NUM*BUCKETS_NUM*PHASE1_EDGE_PER_BUCKET*8;}


MTL::Buffer* MetalContextPhase2::get_collapse_stash_buffer() {
    assert(collapse_stash);
    return collapse_stash;
}

MTL::Buffer* MetalContextPhase2::take_collapse_stash_buffer() {
    // not expected to be called
    assert(false);
    return nullptr;
}

uint32_t* MetalContextPhase2::get_collapse_stash() {
    assert(collapse_stash);
    return (uint32_t *) collapse_stash->contents();
}