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

#include "metal_context_phase1.h"
#include "metal.h"
#include <iostream>
#include "features.h"

void MetalContextPhase1::init_phase1(const MetalOps & metal) {
    BUCKETS_NUM = metal.BUCKETS_NUM;
    MEM_SIZE = metal.MEM_SIZE;

    event = metal.device->newSharedEvent();

    assert(metal.MEM_SIZE%2==0);
    buffer0 = metal.device->newBuffer(MEM_SIZE/2, RES_TYPE);
    buffer1 = metal.device->newBuffer(MEM_SIZE/2, RES_TYPE);
#ifndef NDEBUG
    memset( buffer0->contents(), 0x63, buffer0->length() );
    memset( buffer1->contents(), 0x63, buffer1->length() );
#endif

    mem_pool.reset(2, MEM_SIZE);

#ifndef NDEBUG
#ifdef SHOW_TRACING
    std::cout << "Starting, memory usage: " << mem_pool.getMemoryStatus() << std::endl;
#endif
#endif

    assert(metal.EDGE_PER_BUCKET % 8 == 0);
    uint32_t half_mask_bytes = metal.EDGE_PER_BUCKET / 8 / 2; //
    assert(half_mask_bytes <= 4 * 1024*1024);

    // Partition and allocating primary_data_private
    mask_buffer_offset = 0;
    mask_buffer_size = half_mask_bytes * 3;

    bucket_positions_offset = calc_aligned_offset( mask_buffer_offset + mask_buffer_size );
    bucket_positions_size = BUCKETS_NUM*BUCKETS_NUM * 4;

    collapse_edges_offset = calc_aligned_offset( bucket_positions_offset + bucket_positions_size );
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
            bucket_positions_offset, bucket_positions_size,
            collapse_edges_offset, collapse_edges_size);

    // Params data that will be initialized from CPU
    metrics_offset = 0;
#ifdef USE_METRICS
    metrics_size = sizeof(MetricsData);
#else
    metrics_size = 0;
#endif

    st1_nonce28_params_offset = calc_aligned_offset( metrics_offset + metrics_size );
    st1_nonce28_params_size = calc_aligned_offset(sizeof(Step1Params));

    st2_params_offset = st1_nonce28_params_offset + st1_nonce28_params_size; //
    st2_params_size = calc_aligned_offset(sizeof(Step2Params)); // size for a one item, total BUCKETS_NUM

    trim_params_offset = st2_params_offset + st2_params_size * BUCKETS_NUM;
    trim_params_size = calc_aligned_offset(sizeof(TrimmingParams)); // size per one item. Total BUCKETS_NUM + BUCKETS_NUM*PRIMARY_TRIM_STEPS

    uint32_t buf2_sz = trim_params_offset + trim_params_size * BUCKETS_NUM * (2+1); // u/v for all combinations
    primary_data_shared = metal.device->newBuffer(buf2_sz, MTL::ResourceStorageModeShared);

#ifndef NDEBUG
    memset(primary_data_shared->contents(), 0x69, buf2_sz);
#endif

#ifdef USE_METRICS
    memset( static_cast<char*>(primary_data_shared->contents()) + metrics_offset, 0, sizeof(MetricsData) );
#endif
}

void MetalContextPhase1::reset_phase1(const MetalOps & metal) {
    if (collapse_stash==nullptr) {
        collapse_stash = metal.device->newBuffer( (metal.COLLAPSE_STASH_SIZE*2 + 1)*4, MTL::ResourceStorageModeShared);
    }
#ifndef NDEBUG
    memset(collapse_stash->contents(), 0x69, collapse_stash->length());
#endif
    ((uint32_t *) collapse_stash->contents())[0] = 0;

    mem_pool.reset(2, MEM_SIZE);
}

void MetalContextPhase1::release() {
    RELEASE(collapse_stash);
    RELEASE(buffer0);
    RELEASE(buffer1);
    RELEASE(primary_data_private);
    RELEASE(primary_data_shared);
    RELEASE(event);
}

MemPool * MetalContextPhase1::get_mem_pool()  {return &mem_pool;}
MTL::Buffer* MetalContextPhase1::get_buffer0()  {return buffer0;}
MTL::Buffer* MetalContextPhase1::get_buffer1()  {return buffer1;}

Step1Params * MetalContextPhase1::get_st1_nonce28_params() {
    return reinterpret_cast<Step1Params*>( static_cast<char*>(primary_data_shared->contents()) + st1_nonce28_params_offset );
}
MTL::Buffer* MetalContextPhase1::get_st1_nonce28_params_buffer() {return primary_data_shared;}
uint32_t MetalContextPhase1::get_st1_nonce28_params_offset() {return st1_nonce28_params_offset;}

Step2Params * MetalContextPhase1::get_st2_params(uint32_t bucket_idx) {
    return reinterpret_cast<Step2Params*>( static_cast<char*>(primary_data_shared->contents()) + get_st2_params_offset(bucket_idx) );
}
MTL::Buffer* MetalContextPhase1::get_st2_params_buffer() {return primary_data_shared;}
uint32_t MetalContextPhase1::get_st2_params_offset(uint32_t bucket_idx) const {
    return st2_params_offset + st2_params_size * bucket_idx;
}

#ifdef USE_METRICS
MetricsData * MetalContextPhase1::get_metrics() {
    return reinterpret_cast<MetricsData*>( static_cast<char*>(primary_data_shared->contents()) + get_metrics_offset() );
}
MTL::Buffer*  MetalContextPhase1::get_metrics_buffer() {return primary_data_shared;}
uint32_t  MetalContextPhase1::get_metrics_offset() {return metrics_offset;}
#endif

MTL::Buffer*  MetalContextPhase1::get_trimming_params_buffer() {return primary_data_shared;}
TrimmingParams * MetalContextPhase1::get_trimming_params(bool passU, uint32_t bucket_idx) {
    return reinterpret_cast<TrimmingParams*>( static_cast<char*>(primary_data_shared->contents()) + get_trimming_params_offset(passU, bucket_idx) );
}
uint32_t MetalContextPhase1::get_trimming_params_offset(bool passU, uint32_t bucket_idx) {
    return trim_params_offset + trim_params_size * ((1 + (passU ? 0 : 1)) * BUCKETS_NUM + bucket_idx);
}

TrimmingParams * MetalContextPhase1::get_st2_trimming_params(uint32_t bucket_idx) {
    return reinterpret_cast<TrimmingParams*>( static_cast<char*>(primary_data_shared->contents()) + get_st2_trimming_params_offset(bucket_idx) );
}
uint32_t MetalContextPhase1::get_st2_trimming_params_offset(uint32_t bucket_idx) const  {
    return trim_params_offset + trim_params_size * bucket_idx;
}

void * MetalContextPhase1::get_mask() {
    return static_cast<char*>(primary_data_private->contents()) + mask_buffer_offset;
}
uint32_t MetalContextPhase1::get_mask_size() {return mask_buffer_size;}
MTL::Buffer* MetalContextPhase1::get_mask_buffer() {return primary_data_private;}
uint32_t MetalContextPhase1::get_mask_offset() {return mask_buffer_offset;}

void * MetalContextPhase1::get_collapse_edges() {
    return static_cast<char*>(primary_data_private->contents()) + collapse_edges_offset;
}
uint32_t MetalContextPhase1::get_collapse_edges_size() {return collapse_edges_size;}
MTL::Buffer* MetalContextPhase1::get_collapse_edges_buffer() {return primary_data_private;}
uint32_t MetalContextPhase1::get_collapse_edges_offset() {return collapse_edges_offset;}

#ifdef TRACK_COLLAPSED
void * MetalContextPhase1::get_collapse_edges_counts() {
    return static_cast<char*>(primary_data_private->contents()) + collapse_edges_counts_offset;
}
uint32_t MetalContextPhase1::get_collapse_edges_counts_size() {return collapse_edges_counts_size;}
MTL::Buffer* MetalContextPhase1::get_collapse_edges_counts_buffer() {return primary_data_private;}
uint32_t MetalContextPhase1::get_collapse_edges_counts_offset() {return collapse_edges_counts_offset;}
#endif

void * MetalContextPhase1::get_bucket_positions()  {
    return static_cast<char*>(primary_data_private->contents()) + bucket_positions_offset;
}
uint32_t MetalContextPhase1::get_bucket_positions_size() {return bucket_positions_size;}
MTL::Buffer* MetalContextPhase1::get_bucket_positions_buffer() {return primary_data_private;}
uint32_t MetalContextPhase1::get_bucket_positions_offset() {return bucket_positions_offset;}

MTL::Buffer* MetalContextPhase1::get_collapse_stash_buffer() {
    assert(collapse_stash);
    return collapse_stash;
}

MTL::Buffer* MetalContextPhase1::take_collapse_stash_buffer() {
    assert(collapse_stash);
    MTL::Buffer* buf = collapse_stash;
    collapse_stash = nullptr;
    return buf;
}

uint32_t* MetalContextPhase1::get_collapse_stash() {
    assert(collapse_stash);
    return (uint32_t *) collapse_stash->contents();
}

MTL::SharedEvent* MetalContextPhase1::get_event() {
    return event;
}
int MetalContextPhase1::generate_next_event() {
    return ++event_counter;
}
