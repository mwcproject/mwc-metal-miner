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


#include "metal_context.h"
#include "metal.h"
#include <iostream>
#include "features.h"


uint32_t calc_aligned_offset(uint32_t offset) {
    return ((offset + 63)/64)*64;
}

void prepare_buffers(const MetalOps & metal, MTL::Buffer* primary_data_private,
            uint32_t mask_offset, uint32_t mask_size,
#ifdef TRACK_COLLAPSED
            uint32_t collapse_edges_counts_offset, uint32_t collapse_edges_counts_size,
#endif
            uint32_t bucket_positions_offset, uint32_t bucket_positions_size,
            uint32_t collapse_edges_offset, uint32_t collapse_edges_size) {
    // let's init the data
    MTL::CommandBuffer* command = metal.primary_queue->commandBuffer();

    { // Clean the masks
        MTL::ComputeCommandEncoder* encoder = command->computeCommandEncoder();
        encoder->setComputePipelineState(metal.f_clean_data.pipeline);
        encoder->setBuffer(primary_data_private, mask_offset, 0);

        // Dispatch threads for each buffer (assuming 1 result per buffer)
        assert( mask_size % (4*64) == 0);
        encoder->dispatchThreads(MTL::Size( mask_size/(4*64), 1, 1), MTL::Size(1, 1, 1));
        encoder->endEncoding();
        encoder->retain();
        encoder->release();
    }

#ifdef TRACK_COLLAPSED
    {
        MTL::ComputeCommandEncoder* encoder = command->computeCommandEncoder();
        encoder->setComputePipelineState(metal.f_clean_data.pipeline);
        encoder->setBuffer(primary_data_private, collapse_edges_counts_offset, 0);

        // Dispatch threads for each buffer (assuming 1 result per buffer)
        assert( collapse_edges_counts_size % (4*64) == 0);
        encoder->dispatchThreads(MTL::Size( collapse_edges_counts_size/(4*64), 1, 1), MTL::Size(1, 1, 1));
        encoder->endEncoding();
        encoder->retain();
        encoder->release();
    }
#endif

#ifndef NDEBUG
    if (bucket_positions_size>0)
    {
        MTL::ComputeCommandEncoder* encoder = command->computeCommandEncoder();
        encoder->setComputePipelineState(metal.f_clean_data.pipeline);
        encoder->setBuffer(primary_data_private, bucket_positions_offset, 0);

        // Dispatch threads for each buffer (assuming 1 result per buffer)
        encoder->dispatchThreads(MTL::Size( bucket_positions_size/(64*4) + 1 , 1, 1), MTL::Size(1, 1, 1));
        encoder->endEncoding();
        encoder->retain();
        encoder->release();
    }
    {
        // add init for collapse_edges here. It is a one time job
        MTL::ComputeCommandEncoder* encoder = command->computeCommandEncoder();
        encoder->setComputePipelineState(metal.f_clean_data_FF.pipeline);
        encoder->setBuffer(primary_data_private, collapse_edges_offset, 0);

        // Dispatch threads for each buffer (assuming 1 result per buffer)
        assert( collapse_edges_size % (4*64) == 0);
        encoder->dispatchThreads(MTL::Size( collapse_edges_size/4/64, 1, 1), MTL::Size(1, 1, 1));
        encoder->endEncoding();
        encoder->retain();
        encoder->release();
    }
#endif
    // one time task, can wait
    command->commit();
    command->waitUntilCompleted();
    command->retain();
    command->release();
}