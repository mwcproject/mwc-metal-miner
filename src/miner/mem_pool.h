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

#ifndef MEM_POOL_H
#define MEM_POOL_H

#include <Metal/MTLBuffer.hpp>
#include <Metal/MTLCommandBuffer.hpp>
#include <cassert>
#include <vector>
#include "metal_context.h"

// Memory manage operate by pages. Allocated blocks are expected to be large
#define MEM_POOL_UNITS 4096

struct MetalEncoderManager;

// Memory pool based on buffers. We need to manage memory somehow.
struct MemRange {
    uint32_t buffer; // index of the buffer (0 or 1)
    uint32_t start;  // Starting page
    uint32_t end;    // endign page of the block
    uint64_t length_bytes; // Expected number of bytes, expected to fit into the allocated pages

    MemRange() = default;
    MemRange(uint32_t _buffer, uint32_t _start, uint32_t _end, uint64_t _length_bytes) : buffer(_buffer), start(_start), end(_end), length_bytes(_length_bytes) {
        assert( (length_bytes + MEM_POOL_UNITS-1)/MEM_POOL_UNITS == end-start );
    }
    MemRange(const MemRange &other) = default;
    MemRange &operator=(const MemRange &other) = default;

    inline uint32_t length() const {assert(start<end); return end - start;}
    inline uint64_t get_length_bytes() const {
        assert(length_bytes<=uint64_t(length()) * MEM_POOL_UNITS);
        assert(length_bytes>=uint64_t(length()) * MEM_POOL_UNITS - MEM_POOL_UNITS);
        return length_bytes;
    }

    // Convert to address that is used on Metal side. It is
    // expected 2 buffers, the high bit is set for that.
    inline uint32_t to_allocated_block() const {
        uint32_t res = (buffer==0 ? 0 : 0x80000000) | (start);
        return res;
    }

    // String representaiton for debugging
    std::string to_string() const;

    inline void * get_data_ptr(MetalContext & context) const  {
        assert(buffer < 2);
        MTL::Buffer * data_buf = get_buffer(context);
        uint8_t * ptr = (uint8_t *) data_buf->contents();
        assert(start<end);
        return ptr + uint64_t(start)*MEM_POOL_UNITS;
    }

    // Get buffer that represent this block
    inline MTL::Buffer * get_buffer( MetalContext & context) const {
        return  buffer == 0 ? context.get_buffer0() : context.get_buffer1();
    }

    // Copy data from Metal buffer into trg_data. Needed mostly for stale cache problem.
    void copy_data_into( MetalContext & context, std::vector<uint64_t> & trg_data ) const;
};

// Momory pool for our metal tasks.
// Pool units are pages, expetec that large bloacks are requested and released
class MemPool {
    // Available memory blocks
    std::vector<MemRange> ranges;
    // Deallocated blocks, can be reused if needed
    std::vector<std::pair<MemRange, MTL::CommandBuffer *> > deallocated;

    // Available memory counter. Mostly for debug/monitoring
    uint64_t total_mem = 0;

    // Access mutex
    mutable std::mutex mtx;
public:
    MemPool() {}
    ~MemPool();

    // Reset memory pool - all buffers data is available
    void reset(int buffers, uint64_t total_pool_size_bytes);

    // Allocate a new memroy block. Will crash in case of not enough memory.
    MemRange allocate(uint64_t length_byte, MetalEncoderManager & enc_man);

    // Release block
    void release(const MemRange & range, MetalContext & context, MTL::CommandBuffer *command);
    void release(const std::vector<MemRange> & range, MetalContext & context, MTL::CommandBuffer *command);

    // Dump of available blocks. For debug only
    std::string getAvailableBlockDump() const;

    // Status of available blocks. For debug only
    std::string getMemoryStatus() const;
private:
    MemRange allocate_impl(uint64_t length_byte, bool wait_for_completion, MetalEncoderManager & enc_man);
};



#endif //MEM_POOL_H
