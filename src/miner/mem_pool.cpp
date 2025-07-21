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

#include "mem_pool.h"

#include <iostream>
#include <ostream>
#include "metal.h"

std::string MemRange::to_string() const {
    return  "Range(" + std::to_string(buffer) + "," + std::to_string(start) + ", " + std::to_string(end) + ")";
}

void MemRange::copy_data_into( MetalContext & context, std::vector<uint64_t> & trg_data ) const {
    void * data = get_data_ptr(context);
    uint64_t len = get_length_bytes();
    assert(len % 8 == 0);
    trg_data.resize(len/8);
    memcpy(trg_data.data(), data, len);
}


///////////////////////////////////////////////////////////////////////////////////////////

MemPool::~MemPool() {
}

void MemPool::reset(int buffers, uint64_t total_pool_size_bytes) {
    std::lock_guard<std::mutex> lock(mtx);

    ranges.clear();
    deallocated.clear();
    total_mem = 0;

    assert((total_pool_size_bytes / buffers) % MEM_POOL_UNITS == 0);

    uint32_t end = total_pool_size_bytes / buffers / MEM_POOL_UNITS;

    for (uint32_t i = 0; i < buffers; i++) {
        ranges.push_back(MemRange(i, 0, end, uint64_t(end)*MEM_POOL_UNITS));
        total_mem += uint64_t(end)*MEM_POOL_UNITS;
    }
}


MemRange MemPool::allocate(uint64_t length_byte, MetalEncoderManager & enc_man) {
    std::lock_guard<std::mutex> lock(mtx);
    return allocate_impl(length_byte, false, enc_man);
}

MemRange MemPool::allocate_impl(uint64_t length_byte, bool wait_for_completion, MetalEncoderManager & enc_man) {
    int mem_len = (length_byte + MEM_POOL_UNITS - 1) / MEM_POOL_UNITS;

    // allocating from the first available unit
    for (uint32_t i = 0; i < ranges.size(); i++) {
        MemRange & r = ranges[i];
        if (r.end - r.start >= mem_len) {
            MemRange res;
            res.buffer = r.buffer;
            res.start = r.start;
            res.end = res.start + mem_len;
            res.length_bytes = length_byte;
            r.start += mem_len;
            r.length_bytes -= mem_len * MEM_POOL_UNITS;
            assert(res.get_length_bytes()>0);

            if (r.start == r.end) {
                ranges.erase(ranges.begin() + i);
            }
            else {
                assert(r.get_length_bytes()>0);
                for (uint32_t j = i; j > 0; j--) {
                    if (ranges[j-1].length() > ranges[j].length()) {
                        std::swap(ranges[j], ranges[j-1]);
                    }
                    else {
                        break;
                    }
                }
            }

            return res;
        }
    }

    // nothing was allocated
    if (deallocated.empty()) {
        std::cout << "ERROR!!! Not able to allocate " << length_byte << " of memory." << std::endl;
        std::cout << "Available blocks: " << getAvailableBlockDump() << std::endl;
        exit(1);
    }

    std::vector<MemRange> new_ranges;

    if (wait_for_completion) {
        MTL::CommandBuffer * cmd_buf = enc_man.take_command_buffer();
        if (cmd_buf) {
            // no reason to keep it here
            enc_man.stash_buf(cmd_buf);
        }
    }

    for (std::pair<MemRange, MTL::CommandBuffer *> & deal : deallocated) {
#ifndef METAL_DO_WAITS
        MTL::CommandBuffer * buf = std::get<1>(deal);
        if (buf==nullptr || buf->status() == MTL::CommandBufferStatusCompleted) {
            std::get<1>(deal) = nullptr;
            new_ranges.push_back( std::get<0>(deal));
        }
        else {
            if (wait_for_completion) {
                buf->waitUntilCompleted();
                std::get<1>(deal) = nullptr;
                new_ranges.push_back( std::get<0>(deal));
                break;
            }
            else {
                continue;
            }
        }
#else
        new_ranges.push_back( std::get<0>(deal));
#endif
    }
    deallocated.erase( std::remove_if(deallocated.begin(), deallocated.end(), [](auto deal) {return std::get<1>(deal)==nullptr;} ), deallocated.end() );

    new_ranges.insert(new_ranges.end(), ranges.begin(), ranges.end() );
    std::sort( new_ranges.begin(), new_ranges.end(), [](const MemRange & a, const MemRange & b) {
        if (a.buffer == b.buffer)
            return a.start < b.start;
        else
            return a.buffer < b.buffer;
    } );
    ranges.resize(0);
    ranges.push_back(new_ranges[0]);
    assert(ranges.back().get_length_bytes()>0);
    for (int i = 1; i < new_ranges.size(); i++) {
        if (ranges.back().buffer == new_ranges[i].buffer && ranges.back().end == new_ranges[i].start) {
            ranges.back().end = new_ranges[i].end;
            ranges.back().length_bytes = ranges.back().length() * MEM_POOL_UNITS;
        }
        else {
            ranges.push_back(new_ranges[i]);
        }
        assert(ranges.back().get_length_bytes()>0);
    }

    std::sort( ranges.begin(), ranges.end(), [](const MemRange & a, const MemRange & b) {
        return a.length() < b.length();
    } );

    return allocate_impl(length_byte, true, enc_man);
}

void MemPool::release(const MemRange & range, MetalContext & context, MTL::CommandBuffer *command) {
    std::lock_guard<std::mutex> lock(mtx);

#ifdef METAL_DO_WAITS
#ifndef NDEBUG
    void * dt = range.get_data_ptr(context);
    memset(dt, 0x63, range.get_length_bytes());
#endif
#endif

    deallocated.push_back( std::pair<MemRange, MTL::CommandBuffer *>(range, command));
}

void MemPool::release(const std::vector<MemRange> & range, MetalContext & context, MTL::CommandBuffer *command) {
    std::lock_guard<std::mutex> lock(mtx);

#ifdef METAL_DO_WAITS
#ifndef NDEBUG
    for (const auto r : range) {
        void * dt = r.get_data_ptr(context);
        memset(dt, 0x63, r.get_length_bytes());
    }
#endif
#endif

    for (const auto & r : range) {
        deallocated.push_back( std::pair<MemRange, MTL::CommandBuffer *>(r, command));
    }
}


std::string MemPool::getAvailableBlockDump() const {
    std::lock_guard<std::mutex> lock(mtx);

    std::string res = "";
    for (const MemRange & r : ranges) {
        if (res.size() > 0)
            res += ", ";
        res += r.to_string();
    }
    return res;
}

std::string MemPool::getMemoryStatus() const {
    std::lock_guard<std::mutex> lock(mtx);

    uint64_t avail_mem = 0;
    for (const MemRange & r : ranges) {
        avail_mem += r.get_length_bytes();
    }

    for (const std::pair<MemRange, MTL::CommandBuffer *> & r : deallocated) {
        avail_mem += std::get<0>(r).get_length_bytes();
    }

    return std::string("Available memory (Gb): ") + std::to_string(avail_mem / 1024.0 / 1024.0 / 1024.0) +
        " from " + std::to_string(total_mem / 1024.0 / 1024.0 / 1024.0);
}

