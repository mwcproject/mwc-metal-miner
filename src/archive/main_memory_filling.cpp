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

#include <iostream>
#include <thread>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <sys/sysctl.h>
#include <ctime>
#include <assert.h>

uint64_t get_timestamp_ms() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);

    // Convert to nanoseconds
    return static_cast<uint64_t>(ts.tv_sec) * 1000 + ts.tv_nsec/1000000;
}

// Function to get the number of available CPU cores
int get_num_cores() {
    int num_cores;
    size_t len = sizeof(num_cores);
    sysctlbyname("hw.ncpu", &num_cores, &len, nullptr, 0);
    return num_cores;
}

uint64_t total_sum = 0;

inline void prefetch_to_l1(const void *addr) {
    asm volatile(
        "prfm pldl1keep, [%[addr]]\n"  // Prefetch data into L1 cache
        :
        : [addr] "r"(addr)             // Input operand (memory address)
        :                              // No clobbered registers
    );
}

inline void evict_from_cache(const void *addr) {
    asm volatile(
        "dc cvau, %[addr]\n"  // Clean and invalidate cache line by VA
        :
        : [addr] "r"(addr)    // Input operand (memory address)
        :                      // No clobbered registers
    );
}
// Function to fill a memory segment with zeroes using ARMv9 assembly
void fill_zeroes(void* start, size_t size) {

    /*uint64_t* mem_start = reinterpret_cast<uint64_t*>(start);
    uint64_t* mem_end = mem_start + size/sizeof(uint64_t);

    uint64_t sum = 0;
    int cnt = 0;
    while (mem_start < mem_end) {
        //prefetch_to_l1(mem_start);
        //if (cnt>1)
          //  evict_from_cache(mem_start-64/8);
        //for (int i=0; i<64/8; i++) {
            sum+= *mem_start++;
        //}
        //cnt++;
    }
    total_sum+=sum;*/

    memset(start, 0, size);
/*  asm volatile (
        "mov x0, %0\n"          // x0 = start address
        "mov x1, %1\n"          // x1 = size in bytes
        "mov w2, #0\n"
        "movi v0.16b, #0\n"
        "movi v1.16b, #0\n"
        "1:\n"
        //"stp xzr, xzr, [x0], #16\n" // Store 16 bytes of zeroes and increment address
//        "STP    q0, q1, [x0], #32\n" //
        //"subs x1, x1, #16\n"    // Subtract 16 from size
        "str w2, [x0], #4\n"
        "ADD w2, w2, #1\n"      // Increment counter
        "subs x1, x1, #4\n"
        "b.gt 1b\n"             // If size > 0, repeat
        :
        : "r"(start), "r"(size)
        : "x0", "x1", "x2", "v0", "v1", "memory"
    );*/

    // prefetch and keep in memory
/*    asm volatile (
        "mov x0, %0\n"          // x0 = start address
        "mov x1, %1\n"          // x1 = size in bytes
        "movi v0.16b, #0\n"     // Zero register q0
        "movi v1.16b, #0\n"     // Zero register q1

        // Prefetch the first 32 cache lines (assuming 64 bytes per line)
        "mov x2, x0\n"          // x2 = prefetch address
        "mov x3, #0\n"          // x3 = loop counter
    "prefetch_loop:\n"
        "PRFM PSTL1KEEP, [x2]\n" // for store
        "ADD x2, x2, #64\n"     // Move to the next cache line (64 bytes)
        "ADD x3, x3, #1\n"
        "CMP x3, #4\n"         // Check if 32 cache lines are prefetched
        "B.LT prefetch_loop\n"

    "zeroing_loop:\n"
        "STP q0, q1, [x0, #32]\n"  // Store 32 bytes of zeroes and increment address
        "STP q0, q1, [x0], #64\n"  // Store 32 bytes of zeroes and increment address

        "SUB x4, x0, #32\n"         // x4 = x0 - 32
        "PRFM PSTL1STRM, [x4]\n"    // Evict the processed cache line

        // Prefetch the next 32nd cache line
        "PRFM PSTL1KEEP, [x0, #256]\n"
        //"PRFM PSTL2KEEP, [x0, #1024]\n"

        "subs x1, x1, #64\n"    // Subtract 32 from size
        "b.gt zeroing_loop\n"             // If size > 0, repeat
        :
        : "r"(start), "r"(size)
        : "x0", "x1", "x2", "x3", "x4", "v0", "v1", "memory"
    );*/
}

int main() {
    const size_t memory_size = 4ULL * 1024 * 1024 * 1024; // 4GB

    // Allocate 4GB of memory with 64-byte alignment
    void* memory = aligned_alloc(64, memory_size+4096);
    if (!memory) {
        std::cerr << "Failed to allocate memory!" << std::endl;
        return 1;
    }
    //memory = std::assume_aligned<64>(memory);

    // Get the number of available cores
    int num_cores = 8;//get_num_cores();
    std::cout << "Number of cores: " << num_cores << std::endl;

    // Calculate the segment size for each thread
    size_t segment_size = memory_size / num_cores;

    // Create threads to fill memory with zeroes
    uint64_t startTime = get_timestamp_ms();
    std::vector<std::thread> threads;
    for (int i = 0; i < num_cores; ++i) {
        void* segment_start = static_cast<char*>(memory) + i * segment_size;
        threads.emplace_back(fill_zeroes, segment_start, segment_size);
    }

    // Wait for all threads to finish
    for (auto& thread : threads) {
        thread.join();
    }

    uint64_t finishTime = get_timestamp_ms();
    std::cout << "Memory filled with zeroes successfully!  Timing: " << (finishTime-startTime) << "ms" << std::endl;

    if (total_sum == 0) {
      return 0;
    }

    assert( static_cast<char*>(memory)[12345] == 0 );

    // Free the allocated memory
    free(memory);

    return 0;
}