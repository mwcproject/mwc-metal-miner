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

#include <Metal/Metal.hpp>
#include <iostream>
#include <chrono>
#include <vector>

int main() {
    // Create a Metal device
    MTL::Device* device = MTL::CreateSystemDefaultDevice();
    if (!device) {
        std::cerr << "Metal is not supported on this device." << std::endl;
        return 1;
    }

    // Query device properties
    size_t maxBufferLength = device->maxBufferLength();
    std::cout << "Max buffer length: " << maxBufferLength << " bytes" << std::endl;

    size_t recommendedMaxWorkingSetSize = device->recommendedMaxWorkingSetSize();
    std::cout << "Recommended max working set size: " << recommendedMaxWorkingSetSize << " bytes" << std::endl;

    // Parameters for allocation
    const size_t NUM_BUFFERS = 512;
    const size_t BUFFER_SIZE = 16384*512; // 16 KB per buffer
    std::vector<MTL::Buffer*> buffers(NUM_BUFFERS);

    // Start timing
    auto start = std::chrono::high_resolution_clock::now();

    // Allocate buffers
    for (size_t i = 0; i < NUM_BUFFERS; ++i) {
        buffers[i] = device->newBuffer(BUFFER_SIZE, MTL::ResourceStorageModeShared);
        if (!buffers[i]) {
            std::cerr << "Failed to allocate buffer " << i << std::endl;
            return 1;
        }
    }

    // Stop timing
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate total time
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Time to allocate " << NUM_BUFFERS << " buffers: " << duration << " ms" << std::endl;

    // Calculate total memory used
    size_t totalMemoryUsed = NUM_BUFFERS * BUFFER_SIZE;
    std::cout << "Total memory used: " << totalMemoryUsed << " bytes (" << (totalMemoryUsed / (1024 * 1024)) << " MB)" << std::endl;

    // Release buffers
    for (size_t i = 0; i < NUM_BUFFERS; ++i) {
        buffers[i]->release();
    }

    // Release the device
    device->release();

    auto end2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - end).count();
    std::cout << "Time to free " << NUM_BUFFERS << " buffers: " << duration2 << " ms" << std::endl;

    return 0;
}