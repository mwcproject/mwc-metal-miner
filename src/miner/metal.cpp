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

#include "metal.h"
#include <iostream>

#define ELEMENTS_PER_TASK  64


// NONCE_0 starting value for the nonce, NONCE_K is 2, NONCE_N 0 or 1 for U/V rows
static std::string generate_sip_hash_src( uint32_t EDGE_BITS, uint32_t ELEMENT_SIZE,
		uint32_t BUCKETS_NUM, uint32_t BUCKET_SIZE )
{
	return std::string(R"(
#include <metal_stdlib>
using namespace metal;

	)") +

"#define EDGE_BITS " + std::to_string(EDGE_BITS) + "\n" +
"#define ELEMENT_SIZE " + std::to_string(ELEMENT_SIZE) + "\n" +
"#define BUCKET_SIZE " + std::to_string(BUCKET_SIZE) + "\n" +
"#define BUCKETS_NUM " + std::to_string(BUCKETS_NUM) + "\n" +

	R"(
#define NUMBER_OF_EDGES (static_cast<ulong>(1) << EDGE_BITS)
#define NODE_MASK (NUMBER_OF_EDGES - 1)

// SipRound
inline void sipRound(thread ulong4 &keys) {

	// Perform SipRound on keys
	keys.even += keys.odd;
	keys.odd = rotate(keys.odd, ulong2(13, 16));
	keys.odd ^= keys.even;
	keys.x = rotate(keys.x, static_cast<ulong>(32));
	keys.even += keys.wy;
	keys.odd = rotate(keys.odd, ulong2(17, 21));
	keys.odd ^= keys.zx;
	keys.z = rotate(keys.z, static_cast<ulong>(32));
}


// SipHash-2-4
inline uint sipHash24(ulong4 keys, const ulong nonce) {

	// Perform hash on keys
	keys.w ^= nonce;
	sipRound(keys);
	sipRound(keys);
	keys.even ^= ulong2(nonce, 255);
	sipRound(keys);
	sipRound(keys);
	sipRound(keys);
	sipRound(keys);
	keys.lo ^= keys.hi;

	// Check if edge bits is 32
	#if EDGE_BITS == 32
		// Return node from keys
		return keys.x ^ keys.y;
	// Otherwise
	#else
		// Return node from keys
		return (keys.x ^ keys.y) & NODE_MASK;
	#endif
}

struct SeedParams {
    uint64_t key0;
    uint64_t key1;
    uint64_t key2;
    uint64_t key3;
    uint nonceN;
};

// SipHash Metal kernel
kernel void build_seed(device ulong * result[[buffer(0)]],
			constant SeedParams &params [[buffer(1)]],
			constant uint &buffer_num [[buffer(2)]],
            uint id [[thread_position_in_grid]]) {

	ulong4 keys = {params.key0, params.key1, params.key2, params.key3};

	uint startIdx = id * 64;
    uint endIndex = startIdx + 64;

	// Incremental index for the seed is 1.
	ulong didx = (static_cast<ulong>(1)) << EDGE_BITS;

	ulong nonce = buffer_num * BUCKET_SIZE; // horiz offset
	nonce += (startIdx / BUCKET_SIZE) * BUCKETS_NUM * BUCKET_SIZE; // vert offset
	nonce += startIdx % BUCKET_SIZE; // offset inside the bucket
	nonce *= 2;
	nonce += params.nonceN;

    for (uint i = startIdx; i < endIndex; i+=1) {
        ulong hash = sipHash24(keys, nonce);

		hash |= didx;
		nonce += 2;
		//V2 nonce += BUCKETS_NUM*2;

		// now we need to write 5 byte hash (ELEMENT_SIZE is 5) into the result.

		size_t byte_pos = i * ELEMENT_SIZE;       // Compute byte position
	    size_t ulong_index = byte_pos / 8;        // Find which ulong to store in
	    size_t bit_offset  = (byte_pos % 8) * 8;  // Bit offset inside ulong

	    ulong mask = ((static_cast<ulong>(1) << (ELEMENT_SIZE*8)) - 1);

	    // Fetch existing value
	    ulong lower = result[ulong_index];

	    // Store in one or two ulong slots
	    lower &= ~(mask << bit_offset);  // Clear existing bits
	    lower |= (hash & mask) << bit_offset;  // Set new value
	    result[ulong_index] = lower;

	    // If hash crosses into the next ulong, store remaining bits
	    if (bit_offset > (64 - ELEMENT_SIZE*8)) {
	        ulong upper = result[ulong_index + 1];
	        upper &= ~(mask >> (64 - bit_offset));  // Clear existing bits
	        upper |= (hash & mask) >> (64 - bit_offset);  // Set new value
	        result[ulong_index + 1] = upper;
	    }
    }
}

// SipHash Metal kernel
kernel void nonce2hash(device uint * result[[buffer(0)]],
			constant SeedParams &params [[buffer(1)]],
            uint id [[thread_position_in_grid]]) {

	ulong4 keys = {params.key0, params.key1, params.key2, params.key3};

	uint startIdx = id * 16;
    uint endIndex = startIdx + 16;

    for (uint i = startIdx; i < endIndex; i+=1) {
        if (result[i]==0xffffffff)
            continue;
	    result[i] = sipHash24(keys, result[i] * 2 + params.nonceN);
    }
}

kernel void uv_to_nonces(device uint * U_table[[buffer(0)]],
            device uint * VV_table[[buffer(1)]],
            device uint * NN_table[[buffer(2)]],
			constant SeedParams &params [[buffer(3)]],
            uint id [[thread_position_in_grid]]) {

    ulong4 keys = {params.key0, params.key1, params.key2, params.key3};

	uint startIdx = id * 64;
    uint endIndex = startIdx + 64;

	ulong nonce = startIdx * 2;
	uint U_size = params.nonceN;

    for (uint i = startIdx; i < endIndex; i+=1) {
        uint u = sipHash24(keys, nonce);
        uint ui = u & 0x1;
        u &= 0xFFFFFFFE;
        uint32_t table_idx = u % U_size;

        for (uint i=0;i<U_size; i++, table_idx++) {
            uint32_t idx = table_idx % U_size;
            if (U_table[idx]==0) {
                break;
            }
            if (U_table[idx]==u) {
                // match for U, checking V hash
                uint v = sipHash24(keys, nonce+1);
                if (v == VV_table[idx*2 + ui]) {
                    // If v is matched - solution is found
                    NN_table[idx*2 + ui] = uint(nonce/2);
                }
            }
        }

		nonce += 2;
    }
}

    )";
}

// Init sip_hash calculations
static MTL::Library* init_sip_hash_library(MTL::Device* device, uint32_t EDGE_BITS, uint32_t ELEMENT_SIZE,
					uint32_t BUCKETS_NUM, uint32_t BUCKET_SIZE) {

	NS::Error* error = nullptr;
	MTL::Library* library = device->newLibrary(NS::String::string(generate_sip_hash_src( EDGE_BITS, ELEMENT_SIZE,
		BUCKETS_NUM, BUCKET_SIZE).c_str(),
				NS::UTF8StringEncoding), nullptr, &error);

	if (!library) {
		std::cerr << "Failed to compile Metal sip_hash: " << error->localizedDescription()->utf8String() << std::endl;
		exit(-1);
	}

	return library;
}


MetalOps::MetalOps(uint32_t _EDGE_BITS, uint32_t _ELEMENT_SIZE, uint32_t _BUCKET_BITS) :
		 EDGE_BITS(_EDGE_BITS), ELEMENT_SIZE(_ELEMENT_SIZE), BUCKET_BITS(_BUCKET_BITS)
{
    device = MTL::CreateSystemDefaultDevice();
    if (!device) {
        std::cerr << "Metal is not supported on this device!" << std::endl;
        exit(-1);
    }

	uint32_t bucket_size = 1 << (EDGE_BITS-BUCKET_BITS*2);

	sip_hash_library = init_sip_hash_library(device, EDGE_BITS, ELEMENT_SIZE, 1 << BUCKET_BITS, bucket_size);
}

MetalOps::~MetalOps() {
	sip_hash_library->release();
	device->release();
}

MTL::Buffer* MetalOps::create_buffer(uint buffer_size) {
    assert(buffer_size > 0);
    return device->newBuffer(buffer_size, MTL::ResourceStorageModeShared);
}

std::vector<MTL::Buffer*> MetalOps::allocate_buffers() {
	uint buffer_num = 1 << BUCKET_BITS;

	uint elements_per_buffer = (uint64_t(1) << EDGE_BITS) / buffer_num;

	std::vector<MTL::Buffer*> buffers(buffer_num);
	for (size_t i = 0; i < buffer_num; ++i) {
		buffers[i] = device->newBuffer(elements_per_buffer*ELEMENT_SIZE, MTL::ResourceStorageModeShared);
		if (!buffers[i]) {
			std::cerr << "Failed to allocate buffer " << i << std::endl;
			exit(-1);
		}
	}
	return buffers;
}

struct SeedParams {
    uint64_t key0;
    uint64_t key1;
    uint64_t key2;
    uint64_t key3;
    uint32_t nonceN;
};

// Calculate the SEED hashes for all buffers
void MetalOps::seed_hash(const std::vector<MTL::Buffer*>& buffers,
			const uint64_t v[4],
			uint nonceN ) {
	// Get function reference from the compiled library
	MTL::Function* function = sip_hash_library->newFunction(NS::String::string("build_seed", NS::UTF8StringEncoding));
	if (!function) {
		std::cerr << "Failed to find kernel function 'build_seed'!" << std::endl;
		exit(-1);
	}

	// Create Compute Pipeline
	NS::Error* error = nullptr;
	MTL::ComputePipelineState* pipelineState = device->newComputePipelineState(function, &error);
	if (!pipelineState) {
		std::cerr << "Failed to create compute pipeline: " << error->localizedDescription()->utf8String() << std::endl;
		exit(-1);
	}

	// Create Command Queue
	MTL::CommandQueue* commandQueue = device->newCommandQueue();

	// Create a list of command buffers for parallel execution
	std::vector<MTL::CommandBuffer*> commandBuffers;

    SeedParams params;
    params.key0 = v[0];
    params.key1 = v[1];
    params.key2 = v[2];
    params.key3 = v[3];
    params.nonceN = nonceN;

    MTL::Buffer* paramBuffer = device->newBuffer(&params, sizeof(SeedParams), MTL::ResourceStorageModeShared);

	for (uint32_t i = 0; i < buffers.size(); ++i) {
		// Unique nonce values per buffer

		// Command Buffer for parallel execution
		MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer();
		commandBuffers.push_back(commandBuffer);

		// Create Compute Encoder
		MTL::ComputeCommandEncoder* encoder = commandBuffer->computeCommandEncoder();
		encoder->setComputePipelineState(pipelineState);

		// Define keys (same for all, can be modified per buffer)

		// Pass small arguments via setBytes

		assert( buffers[i]->length() > ELEMENTS_PER_TASK * ELEMENT_SIZE );
		assert( buffers[i]->length() % (ELEMENTS_PER_TASK * ELEMENT_SIZE) == 0 );
		assert( ELEMENTS_PER_TASK * ELEMENT_SIZE % 64 == 0);

		uint tasks = buffers[i]->length() / (ELEMENTS_PER_TASK * ELEMENT_SIZE);

		// Bind buffer
		encoder->setBuffer(buffers[i], 0, 0);
        encoder->setBuffer(paramBuffer, 0, 1);
        encoder->setBytes(&i, sizeof(i), 2);

		// Dispatch threads for each buffer (assuming 1 result per buffer)
		encoder->dispatchThreads(MTL::Size(tasks, 1, 1), MTL::Size(1, 1, 1));
		encoder->endEncoding();

		// Submit task, no waiting
		commandBuffer->commit();
	}

	// wait for all to complete if needed
	for (auto* cmdBuffer : commandBuffers) {
		cmdBuffer->waitUntilCompleted();
		cmdBuffer->release();
	}

	// Buffers are not cleaned, we will need them for further processing

    // params are cleaned
	paramBuffer->release();

	commandQueue->release();
	pipelineState->release();
	function->release();
}

// Convert nonces into hashes
void MetalOps::nonces_to_hash(const std::vector<MTL::Buffer*>& buffers, const uint64_t v[4], uint nonceN) {
    // Get function reference from the compiled library
	MTL::Function* function = sip_hash_library->newFunction(NS::String::string("nonce2hash", NS::UTF8StringEncoding));
	if (!function) {
		std::cerr << "Failed to find kernel function 'nonce2hash'!" << std::endl;
		exit(-1);
	}

	// Create Compute Pipeline
	NS::Error* error = nullptr;
	MTL::ComputePipelineState* pipelineState = device->newComputePipelineState(function, &error);
	if (!pipelineState) {
		std::cerr << "Failed to create compute pipeline: " << error->localizedDescription()->utf8String() << std::endl;
		exit(-1);
	}

	// Create Command Queue
	MTL::CommandQueue* commandQueue = device->newCommandQueue();

	// Create a list of command buffers for parallel execution
	std::vector<MTL::CommandBuffer*> commandBuffers;

    SeedParams params;
    params.key0 = v[0];
    params.key1 = v[1];
    params.key2 = v[2];
    params.key3 = v[3];
    params.nonceN = nonceN;

    MTL::Buffer* paramBuffer = device->newBuffer(&params, sizeof(SeedParams), MTL::ResourceStorageModeShared);

	for (uint32_t i = 0; i < buffers.size(); ++i) {
		// Unique nonce values per buffer

		// Command Buffer for parallel execution
		MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer();
		commandBuffers.push_back(commandBuffer);

		// Create Compute Encoder
		MTL::ComputeCommandEncoder* encoder = commandBuffer->computeCommandEncoder();
		encoder->setComputePipelineState(pipelineState);

		// Define keys (same for all, can be modified per buffer)

		// Pass small arguments via setBytes

		assert( buffers[i]->length() > 0 );
		assert( buffers[i]->length() % 64 == 0 ); // processing 16 hashes per batch

		uint tasks = buffers[i]->length() / 64;

		// Bind buffer
		encoder->setBuffer(buffers[i], 0, 0);
        encoder->setBuffer(paramBuffer, 0, 1);

		// Dispatch threads for each buffer (assuming 1 result per buffer)
		encoder->dispatchThreads(MTL::Size(tasks, 1, 1), MTL::Size(1, 1, 1));
		encoder->endEncoding();

		// Submit task, no waiting
		commandBuffer->commit();
	}

	// wait for all to complete if needed
	for (auto* cmdBuffer : commandBuffers) {
		cmdBuffer->waitUntilCompleted();
		cmdBuffer->release();
	}

	// Buffers are not cleaned, we will need them for further processing

    // params are cleaned
	paramBuffer->release();

	commandQueue->release();
	pipelineState->release();
	function->release();
}

void MetalOps::uv_to_nonces(const uint64_t v[4], MTL::Buffer* U_hash_table, MTL::Buffer* VV_table, MTL::Buffer* NN_table) {
    // Get function reference from the compiled library
	MTL::Function* function = sip_hash_library->newFunction(NS::String::string("uv_to_nonces", NS::UTF8StringEncoding));
	if (!function) {
		std::cerr << "Failed to find kernel function 'uv_to_nonces'!" << std::endl;
		exit(-1);
	}

	// Create Compute Pipeline
	NS::Error* error = nullptr;
	MTL::ComputePipelineState* pipelineState = device->newComputePipelineState(function, &error);
	if (!pipelineState) {
		std::cerr << "Failed to create compute pipeline: " << error->localizedDescription()->utf8String() << std::endl;
		exit(-1);
	}

	// Create Command Queue
	MTL::CommandQueue* commandQueue = device->newCommandQueue();

    SeedParams params;
    params.key0 = v[0];
    params.key1 = v[1];
    params.key2 = v[2];
    params.key3 = v[3];
    params.nonceN = U_hash_table->length() / sizeof(uint32_t);

    MTL::Buffer* paramBuffer = device->newBuffer(&params, sizeof(SeedParams), MTL::ResourceStorageModeShared);

    // Command Buffer for parallel execution
	MTL::CommandBuffer* commandBuffer = commandQueue->commandBuffer();

	// Create Compute Encoder
	MTL::ComputeCommandEncoder* encoder = commandBuffer->computeCommandEncoder();
	encoder->setComputePipelineState(pipelineState);

	uint tasks = (uint64_t(1) << EDGE_BITS) / ELEMENTS_PER_TASK;

    // Bind buffer
	encoder->setBuffer(U_hash_table, 0, 0);
	encoder->setBuffer(VV_table, 0, 1);
	encoder->setBuffer(NN_table, 0, 2);
    encoder->setBuffer(paramBuffer, 0, 3);

	// Dispatch threads for each buffer (assuming 1 result per buffer)
	encoder->dispatchThreads(MTL::Size(tasks, 1, 1), MTL::Size(1, 1, 1));
	encoder->endEncoding();

	commandBuffer->commit();

	// wait for all to complete if needed
	commandBuffer->waitUntilCompleted();
    commandBuffer->release();

    // params are cleaned, other buffers are not
	paramBuffer->release();

	commandQueue->release();
	pipelineState->release();
	function->release();
}



