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

#include "../miner/metal.h"
#include "metal_test.h"
#include "../miner/data_storage.h"
#include "../miner/sip_hash.h"

void test_metal_4b() {
	const uint32_t EDGE_BITS = 12; // 128 elements per bucket
	const uint32_t ELEMENT_SIZE = 4;
	const uint32_t BUCKET_BITS = 2; // Buckets 4*4=16
	MetalOps metal(EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS);

	const uint32_t HASH_MASK = (1 << EDGE_BITS) - 1;

	std::vector<MTL::Buffer*> buffers = metal.allocate_buffers();
	assert(buffers.size()== 1<<BUCKET_BITS);

	for (uint32_t i = 0; i < buffers.size(); ++i) {
	    MTL::Buffer* buffer = buffers[i];
	    assert( buffer->length() == ELEMENT_SIZE* (1<<(EDGE_BITS-BUCKET_BITS)) );
	}

	uint64_t v[4] = { 0x736f6d6570736575ULL, 0x646f72616e646f6dULL,0x6c7967656e657261ULL,0x7465646279746573ULL};

	metal.seed_hash(buffers, v, 0);


	const int bucket_sz = 1<<(EDGE_BITS-BUCKET_BITS-BUCKET_BITS);
	const int bucket_num = 1 << BUCKET_BITS;
	// Validate what we get here
    for (uint32_t i = 0; i < buffers.size(); ++i) {
	    MTL::Buffer* buffer = buffers[i];

    	int nonce = i * bucket_sz * 2;

    	assert(ELEMENT_SIZE==4);
	    assert( buffer->length() == ELEMENT_SIZE*bucket_sz*bucket_num );
    	uint32_t * hashes = (uint32_t *) buffer->contents();

    	int k=0;
		for (uint32_t j = 0; j < bucket_sz*bucket_num; ++j) {
			uint32_t metalHash = hashes[j] & HASH_MASK;
			uint32_t calcHash = sip_hash(v, nonce) & HASH_MASK; // nonce & HASH_MASK;
			nonce += 2;
			k++;
			if (k==bucket_sz) {
				k=0;
				nonce+=bucket_sz*(bucket_num-1)*2;
			}
			assert( metalHash == calcHash );
		}
    }

	for (auto* buffer : buffers) {
		buffer->release();
	}
}

void test_metal_5b() {
	const uint32_t EDGE_BITS = 12; // 128 elements per bucket
	const uint32_t ELEMENT_SIZE = 5;
	const uint32_t BUCKET_BITS = 2; // Buckets 4*4=16
	MetalOps metal(EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS);

	const uint32_t HASH_MASK = (1 << EDGE_BITS) - 1;

	std::vector<MTL::Buffer*> buffers = metal.allocate_buffers();
	assert(buffers.size()== 1<<BUCKET_BITS);

	for (uint32_t i = 0; i < buffers.size(); ++i) {
		MTL::Buffer* buffer = buffers[i];
		assert( buffer->length() == ELEMENT_SIZE* (1<<(EDGE_BITS-BUCKET_BITS)) );
	}

	uint64_t v[4] = { 0x736f6d6570736575ULL, 0x646f72616e646f6dULL,0x6c7967656e657261ULL,0x7465646279746573ULL};

	metal.seed_hash(buffers, v, 0);


	const int bucket_sz = 1<<(EDGE_BITS-BUCKET_BITS-BUCKET_BITS);
	const int bucket_num = 1 << BUCKET_BITS;
	// Validate what we get here
	for (uint32_t i = 0; i < buffers.size(); ++i) {
		MTL::Buffer* buffer = buffers[i];

		int nonce = i * bucket_sz * 2;

		assert( buffer->length() == ELEMENT_SIZE*bucket_sz*bucket_num );
		void * hashes = buffer->contents();

		int k=0;
		for (uint32_t j = 0; j < bucket_sz*bucket_num; ++j) {
		    HashNonce hn = extract_hash_nonce<ELEMENT_SIZE, EDGE_BITS>((uint64_t *)hashes, j);
		    assert(hn.nonce == 1);
			uint32_t metalHash = hn.hash;
			uint32_t calcHash = sip_hash(v, nonce) & HASH_MASK; // nonce & HASH_MASK;
			nonce += 2;
			k++;
			if (k==bucket_sz) {
				k=0;
				nonce+=bucket_sz*(bucket_num-1)*2;
			}
			assert( metalHash == calcHash );
		}
	}

	for (auto* buffer : buffers) {
		buffer->release();
	}
}

const uint32_t SZ1 = 16;
const uint32_t SZ2 = 128;

static void init_bufs(std::vector<MTL::Buffer*> & buffers) {
    assert(buffers.size()==2);
    // Generate the nonces
    uint32_t * buf1 = (uint32_t *) buffers[0]->contents();
    for (uint i=0; i<SZ1; i++) {
        buf1[i] = 17*i + 3;
    }
    uint32_t * buf2 = (uint32_t *) buffers[1]->contents();
    for (uint i=0; i<SZ2; i++) {
        buf2[i] = 11*i + 100;
    }
}

void test_nonces_to_hash() {
    const uint32_t EDGE_BITS = 25; // 128 elements per bucket
    const uint32_t ELEMENT_SIZE = 5;
    const uint32_t BUCKET_BITS = 2; // Buckets 4*4=16
    MetalOps metal(EDGE_BITS, ELEMENT_SIZE, BUCKET_BITS);

    std::vector<MTL::Buffer*> buffers;
    buffers.push_back( metal.create_buffer(SZ1*4) );
    buffers.push_back( metal.create_buffer(SZ2*4) );

    uint64_t v[4] = { 0x736f6d6570736575ULL, 0x646f72616e646f6dULL,0x6c7967656e657261ULL,0x7465646279746573ULL};

    init_bufs(buffers);
    metal.nonces_to_hash(buffers, v, 1);

    uint32_t * buf1 = (uint32_t *) buffers[0]->contents();
    uint32_t * buf2 = (uint32_t *) buffers[1]->contents();
    const uint HASH_MASK = (1 << EDGE_BITS) - 1;

    for (uint i=0; i<SZ1; i++) {
        uint nonce = 17*i + 3;
        uint hash = sip_hash(v, nonce*2+1) & HASH_MASK;
        assert(buf1[i] == hash);
    }

    for (uint i=0; i<SZ2; i++) {
        uint nonce = 11*i + 100;
        uint hash = sip_hash(v, nonce*2+1) & HASH_MASK;
        assert(buf2[i] == hash);
    }

    init_bufs(buffers);
    metal.nonces_to_hash(buffers, v, 0);

    for (uint i=0; i<SZ1; i++) {
        uint nonce = 17*i + 3;
        uint hash = sip_hash(v, nonce*2+0) & HASH_MASK;
        assert(buf1[i] == hash);
    }

    for (uint i=0; i<SZ2; i++) {
        uint nonce = 11*i + 100;
        uint hash = sip_hash(v, nonce*2+0) & HASH_MASK;
        assert(buf2[i] == hash);
    }

    for (auto* buffer : buffers) {
        buffer->release();
    }
}


