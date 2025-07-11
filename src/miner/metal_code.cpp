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

#include "metal_code.h"
#include <assert.h>
#include "mem_pool.h"

#include "features.h"

std::string generate_metal_sources( uint32_t EDGE_BITS, uint32_t CYCLE_LEN, uint32_t BUCKETS_NUM, uint32_t BUCKETS_8B_NUM,
    const uint32_t EDGE_PER_BUCKET_PREV_PRIME, const uint32_t PHASE1_EDGE_PER_BUCKET,
    uint32_t COLLAPSE_STASH_SIZE, uint32_t COLLAPSE_STASH_COUNTER )
{
    assert((uint32_t(1) << EDGE_BITS) % BUCKETS_NUM == 0);
    const uint32_t EDGE_PER_BUCKET = (uint32_t(1) << EDGE_BITS) / BUCKETS_NUM;


    assert( EDGE_PER_BUCKET_PREV_PRIME <= EDGE_PER_BUCKET/2 );
    assert( EDGE_PER_BUCKET_PREV_PRIME > EDGE_PER_BUCKET/2 * 0.99 );
    assert(EDGE_PER_BUCKET_PREV_PRIME%2!=0);
    assert(EDGE_PER_BUCKET_PREV_PRIME%3!=0);
    assert(EDGE_PER_BUCKET_PREV_PRIME%5!=0);

	return std::string(R"(
#include <metal_stdlib>
using namespace metal;

	)") +

"#define EDGE_BITS " + std::to_string(EDGE_BITS) + "\n" +
"#define CYCLE_LEN " + std::to_string(CYCLE_LEN) + "\n" +
"#define BUCKETS_NUM " + std::to_string(BUCKETS_NUM) + "\n" +
"#define BUCKETS_8B_NUM " + std::to_string(BUCKETS_8B_NUM) + "\n" +
"#define MEM_BLOCK_SIZE " + std::to_string(MEM_POOL_UNITS) + "\n" +
"#define EDGE_PER_BUCKET " + std::to_string(EDGE_PER_BUCKET) + "\n" +
"#define EDGE_PER_BUCKET_PREV_PRIME " + std::to_string(EDGE_PER_BUCKET_PREV_PRIME) + "\n" +
"#define PHASE1_EDGE_PER_BUCKET " + std::to_string(PHASE1_EDGE_PER_BUCKET) + "\n" +
"#define COLLAPSE_STASH_SIZE " + std::to_string(COLLAPSE_STASH_SIZE) + "\n" +
"#define COLLAPSE_STASH_COUNTER " + std::to_string(COLLAPSE_STASH_COUNTER) + "\n" +

#ifdef USE_METRICS
"#define USE_METRICS\n" +
#endif

#ifdef TRACK_COLLAPSED
"#define TRACK_COLLAPSED\n" +
#endif

	R"(

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

    return keys.x ^ keys.y;
}

struct Step1Params {
    ulong key0;
    ulong key1;
    ulong key2;
    ulong key3;

    uint bucket_size;
    uint bucket_blocks[128];
};

// it is a kernel data that is responsible for:
// Hash keys
// Memory allocated chunks (every chunk size of 1/(BUCKET_NUM^2))
struct Step2Params {
    ulong key0;
    ulong key1;
    ulong key2;
    ulong key3;

    uint row_idx;

    uint in_block_size;
    uint out_block_size;

    // indexes of allocated blocks by bucket, high 16 bits are reserved for the block size
    uint in_blocks[128];

    // Indexes for the source data
    uint out_blocks[128];
};

// it is a kernel data for the trimming
// Has memory allocated buckets addresses.
struct TrimmingParams {
    uint mask_scale_k;
    uint mask_scale_mul;
    uint block_size;
    uint pass_bucket;

    // addresses of the blocks. Must be LAST item because of caller c++ code
    uint blocks[BUCKETS_NUM * BUCKETS_NUM];
};

struct DiscoverParams {
    ulong key0;
    ulong key1;
    ulong key2;
    ulong key3;

    uint hh_len;
    uint hh_k;
    uint nonce_sz;
    atomic_uint nonce_pos;

    bool upass;
};

#ifdef USE_METRICS
struct alignas(16) MetricsData {
    atomic_uint skipped_edges;
};
#endif

template <typename T>
inline device T * get_address(uint addr, device uint * buffer0, device uint * buffer1) {
    if (addr & 0x80000000) {
        return reinterpret_cast<device T*>(buffer1 + ulong((addr & 0x7FFFFFFF))*ulong(MEM_BLOCK_SIZE/4));
    }
    else {
        return reinterpret_cast<device T*>(buffer0 + ulong(addr)*ulong(MEM_BLOCK_SIZE/4));
    }
}

inline uint calc_bckt_index(bool uPass, uint row_idx, uint bucket ) {
    if (uPass)
        return row_idx*BUCKETS_NUM + bucket;
    else
        return bucket*BUCKETS_NUM + row_idx;
}

// step 1, shuffle nonces by the buckets. Number of expected buckets is BUCKETS_NUM.
// Every bucket will be used to build the mask after, then move the dat aback and forth
kernel void st1_build_buckets(device uint * buffer0[[buffer(0)]],
            device uint * buffer1[[buffer(1)]],
            constant Step1Params & params[[buffer(2)]],
#ifdef USE_METRICS
            device MetricsData * metric[[buffer(3)]],
#endif
            uint id [[thread_position_in_grid]],
            uint tid [[thread_index_in_threadgroup]],
            uint tgid [[threadgroup_position_in_grid]],
            uint tg_count [[threadgroups_per_grid]],
            uint threads_per_threadgroup [[threads_per_threadgroup]]) {

    const ulong4 keys = {params.key0, params.key1, params.key2, params.key3};
    const uint hash_mask = (uint(1) << EDGE_BITS) - 1;

    threadgroup atomic_uint group_indexes[BUCKETS_NUM];

    const uint job_size = params.bucket_size;
    const uint tg_limit = ulong(job_size) * (tgid+1) / tg_count;

    if (tid<BUCKETS_NUM) {
        const uint tg_offset = ulong(job_size) * tgid / tg_count;
        atomic_store_explicit( &group_indexes[tid], tg_offset, memory_order_relaxed);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    const uint tasks_per_threads = (uint(1) << EDGE_BITS) / tg_count / threads_per_threadgroup;
    const uint startNonce = id * tasks_per_threads;
    const uint endNonce = startNonce + tasks_per_threads;

    for (uint nonce=startNonce; nonce<endNonce; nonce++) {
        uint hash = sipHash24(keys, ulong(nonce)*2) & hash_mask;
        uint bucket_idx = hash / EDGE_PER_BUCKET;

        uint idx = atomic_fetch_add_explicit(&group_indexes[bucket_idx], 1, memory_order_relaxed);

        if ( idx < tg_limit ) {
            if (bucket_idx < BUCKETS_8B_NUM) {
                // 8 byte output...
                ulong res = ulong(hash) | ( ulong(nonce) << 32);
                get_address<ulong>(params.bucket_blocks[bucket_idx], buffer0, buffer1)[idx] = res;
            }
            else {
                // 4 byte nonce output
                get_address<uint>(params.bucket_blocks[bucket_idx], buffer0, buffer1)[idx] = nonce;
            }
        }
#ifdef USE_METRICS
        else {
            atomic_fetch_add_explicit(&metric->skipped_edges, 1, memory_order_relaxed);
        }
#endif
    }

   threadgroup_barrier(mem_flags::mem_threadgroup);

    // Fill with empty unused data
    const uint my_bucket = tid % BUCKETS_NUM;
    if (my_bucket<BUCKETS_8B_NUM) {
        device ulong * data = get_address<ulong>(params.bucket_blocks[my_bucket], buffer0, buffer1);
        while (true) {
            uint idx = atomic_fetch_add_explicit(&group_indexes[my_bucket], 8, memory_order_relaxed);
            if (idx >= tg_limit)
                break;

            for (int i=0; i<8 && idx+i<tg_limit; i++) {
                data[idx+i] = 0;
            }
            data[idx] = 0;
        }
    }
    else {
        device uint * data = get_address<uint>(params.bucket_blocks[my_bucket], buffer0, buffer1);
        while (true) {
            uint idx = atomic_fetch_add_explicit(&group_indexes[my_bucket], 16, memory_order_relaxed);
            if (idx >= tg_limit)
                break;

            for (int i=0; i<16 && idx+i<tg_limit; i++) {
                data[idx+i] = 0xFFFFFFFF;
            }
        }
    }
}

// TODO MAKE IT CONFIGURABLE
#define NONCE_TO_8B_RW_THREADS  1

// Decoding 1b data into 8b
// Currently 1b - it is a 4byte nonces. FF - empty data to skip
kernel void nonce_to_8b(device uint * nonce_data[[buffer(0)]],
                        device ulong * b8_data[[buffer(1)]],
                        constant Step1Params & params[[buffer(2)]],
                        uint id [[thread_position_in_grid]],
                        uint tid [[thread_index_in_threadgroup]],
                        uint tgid [[threadgroup_position_in_grid]],
                        uint tg_count [[threadgroups_per_grid]]) {

    const ulong4 keys = {params.key0, params.key1, params.key2, params.key3};
    const uint hash_mask = (uint(1) << EDGE_BITS) - 1;

    //const uint startIdx = id * 16;
    //const uint endIdx = startIdx + 16;

    threadgroup atomic_uint rw_pos[NONCE_TO_8B_RW_THREADS];

    if (tid<NONCE_TO_8B_RW_THREADS) {
        atomic_store_explicit( &rw_pos[tid], ulong(params.bucket_size) * tgid / tg_count, memory_order_relaxed);
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    const uint rw_pos_limit = ulong(params.bucket_size) * (tgid+1) / tg_count;
    const uint my_rw_idx = tid % NONCE_TO_8B_RW_THREADS;

    //for (uint i=startIdx; i<endIdx; i++) {
    while(true) {
        uint i = atomic_fetch_add_explicit(&rw_pos[my_rw_idx], 1, memory_order_relaxed);
        if (i>=rw_pos_limit)
            break;

        uint nonce = nonce_data[i];
        if (nonce==0xFFFFFFFF) {
            b8_data[i] = 0;
        }
        else {
            uint hash = sipHash24(keys, ulong(nonce)*2) & hash_mask;
            ulong res = ulong(hash) | ( ulong(nonce) << 32);
            b8_data[i] = res;
        }
    }
}

kernel void clean_data(device uint * mask[[buffer(0)]],
            uint id [[thread_position_in_grid]]) {

	uint startIdx = id * 64;
    //uint endIndex = startIdx + 64;

    device simd::ulong4 *mask_vec = reinterpret_cast<device simd::ulong4*>(&mask[startIdx]);
    mask_vec[0] = simd::ulong4(0);
    mask_vec[1] = simd::ulong4(0);
    mask_vec[2] = simd::ulong4(0);
    mask_vec[3] = simd::ulong4(0);
    mask_vec[4] = simd::ulong4(0);
    mask_vec[5] = simd::ulong4(0);
    mask_vec[6] = simd::ulong4(0);
    mask_vec[7] = simd::ulong4(0);
}

kernel void clean_data_FF(device uint * mask[[buffer(0)]],
            uint id [[thread_position_in_grid]]) {

	uint startIdx = id * 64;
    //uint endIndex = startIdx + 64;

    device simd::ulong4 *mask_vec = reinterpret_cast<device simd::ulong4*>(&mask[startIdx]);
    mask_vec[0] = simd::ulong4(0xFFFFFFFFFFFFFFFF);
    mask_vec[1] = simd::ulong4(0xFFFFFFFFFFFFFFFF);
    mask_vec[2] = simd::ulong4(0xFFFFFFFFFFFFFFFF);
    mask_vec[3] = simd::ulong4(0xFFFFFFFFFFFFFFFF);
    mask_vec[4] = simd::ulong4(0xFFFFFFFFFFFFFFFF);
    mask_vec[5] = simd::ulong4(0xFFFFFFFFFFFFFFFF);
    mask_vec[6] = simd::ulong4(0xFFFFFFFFFFFFFFFF);
    mask_vec[7] = simd::ulong4(0xFFFFFFFFFFFFFFFF);
}


kernel void trim_mask(device simd::uint4 * mask_even[[buffer(0)]],
            device simd::uint4 * mask_odd[[buffer(1)]],
            uint id [[thread_position_in_grid]]) {
   mask_even[id] &= mask_odd[id];
}

// TODO MAKE IT CONFIGURABLE
// So far 1 is the best
#define MASK_RW_THREADS  1

kernel void build_mask_simple(device uint * buffer0[[buffer(0)]],
                        device uint * buffer1[[buffer(1)]],
                        device atomic_uint * mask[[buffer(2)]],
                        constant Step2Params & params[[buffer(3)]],
                        uint id [[thread_position_in_grid]],
                        uint tid [[thread_index_in_threadgroup]],
                        uint tgid [[threadgroup_position_in_grid]]) {

    device ulong * data = get_address<ulong>(params.in_blocks[tgid], buffer0, buffer1);
    const uint block_len = params.in_block_size;

    const uint mask_sz_init = EDGE_PER_BUCKET/2;
    const uint mask_sz2 = (mask_sz_init+31)/32;

    threadgroup atomic_uint read_pos[MASK_RW_THREADS];

    if (tid<MASK_RW_THREADS) {
        atomic_store_explicit( &read_pos[tid], block_len * tid / MASK_RW_THREADS, memory_order_relaxed);
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    const uint read_idx = tid % MASK_RW_THREADS;
    const uint read_pos_limit = block_len * (read_idx+1) / MASK_RW_THREADS;

    // setting up the mask
    while(true) {
        uint i = atomic_fetch_add_explicit(&read_pos[read_idx], 1, memory_order_relaxed);
        if (i>=read_pos_limit)
            break;

        const uint64_t dt = data[i];
        if (!dt) {
            continue;
        }

        const uint hash = (uint) dt;
        const uint mask_bit_pos = (hash >> 1) % mask_sz_init;

        const uint mask_int = (mask_bit_pos / 32) + ( hash % 2 == 0 ? 0 : mask_sz2 );
        const uint mask_bit = 1 << ((mask_bit_pos) % 32);

        //sum++;

        uint prev = atomic_fetch_or_explicit(&mask[mask_int], mask_bit, memory_order_relaxed);
        if (prev & mask_bit) {
            atomic_fetch_or_explicit(&mask[mask_sz2*2 + mask_bit_pos / 32], mask_bit, memory_order_relaxed);
        }
    }
}

/*kernel void finish_secondary_mask(device atomic_uint * mask[[buffer(0)]],
                        device uint * mask2bits[[buffer(1)]],
                        device uint * mask2bitsPositions[[buffer(2)]],
                        uint id [[thread_position_in_grid]],
                        uint tid [[thread_index_in_threadgroup]],
                        uint tgid [[threadgroup_position_in_grid]]) {

    const uint mask2bitsPositionSize = EDGE_PER_BUCKET/2; // 1/ of full data size
    const uint mask2bitsPositionGrSz = mask2bitsPositionSize/BUCKETS_NUM;

    threadgroup atomic_uint read_pos[MASK_RW_THREADS];

    if (tid<MASK_RW_THREADS) {
        atomic_store_explicit( &read_pos[tid], mask2bitsPositionGrSz*tgid + mask2bitsPositionGrSz * tid / MASK_RW_THREADS, memory_order_relaxed);
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    const uint read_idx = tid % MASK_RW_THREADS;
    const uint read_pos_limit = mask2bitsPositions[tgid*MASK_RW_THREADS + read_idx];

    while(true) {
        uint i = atomic_fetch_add_explicit(&read_pos[read_idx], 1, memory_order_relaxed);

        if (i>=read_pos_limit)
            break;

        uint mask_bit_pos = mask2bits[i];

        const uint mask_int = (mask_bit_pos / 32);
        const uint mask_bit = 1 << ((mask_bit_pos) % 32);

        atomic_fetch_or_explicit(&mask[mask_int], mask_bit, memory_order_relaxed);
    }
}*/

inline uint to_bckt_index(uint ubkt, uint vbkt) {
    return ubkt * BUCKETS_NUM + vbkt;
}

kernel void build_mask(device uint * buffer0[[buffer(0)]],
                        device uint * buffer1[[buffer(1)]],
                        device atomic_uint * mask[[buffer(2)]],
                        constant TrimmingParams & params[[buffer(3)]],
                        constant uint * bucket_positions[[buffer(4)]],
                        uint id [[thread_position_in_grid]],
                        uint tid [[thread_index_in_threadgroup]],
                        uint tgid [[threadgroup_position_in_grid]]) {

    // We want to scale
    const uint mask_scale_k = params.mask_scale_k;
    const uint mask_mul = params.mask_scale_mul;
    const bool passU = (params.pass_bucket & 0x01)!=0;

    const uint row_idx = (params.pass_bucket >> 4) & 0xFF;

    const uint bkt_idx = calc_bckt_index(passU, row_idx, tgid);
    device ulong * data = get_address<ulong>( params.blocks[bkt_idx], buffer0, buffer1);
    uint block_edge_num = bucket_positions[bkt_idx];

    const uint mask_sz_init = EDGE_PER_BUCKET/2;
    const uint mask_sz_final = mask_sz_init / mask_scale_k;

    const uint mask_sz2 = (mask_sz_init+31)/32;

    const uint UBit = uint(1) << 31;

    threadgroup atomic_uint read_pos[MASK_RW_THREADS];

    if (tid<MASK_RW_THREADS) {
        const uint idx = tid;
        atomic_store_explicit( &read_pos[idx], block_edge_num * idx / MASK_RW_THREADS, memory_order_relaxed);
    }

    const uint read_idx = tid % MASK_RW_THREADS;
    const uint read_pos_limit = block_edge_num * (read_idx+1) / MASK_RW_THREADS;

    threadgroup_barrier(mem_flags::mem_threadgroup);

    // setting up the mask
    while(true) {
            uint i = atomic_fetch_add_explicit(&read_pos[read_idx], 1, memory_order_relaxed);
            if (i>=read_pos_limit)
                break;

            const ulong dt = data[i];
            if (!dt) {
                continue;
            }

            uint hash;
            if (passU) {
                hash = (uint) dt;
                if ((hash & UBit) == 0)
                    continue; // V - V pair

                hash ^= UBit;
            }
            else {
                hash = (uint) (dt>>32);
                if ((hash & UBit) != 0)
                    continue; // U - U pair
            }

            uint mask_bit_pos = (hash >> 1) % mask_sz_init;

            mask_bit_pos =  ((uint)((ulong(mask_bit_pos)*mask_mul) % EDGE_PER_BUCKET_PREV_PRIME)) % mask_sz_final;

            uint mask_int = (mask_bit_pos / 32) + ( hash % 2 == 0 ? 0 : mask_sz2 );
            uint mask_bit = 1 << ((mask_bit_pos) % 32);

            uint prev = atomic_fetch_or_explicit(&mask[mask_int], mask_bit, memory_order_relaxed);
            if (prev & mask_bit) {
                const uint mask2_int = (mask_bit_pos / 32) + mask_sz2*2;
                atomic_fetch_or_explicit(&mask[mask2_int], mask_bit, memory_order_relaxed);
            }

            // Checking if it is a pair, the second item
            if (passU) {
                hash = (uint) (dt>>32);
                if ((hash & UBit) == 0)
                    continue; // Normal case

                // Got U - U pair
                hash ^= UBit;
            }
            else {
                hash = (uint) dt;
                if ((hash & UBit) != 0)
                    continue; // Normal case
                // got V - V pair
            }

            uint secondary_bucket = hash / EDGE_PER_BUCKET;
            if (secondary_bucket!=row_idx)
                continue; // skipping here. Will move during next trimming step

            mask_bit_pos = (hash >> 1) % mask_sz_init;

            mask_bit_pos =  ((uint)((ulong(mask_bit_pos)*mask_mul) % EDGE_PER_BUCKET_PREV_PRIME)) % mask_sz_final;

            mask_int = (mask_bit_pos / 32) + ( hash % 2 == 0 ? 0 : mask_sz2 );
            mask_bit = 1 << ((mask_bit_pos) % 32);

            prev = atomic_fetch_or_explicit(&mask[mask_int], mask_bit, memory_order_relaxed);
            if (prev & mask_bit) {
                const uint mask2_int = (mask_bit_pos / 32) + mask_sz2*2;
                atomic_fetch_or_explicit(&mask[mask2_int], mask_bit, memory_order_relaxed);
            }
    }
}

// Splitting 8b data into buckets.
// arg1 - resulting buckets bits
// arg2 - shift for hash to cacluate the mask index.
kernel void trimmed_to_next_buckets_st2(device uint * buffer0[[buffer(0)]],
                        device uint * buffer1[[buffer(1)]],
                        constant uint * mask[[buffer(2)]],
                        constant Step2Params & params[[buffer(3)]],
                        device uint * collapse_edges[[buffer(4)]],
                        device atomic_uint * collapse_count[[buffer(5)]], // Can be null if TRACK_COLLAPSED not defined!
                        device uint * bucket_positions[[buffer(6)]],
#ifdef USE_METRICS
                        device MetricsData & metric[[buffer(7)]],
#endif
                        uint id [[thread_position_in_grid]],
                        uint tid [[thread_index_in_threadgroup]],
                        uint tgid [[threadgroup_position_in_grid]],
                        uint tg_count [[threadgroups_per_grid]])
{
    const ulong4 keys = {params.key0, params.key1, params.key2, params.key3};

    device ulong * data = get_address<ulong>(params.in_blocks[tgid], buffer0, buffer1);
    const uint row_idx = params.row_idx;
    const uint in_items_per_block = params.in_block_size;

    const uint out_items_per_bucket = params.out_block_size;
    const uint b8_limit = uint( ulong(out_items_per_bucket) * ulong(tgid+1) / tg_count );

    const uint mask_sz_init = EDGE_PER_BUCKET/2;
    const uint mask_chunk_sz = mask_sz_init/32;
    const uint mask_secondary_offset = mask_chunk_sz*2;

    const uint sip_hash_mask = (uint(1) << EDGE_BITS) - 1;

    // TODO ADD INTO CONFIG. So far 1 is the best
    const uint TRIM_ST2_READ_THREADS = 1;
    threadgroup atomic_uint read_pos[TRIM_ST2_READ_THREADS];

    threadgroup atomic_uint group_indexes[BUCKETS_NUM];

    if (tid<TRIM_ST2_READ_THREADS) {
        atomic_store_explicit( &read_pos[tid], out_items_per_bucket * tid / TRIM_ST2_READ_THREADS, memory_order_relaxed);
    }

    if (tid<BUCKETS_NUM) {
        const uint b8_offset = uint( ulong(out_items_per_bucket) * ulong(tgid) / tg_count );
        atomic_store_explicit( &group_indexes[tid], b8_offset, memory_order_relaxed);
    }

    // !!!!!! ID, not tid. Just need to update bucket_positions once
    if (id>=BUCKETS_NUM*2 && id<BUCKETS_NUM*3) {
        const uint idx = id - BUCKETS_NUM*2;
        bucket_positions[row_idx*BUCKETS_NUM + idx] = out_items_per_bucket;
    }


    threadgroup_barrier(mem_flags::mem_threadgroup);

    const uint read_idx = tid % TRIM_ST2_READ_THREADS;
    const uint read_pos_limit = in_items_per_block * (read_idx+1) / TRIM_ST2_READ_THREADS;

    // Note, reading both primary and secondary masks are NOT OK. 19 ms - both masks. 11ms - any of them
    //for (uint i = startIdx; i < endIdx; i++) {
    while(true) {
        uint i = atomic_fetch_add_explicit(&read_pos[read_idx], 1, memory_order_relaxed);
        if (i>=read_pos_limit)
            break;

        // Hash is a low int of the data.
        const ulong dt = data[i];
        if (!dt)
            continue;

        uint hash1 = (uint) dt;

        uint mask_bit_pos = (hash1>>1) % mask_sz_init;

        uint mask_int = (mask_bit_pos / 32);
        uint mask_bit = 1 << ((mask_bit_pos) % 32);

        if ( (mask[mask_int] & mask_bit) == 0) {
            continue;
        }

        if ( (mask[mask_secondary_offset+mask_int] & mask_bit) == 0) {
            uint hash2 = (uint) (dt >> 32);
            hash2 = sipHash24(keys, ulong(hash2)*2+1) & sip_hash_mask;

            collapse_edges[mask_bit_pos*2 + hash1 % 2] = hash2;
            // no updated for counts because count is 0. No needs for that
#ifdef TRACK_COLLAPSED
            atomic_fetch_add_explicit( &collapse_count[ mask_bit_pos/4 ], uint(1) << (8 * (mask_bit_pos % 4)), memory_order_relaxed);
#endif
        }
        else {
            // Regular U-V hash
            uint hash2 = (uint) (dt >> 32);
            hash2 = sipHash24(keys, ulong(hash2)*2+1) & sip_hash_mask;

            uint bucket_idx = hash2 / EDGE_PER_BUCKET;

            uint b8_idx = atomic_fetch_add_explicit(&group_indexes[bucket_idx], 1, memory_order_relaxed);

            if ( b8_idx < b8_limit ) {
#ifdef TRACK_COLLAPSED
                // let's first hash contain the collapse index. It is 0 for emitted edges
                get_address<ulong>( params.out_blocks[bucket_idx], buffer0, buffer1 )[b8_idx] = ulong(hash1 % EDGE_PER_BUCKET) | (ulong(1)<<31) | (ulong(hash2) << 32);
#else
                get_address<ulong>( params.out_blocks[bucket_idx], buffer0, buffer1 )[b8_idx] = ulong(hash1) | (ulong(1)<<31) | (ulong(hash2) << 32);
#endif
            }
#ifdef USE_METRICS
            else {
                 atomic_fetch_add_explicit(&metric.skipped_edges, 1, memory_order_relaxed);
            }
#endif
        }
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Fill the rest with empty unused data, expected that small amount of the data will be used
    const uint my_bucket = tid % BUCKETS_NUM;
    device ulong * fill_data = get_address<ulong>( params.out_blocks[my_bucket], buffer0, buffer1 );
    while (true) {
        uint b8_idx = atomic_fetch_add_explicit(&group_indexes[my_bucket], 16, memory_order_relaxed);
        if (b8_idx >= b8_limit)
            break;
        for (int i=0; i<16 && b8_idx+i<b8_limit; i++) {
            fill_data[b8_idx+i] = 0;
        }
    }
}

kernel void compact_zeroes(device uint * buffer0[[buffer(0)]],
                        device uint * buffer1[[buffer(1)]],
                        constant TrimmingParams & params[[buffer(2)]],
                        device uint * bucket_positions[[buffer(3)]],
                        uint id [[thread_position_in_grid]],
                        uint tid [[thread_index_in_threadgroup]],
                        uint tgid [[threadgroup_position_in_grid]],
                        uint tg_count [[threadgroups_per_grid]])
{
    const bool uPass = (params.pass_bucket & 1)!=0;
    const uint row_idx = (params.pass_bucket >> 4) & 0xFF;

    // TODO ADD INTO CONFIG
    const uint COMPACT_READ_THREADS = 1;
    threadgroup atomic_int read_pos_left[COMPACT_READ_THREADS];
    threadgroup atomic_int read_pos_right[COMPACT_READ_THREADS];
    threadgroup atomic_int resulting_max_pos;
    threadgroup atomic_int resulting_min_pos;

    uint bkt_index = calc_bckt_index(uPass, row_idx, tgid);
    device ulong * data = get_address<ulong>(params.blocks[bkt_index], buffer0, buffer1);

    int rightPosLimit = bucket_positions[bkt_index];
    rightPosLimit--; // We don't want next, we need current


    if (tid<COMPACT_READ_THREADS) {
        atomic_store_explicit( &read_pos_left[tid], tid, memory_order_relaxed);
    }
    else if (tid<COMPACT_READ_THREADS*2) {
        uint idx = tid - COMPACT_READ_THREADS;
        int rpos = rightPosLimit - (rightPosLimit % COMPACT_READ_THREADS) + idx;
        if (rpos>rightPosLimit) {
            rpos -= COMPACT_READ_THREADS;
        }
        atomic_store_explicit( &read_pos_right[idx], rpos, memory_order_relaxed);
    }
    else if (tid<COMPACT_READ_THREADS*2 + 1) {
        atomic_store_explicit( &resulting_max_pos, 0, memory_order_relaxed);
    }
    else if (tid<COMPACT_READ_THREADS*2 + 2) {
        atomic_store_explicit( &resulting_min_pos, 0x7FFFFFFF, memory_order_relaxed);
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    const uint atomic_idx = tid % COMPACT_READ_THREADS;

    // Now we want to do some compaction
    int rightIdx = rightPosLimit;
    while(true) {
        int leftIdx = atomic_fetch_add_explicit(&read_pos_left[atomic_idx], COMPACT_READ_THREADS, memory_order_relaxed);
        while(leftIdx < rightIdx) {
            if (!data[leftIdx]) // found zero
                break;
            leftIdx = atomic_fetch_add_explicit(&read_pos_left[atomic_idx], COMPACT_READ_THREADS, memory_order_relaxed);
        }

        rightIdx = atomic_fetch_add_explicit(&read_pos_right[atomic_idx], int(-COMPACT_READ_THREADS), memory_order_relaxed);
        while (rightIdx >= leftIdx) {
            if (data[rightIdx]) // found non zero
                break;
            rightIdx = atomic_fetch_add_explicit(&read_pos_right[atomic_idx], int(-COMPACT_READ_THREADS), memory_order_relaxed);
        }

        if (leftIdx < rightIdx) {
            data[leftIdx] = data[rightIdx];
            data[rightIdx] = 0;
        }
        else {
            atomic_fetch_min_explicit(&resulting_min_pos, max(0,min(leftIdx, rightIdx)), memory_order_relaxed);
            if (leftIdx > rightPosLimit) {
                leftIdx = rightPosLimit+1;
            }
            atomic_fetch_max_explicit(&resulting_max_pos, leftIdx, memory_order_relaxed);
            break;
        }
    }

    // preparing for the second step. According our test we never had more then 4000 spread per bucket. Let's compact that as well
    // 2000 points should be enough to compress extra zeroes
    const uint END_BUFFER_SIZE = 2000;
    threadgroup atomic_int end_buffer_pos;
    threadgroup ulong end_buffer[END_BUFFER_SIZE];

    if (tid==0) {
        atomic_store_explicit( &end_buffer_pos, 0, memory_order_relaxed);
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    // read from resulting_max_pos back to the resulting_min_pos
    int read_min_pos = atomic_load_explicit(&resulting_min_pos, memory_order_relaxed);

    while(true) {
        int idx = atomic_fetch_add_explicit(&resulting_max_pos, -1, memory_order_relaxed);
        if (idx<read_min_pos)
            break;
        ulong dt = data[idx];
        if (!dt) // skip 0
            continue;
        // dt is a data, allocating end_buffer position
        int eb_idx = atomic_fetch_add_explicit(&end_buffer_pos, 1, memory_order_relaxed);
        if (eb_idx>=int(END_BUFFER_SIZE))
            break;
        end_buffer[eb_idx] = dt;
        data[idx] = 0;
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    while(true) {
        // -1 because end_buffer_pos is one index ahead (fetch first, operation after), compensating that
        int eb_idx = atomic_fetch_add_explicit(&end_buffer_pos, -1, memory_order_relaxed) - 1;
        if (eb_idx>=int(END_BUFFER_SIZE))
            continue;
        if (eb_idx<0)
            break;

        ulong dt = end_buffer[eb_idx];
        // Search for a location to insert
        while(true) {
            int dt_idx = atomic_fetch_add_explicit(&resulting_min_pos, 1, memory_order_relaxed);
            if (dt_idx>rightPosLimit)
                break;
            if (data[dt_idx])
                continue;
            data[dt_idx] = dt;
            break;
        }
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    // Updating bucket_positions with a new value.
    // Note, there are probably some zeroes becuase of race conditions, but it should be not big number, so we can accept that
    if (tid==0) {
        int max_pos = max(0, int(atomic_load_explicit(&resulting_max_pos, memory_order_relaxed)) + 1);
        int min_pos = atomic_load_explicit(&resulting_min_pos, memory_order_relaxed);
        bucket_positions[ bkt_index ] = max( min_pos, max_pos );
    }
}

kernel void apply_collapsed_data(device uint * buffer0[[buffer(0)]],
                        device uint * buffer1[[buffer(1)]],
                        device ulong * mask[[buffer(2)]],
                        constant TrimmingParams & params[[buffer(3)]],
                        device atomic_uint * bucket_positions[[buffer(4)]],
                        device uint * collapse_edges[[buffer(5)]],
                        device uint * collapse_count[[buffer(6)]], // Can be null if TRACK_COLLAPSED not defined!
#ifdef USE_METRICS
                        device MetricsData & metric[[buffer(7)]],
#endif
                        uint id [[thread_position_in_grid]])
{

    // indexes are for mask uint items
    const uint mask_sz_init = EDGE_PER_BUCKET/2;
    const uint mask_chunk_sz = mask_sz_init/64;
    const uint mask_odd_offset = mask_chunk_sz;
    const uint mask_secondary_offset = mask_chunk_sz*2;

    const uint idx = id;

    ulong collapsed_mask = mask[idx] & (~mask[mask_secondary_offset + idx]);
    // we can reset the mask here, so no
    mask[idx] = 0;
    mask[mask_secondary_offset + idx] = 0;
    mask[mask_odd_offset + idx] = 0;

    if (collapsed_mask != 0) {
        const bool passU = (params.pass_bucket & 0x01)!=0;
        const uint row_idx = (params.pass_bucket >> 4) & 0xFF;
        const uint UBit = uint(1) << 31;

#ifdef TRACK_COLLAPSED
        uint need_reset_counters = 0;
#endif
        for (uint i=0; i<64; i++) {
            if (collapsed_mask & (ulong(1) << i)) {
                int base_idx = (idx * 64 + i)*2;
                uint hash1 = collapse_edges[base_idx];
                uint hash2 = collapse_edges[base_idx+1];
                collapse_edges[base_idx] = 0xFFFFFFFF; // restoring data back
                collapse_edges[base_idx+1] = 0xFFFFFFFF;

                // need to move now
                bool u1 = (hash1 & UBit)!=0;
                bool u2 = (hash2 & UBit)!=0;

                uint ubkt;
                uint vbkt;

                // we want to keep it in normal order for U/V.
                if (u1 != u2) {
                    if (!u1) {
                        // swapping hashes
                        uint tmp = hash1;
                        hash1 = hash2;
                        hash2 = tmp;
                    }
                    ubkt = (hash1 & 0x7FFFFFFF) / EDGE_PER_BUCKET;
                    vbkt = hash2 / EDGE_PER_BUCKET;
                }
                else {
                    // U-U or V-V case
                    if (u1) {
                        // U-U
                        if (hash1>hash2) {
                            uint tmp = hash1;
                            hash1 = hash2;
                            hash2 = tmp;
                        }
                        // hash1 < hash2, hash1 is a target
                        if (passU) {
                            uint ub1 = (hash1 & 0x7FFFFFFF) / EDGE_PER_BUCKET;
                            uint ub2 = (hash2 & 0x7FFFFFFF) / EDGE_PER_BUCKET;

                            if (ub1 > row_idx) {
                                ubkt = ub1;
                            }
                            else if (ub2>row_idx) {
                                ubkt = ub2;
                                // also swapping hashes
                                uint tmp = hash1;
                                hash1 = hash2;
                                hash2 = tmp;
                            }
                            else {
                                // both are late, using min bucket for the next iteration
                                ubkt = ub1;
                            }
                        }
                        else {
                            // V pass will skip this edge, so moving into the low U bucket
                            ubkt = (hash1 & 0x7FFFFFFF) / EDGE_PER_BUCKET;
                        }

                        vbkt = row_idx; // Can be any bucket, lets use current
                    }
                    else {
                        // V-V, second must be smaller
                        if (hash2>hash1) {
                            uint tmp = hash1;
                            hash1 = hash2;
                            hash2 = tmp;
                        }
                        // hash2 < hash1. hash2 is a target
                        ubkt = row_idx; // Can be any bucket, let use current

                        if (passU) {
                            // U pass will skip this edge, so moving into the low V bucket
                            // Lowest for V is the second one
                            vbkt = hash2 / EDGE_PER_BUCKET;
                        }
                        else {
                            uint vb1 = (hash2 & 0x7FFFFFFF) / EDGE_PER_BUCKET;
                            uint vb2 = (hash1 & 0x7FFFFFFF) / EDGE_PER_BUCKET;

                            if (vb1 > row_idx) {
                                vbkt = vb1;
                            }
                            else if (vb2>row_idx) {
                                vbkt = vb2;
                                // also swapping hashes
                                uint tmp = hash1;
                                hash1 = hash2;
                                hash2 = tmp;
                            }
                            else {
                                // both are late, using min bucket for the next iteration
                                vbkt = vb1;
                            }
                        }
                    }
                }

#ifdef TRACK_COLLAPSED
                // collapse_count is written in every byte.
                // -2 because counters are -1 based.
                // +1 to cover this one remove point.
                int clps_value = int( ( collapse_count[base_idx/8] >> (8 * ((base_idx/2) % 4)) ) & 0xFF) - 1; // -2 + 1

                need_reset_counters |= 1 << (i/4);

                if (clps_value<0 || clps_value > CYCLE_LEN-1)
                    continue; // It is a very long path, cycle will be longer than 41 (CYCLE_LEN-1). We can delete it by skipping

                if (hash1/2 == hash2/2) {
                    // found a perfect cycle.
                    if (clps_value != CYCLE_LEN-1)   // cycle 41 - is a solution
                        continue;
                }

                ulong edge_data;
                if (u1==false && u2==false) {
                    // V-V case, Hash2 should have index data
                    edge_data = (ulong(hash2 % EDGE_PER_BUCKET + uint(clps_value) * EDGE_PER_BUCKET ) << 32) | ulong(hash1);
                }
                else {
                    const uint h1U = hash1 & 0x80000000;
                    edge_data = (ulong(hash2) << 32) | h1U | ulong((hash1 & 0x7FFFFFFF) % EDGE_PER_BUCKET + uint(clps_value) * EDGE_PER_BUCKET );
                }
#else
                // u/v bkt are known
                ulong edge_data = (ulong(hash2) << 32) | ulong(hash1);
#endif

                uint bkt_idx = ubkt*BUCKETS_NUM + vbkt;
                uint move_idx = atomic_fetch_add_explicit( &bucket_positions[bkt_idx], 1, memory_order_relaxed );

                if (move_idx < params.block_size) {
                    device ulong * bkt_data = get_address<ulong>( params.blocks[bkt_idx], buffer0, buffer1);
                    bkt_data[move_idx] = edge_data;
                }
#ifdef USE_METRICS
                else {
                    atomic_fetch_add_explicit(&metric.skipped_edges, 1, memory_order_relaxed);
                }
#endif
            }
        }
#ifdef TRACK_COLLAPSED
        const uint idx0 = idx*64/4;
        for (int r=0;r<16; r++) {
            if (need_reset_counters & (1<<r) ) {
                collapse_count[idx0 + r] = 0;
            }
        }
#endif
    }
}

inline void stash_collapsed(uint collapse_counter, uint hash2keep, uint hash_meet_point, device atomic_uint * collapse_stash) {
    if ( collapse_counter < COLLAPSE_STASH_COUNTER)
        return;

    uint stash_idx = atomic_fetch_add_explicit( &collapse_stash[0], 1, memory_order_relaxed);
    if (stash_idx < COLLAPSE_STASH_SIZE) {
        atomic_store_explicit( &collapse_stash[1 + stash_idx*2],
                (hash2keep & 0x07FFFFFF) | ( uint(collapse_counter-COLLAPSE_STASH_COUNTER) << (32-5) ), memory_order_relaxed);
        atomic_store_explicit( &collapse_stash[1 + stash_idx*2+1], hash_meet_point, memory_order_relaxed);
    }
}


// Splitting 8b data into buckets.
// arg1 - resulting buckets bits
// arg2 - shift for hash to cacluate the mask index.
kernel void trimmed_to_next_buckets_st_trim(device uint * buffer0[[buffer(0)]],
                        device uint * buffer1[[buffer(1)]],
                        device uint * mask[[buffer(2)]],
                        constant TrimmingParams & params[[buffer(3)]],
                        device atomic_uint * bucket_positions[[buffer(4)]],
                        device uint * collapse_edges[[buffer(5)]],
                        device atomic_uint * collapse_count[[buffer(6)]], // Can be null if TRACK_COLLAPSED not defined!
                        device atomic_uint * collapse_stash[[buffer(7)]],
#ifdef USE_METRICS
                        device MetricsData & metric[[buffer(8)]],
#endif
                        uint id [[thread_position_in_grid]],
                        uint tid [[thread_index_in_threadgroup]],
                        uint tgid [[threadgroup_position_in_grid]],
                        uint tg_count [[threadgroups_per_grid]])
{
    const uint mask_scale_k = params.mask_scale_k;
    const uint mask_mul = params.mask_scale_mul;
    const bool passU = (params.pass_bucket & 0x01)!=0;

    const uint row_idx = (params.pass_bucket >> 4) & 0xFF;

    const uint bkt_idx = calc_bckt_index(passU, row_idx, tgid);

    device ulong * data = get_address<ulong>( params.blocks[bkt_idx], buffer0, buffer1);
    uint edges_in_block = atomic_load_explicit( &bucket_positions[bkt_idx], memory_order_relaxed);

    const uint mask_sz_init = EDGE_PER_BUCKET/2;
    const uint mask_chunk_sz = mask_sz_init/32;
    const uint mask_secondary_offset = mask_chunk_sz*2;
    const uint mask_sz_final = mask_sz_init / mask_scale_k;

    const uint UBit = uint(1) << 31;
    const uint UMask = UBit-1;

    // TODO MOVE INTO PARAMS
    const uint TRIM_READ_THREADS = 1;
    threadgroup atomic_uint read_pos[TRIM_READ_THREADS];

    if (tid<TRIM_READ_THREADS) {
        atomic_store_explicit( &read_pos[tid], edges_in_block * tid / TRIM_READ_THREADS, memory_order_relaxed);
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    const uint read_idx = tid % TRIM_READ_THREADS;
    const uint read_pos_limit = edges_in_block * (read_idx+1) / TRIM_READ_THREADS;

    while(true) {
        uint idx = atomic_fetch_add_explicit(&read_pos[read_idx], 1, memory_order_relaxed);
        if (idx>=read_pos_limit)
                break;

        // Hash is a low int of the data.
        const ulong dt = data[idx];
        if (!dt) {
            continue;
        }

        const uint hash1 = (uint) dt;
        const uint hash2 = (uint) (dt>>32);

        const bool u1 = (hash1 & UBit) != 0;
        const bool u2 = (hash2 & UBit) != 0;

#ifdef TRACK_COLLAPSED
        // hash1 has collapse counter, let's retrieve it
        uint collapse_counter;
        if (u1==false && u2==false) {
            collapse_counter = (hash2 & UMask) / EDGE_PER_BUCKET;
        }
        else {
            collapse_counter = (hash1 & UMask) / EDGE_PER_BUCKET;
        }
#endif

        if (u1!=u2) {
            // U-V case
            const uint hash = passU ? (hash1 & UMask) : hash2;

            uint mask_bit_pos = (hash>>1) % mask_sz_init;
            mask_bit_pos =  ((uint)((ulong(mask_bit_pos)*mask_mul) % EDGE_PER_BUCKET_PREV_PRIME)) % mask_sz_final;

            const uint mask_int = (mask_bit_pos / 32);
            const uint mask_bit = 1 << ((mask_bit_pos) % 32);

            if ((mask[mask_int] & mask_bit) == 0) {
                data[idx] = 0; // inplace, resetting to zero, will compact after
                continue;
            }

            // Check if can be collapsed this hash
            if ((mask[mask_secondary_offset+mask_int] & mask_bit) == 0) {
                // collapsed will be resetted in any case. Since we merge opposite edges, not likely result will be in this bucket
                data[idx] = 0;

#ifdef TRACK_COLLAPSED
                // hash1 is U because it is U/V case
                const uint hash1mask = (hash1 & UMask) % EDGE_PER_BUCKET;
                const uint rec_hash1 = UBit | (uint(EDGE_PER_BUCKET) * tgid + hash1mask);
                atomic_fetch_add_explicit( &collapse_count[ mask_bit_pos/4 ], (collapse_counter+1) << (8 * (mask_bit_pos % 4)), memory_order_relaxed);
#endif

                if (passU) {
                    // Hash2 has a full info, including the bucket
                    collapse_edges[mask_bit_pos*2 + hash%2] = hash2;
#ifdef TRACK_COLLAPSED
                    stash_collapsed(collapse_counter, hash2, rec_hash1, collapse_stash);
#endif
                }
                else
                {
                    // hash1 doesn't have bucket info, we need to put it back.
#ifdef TRACK_COLLAPSED
                    stash_collapsed(collapse_counter, rec_hash1, hash2, collapse_stash);
                    collapse_edges[mask_bit_pos*2 + hash%2] = rec_hash1;
#else
                    collapse_edges[mask_bit_pos*2 + hash%2] = hash1;
#endif
                }
                continue;
            }
        }
        else {
            // U-U or V-V case
            if (u1 != passU)
                continue;

            uint h1,h2, case_UBit;
            if (passU) {
                h1 = hash1 & UMask;
                h2 = hash2 & UMask;
                case_UBit = UBit;
            }
            else {
                h1 = hash2;
                h2 = hash1;
                case_UBit = 0;
            }

            // h1 is expected to be in the current bucket, checking if both hashes are - it is an edge case
            // With TRACK_COLLAPSED we can't validate that.

            uint secondary_bucket_idx = h2 / EDGE_PER_BUCKET;
            if (secondary_bucket_idx != row_idx) {
                // Normal use case for U-U or V-V

                uint mask_bit_pos = (h1>>1) % mask_sz_init;
                mask_bit_pos =  ((uint)((ulong(mask_bit_pos)*mask_mul) % EDGE_PER_BUCKET_PREV_PRIME)) % mask_sz_final;

                const uint mask_int = (mask_bit_pos / 32);
                const uint mask_bit = 1 << ((mask_bit_pos) % 32);

                data[idx] = 0; // inplace, resetting to zero, will compact after. Reset happens in ALL branches

                if ((mask[mask_int] & mask_bit) == 0) {
                    continue;
                }

                // Check if can be collapsed this hash
                if ((mask[mask_secondary_offset+mask_int] & mask_bit) == 0) {
                    // collapsed will be resetted in any case. Since we merge opposite edges, not likely result will be in this bucket
#ifdef TRACK_COLLAPSED
                    atomic_fetch_add_explicit( &collapse_count[ mask_bit_pos/4 ], (collapse_counter+1) << (8 * (mask_bit_pos % 4)), memory_order_relaxed);
#endif
                    if (passU) {
                        // U-U case, Hash2 has a full info, including the bucket
                        collapse_edges[mask_bit_pos*2 + h1%2] = hash2;
#ifdef TRACK_COLLAPSED
                        const uint hash1mask = (hash1 & UMask) % EDGE_PER_BUCKET;
                        const uint rec_hash1 = UBit | (uint(EDGE_PER_BUCKET) * row_idx + hash1mask);
                        stash_collapsed(collapse_counter, hash2, rec_hash1, collapse_stash);
#endif
                    }
                    else
                    {
                        // V-V case, hash1 has a full info, including the bucket
                        collapse_edges[mask_bit_pos*2 + h1%2] = hash1;
#ifdef TRACK_COLLAPSED
                        const uint hash2mask = hash2 % EDGE_PER_BUCKET;
                        const uint cnt_hash2 = (uint(EDGE_PER_BUCKET) * row_idx + hash2mask);
                        stash_collapsed(collapse_counter, hash1, cnt_hash2, collapse_stash);
#endif
                    }
                    continue;
                }

                // migrating to the other bucket. Any of the line/row will be fine, but let's prefer the same bucket:
                const uint other_bkt_idx = calc_bckt_index(passU, secondary_bucket_idx, tgid);
                device ulong * other_data = get_address<ulong>( params.blocks[other_bkt_idx], buffer0, buffer1);

                int move_idx = atomic_fetch_add_explicit( &bucket_positions[other_bkt_idx], 1, memory_order_relaxed );
                // Legit compare with endIdx, even it is different bucket. Ranges are the same because tid is the same
                if (move_idx<int(params.block_size)) {
                    // need to swap hashed on move, the rpimary one must match the current bucket
#ifdef TRACK_COLLAPSED
                    // Collapse counter needs to migrate from hash1 to hash2.  Hash1 should get back bucket related data.
                    if (u1) {
                        // U-U case
                        const uint hash1mask = (hash1 & UMask) % EDGE_PER_BUCKET;
                        const uint rec_hash1 = UBit | (uint(EDGE_PER_BUCKET) * row_idx + hash1mask);
                        const uint hash2mask = (hash2 & UMask) % EDGE_PER_BUCKET;
                        const uint cnt_hash2 = UBit | (uint(EDGE_PER_BUCKET) * collapse_counter + hash2mask);
                        other_data[move_idx] = (ulong(rec_hash1) << 32) | cnt_hash2;
                    }
                    else {
                        // V-V case
                        const uint hash1mask = hash1 % EDGE_PER_BUCKET;
                        const uint rec_hash1 = (uint(EDGE_PER_BUCKET) * collapse_counter + hash1mask);
                        const uint hash2mask = hash2 % EDGE_PER_BUCKET;
                        const uint cnt_hash2 = (uint(EDGE_PER_BUCKET) * row_idx + hash2mask);

                        other_data[move_idx] = (ulong(rec_hash1) << 32) | cnt_hash2;
                    }
#else
                    other_data[move_idx] = (ulong(hash1) << 32) | hash2;
#endif
                }
#ifdef USE_METRICS
                else {
                    // Note, on else we still want to reset it, because without move - it is likely what will happens in any case
                    // At least we need to report that
                    atomic_fetch_add_explicit(&metric.skipped_edges, 1, memory_order_relaxed);
                }
#endif
            }
            else {
                // U-U or V-V with both edges in current bucket. In case of deletion/collapsing we can't do both, the second one needs to go as stub

                uint mask_bit_pos1 = (h1>>1) % mask_sz_init;
                mask_bit_pos1 =  ((uint)((ulong(mask_bit_pos1)*mask_mul) % EDGE_PER_BUCKET_PREV_PRIME)) % mask_sz_final;

                const uint mask_int1 = (mask_bit_pos1 / 32);
                const uint mask_bit1 = 1 << ((mask_bit_pos1) % 32);

                bool del1 = (mask[mask_int1] & mask_bit1) == 0;
                bool collapse1 = !del1 && ((mask[mask_secondary_offset+mask_int1] & mask_bit1) == 0);

                //processing secondary edge
                uint mask_bit_pos2 = (h2>>1) % mask_sz_init;
                mask_bit_pos2 =  ((uint)((ulong(mask_bit_pos2)*mask_mul) % EDGE_PER_BUCKET_PREV_PRIME)) % mask_sz_final;

                const uint mask_int2 = (mask_bit_pos2 / 32);
                const uint mask_bit2 = 1 << ((mask_bit_pos2) % 32);

                bool del2 = (mask[mask_int2] & mask_bit2) == 0;
                bool collapse2 = !del2 && ((mask[mask_secondary_offset+mask_int2] & mask_bit2) == 0);

                if ((del1 || del2 || collapse1 || collapse2) == false) {
                    // nothing to do, no migration is needed, can stay in the same bucket
                    continue;
                }

                // inplace, resetting to zero, will compact after. Reset happens in ALL branches
                data[idx] = 0;

                // total delete - easy, normal case
                if (del1 && del2) {
                    continue;
                }

                // below we want to use h1 & h1 as full hashes. Let's recover the counter data
#ifdef TRACK_COLLAPSED
                if (collapse1 || collapse2) {
                    // For V-V h1/h2 is already swapped, so H1 is one with counter
                    const uint h1mask = h1 % EDGE_PER_BUCKET;
                    h1 = h1mask + EDGE_PER_BUCKET * row_idx;
                }
#endif

                if (del1 || del2) {
                    if (collapse1 || collapse2) {
                        // collapse needs to be handled by providing stub
                        // collapse_count no need to update, we can use 0 aka -1 because this section is deleted.
                        // Note, h1,h2 are already restored to the bucket data, no collapse counter
                        if (collapse1) {
                            collapse_edges[mask_bit_pos1*2 + h1%2] = ((h1 ^ 0x01) | case_UBit); // writing the same Hash
                        }
                        if (collapse2) {
                            collapse_edges[mask_bit_pos2*2 + h2%2] = ((h2 ^ 0x01) | case_UBit); // writing the same Hash
                        }
                        continue;
                    }
                    else {
                        // just delete without collapse. Normal case
                        continue;
                    }
                }

                if (collapse1 && collapse2) {
                    // can satisfy one but not another
                    // We want randomly satisfy one and drop another one
                    if (h1 % mask_mul > mask_mul/2) {
                        // Satisfy first
                        collapse_edges[mask_bit_pos1*2 + h1%2] = (h2 | case_UBit);
#ifdef TRACK_COLLAPSED
                        stash_collapsed(collapse_counter, h2|case_UBit, h1|case_UBit, collapse_stash);
                        atomic_fetch_add_explicit( &collapse_count[ mask_bit_pos1/4 ], (collapse_counter+1) << (8 * (mask_bit_pos1 % 4)), memory_order_relaxed);
#endif
                        // Put stub into the second
                        collapse_edges[mask_bit_pos2*2 + h2%2] = ((h2 ^ 0x01) | case_UBit);
                        // Zero to the second
                    }
                    else {
                        // Stub in the first
                        collapse_edges[mask_bit_pos1*2 + h1%2] = ( (h1 ^ 0x01) | case_UBit);
                        // Satisfy the second
                        collapse_edges[mask_bit_pos2*2 + h2%2] = (h1 | case_UBit);
#ifdef TRACK_COLLAPSED
                        stash_collapsed(collapse_counter, h1|case_UBit, h2|case_UBit, collapse_stash);
                        atomic_fetch_add_explicit( &collapse_count[ mask_bit_pos2/4 ], (collapse_counter+1) << (8 * (mask_bit_pos2 % 4)), memory_order_relaxed);
#endif
                    }
                }
                else if (collapse1 || collapse2) {  // Must be collapse1 || collapse2
                    if (collapse1) {
                        collapse_edges[mask_bit_pos1*2 + h1%2] = (h2 | case_UBit);
#ifdef TRACK_COLLAPSED
                        stash_collapsed(collapse_counter, h2|case_UBit, h1|case_UBit, collapse_stash);
                        atomic_fetch_add_explicit( &collapse_count[ mask_bit_pos1/4 ], (collapse_counter+1) << (8 * (mask_bit_pos1 % 4)), memory_order_relaxed);
#endif
                    }
                    else {
                        collapse_edges[mask_bit_pos2*2 + h2%2] = (h1 | case_UBit);
#ifdef TRACK_COLLAPSED
                        stash_collapsed(collapse_counter, h1|case_UBit, h2|case_UBit, collapse_stash);
                        atomic_fetch_add_explicit( &collapse_count[ mask_bit_pos2/4 ], (collapse_counter+1) << (8 * (mask_bit_pos2 % 4)), memory_order_relaxed);
#endif
                    }
                }
#ifdef USE_METRICS
                else {
                    // It is ASSERT case, should never happen
                    atomic_fetch_add_explicit(&metric.skipped_edges, 10000, memory_order_relaxed);
                }
#endif
            }
        }
    }
}

kernel void copy_resulting_data(device uint * buffer0[[buffer(0)]],
                        device uint * buffer1[[buffer(1)]],
                        constant TrimmingParams & params[[buffer(2)]],
                        device uint * bucket_positions[[buffer(3)]],
                        device ulong * resulting_data[[buffer(4)]],
                        device uint * resulting_bucket_positions[[buffer(5)]],
                        uint id [[thread_position_in_grid]],
                        uint tid [[thread_index_in_threadgroup]],
                        uint tgid [[threadgroup_position_in_grid]])
{

    // TODO ADD INTO CONFIG
    threadgroup atomic_int read_pos;

    const uint bkt_index = tgid;
    device ulong * data = get_address<ulong>(params.blocks[bkt_index], buffer0, buffer1);

    if (tid==0) {
        atomic_store_explicit( &read_pos, tid, memory_order_relaxed);
    }
    else if (tid==1) {
        resulting_bucket_positions[bkt_index] = min(uint(PHASE1_EDGE_PER_BUCKET), bucket_positions[bkt_index]);
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    while(true) {
        int pos = atomic_fetch_add_explicit(&read_pos, 1, memory_order_relaxed);
        if (pos >= PHASE1_EDGE_PER_BUCKET)
            break;

        resulting_data[bkt_index*PHASE1_EDGE_PER_BUCKET + pos] = data[pos];
    }
}

kernel void discover_nonces(device uint * hh_keys[[buffer(0)]],
                        device uint * nonces[[buffer(1)]],
                        device DiscoverParams & params[[buffer(2)]],
                        uint id [[thread_position_in_grid]])
{
    // TODO put into config
    const uint DISCOVER_NONCES_PER_JOB = 1024;

    const ulong4 keys = {params.key0, params.key1, params.key2, params.key3};
    const uint hash_mask = (uint(1) << EDGE_BITS) - 1;

    uint start_nonce = id*DISCOVER_NONCES_PER_JOB;
    uint end_nonce = start_nonce + DISCOVER_NONCES_PER_JOB;

    const uint nonce_add = params.upass ? 0 : 1;
    const uint hh_len = params.hh_len;
    const uint hh_k = params.hh_k;

    for ( uint nonce=start_nonce; nonce<end_nonce; nonce++ ) {
        uint hash = sipHash24(keys, ulong(nonce)*2+nonce_add) & hash_mask;
        uint hh = hash;
        if ( hh == hh_keys[ (hh*hh_k)%hh_len ] ) {
            // detected
            uint hash1, hash2;
            if (params.upass) {
                hash1 = hash;
                hash2 = sipHash24(keys, ulong(nonce)*2+1) & hash_mask;
            }
            else {
                hash1 = sipHash24(keys, ulong(nonce)*2) & hash_mask;
                hash2 = hash;
            }
            uint nonce_pos = atomic_fetch_add_explicit(&params.nonce_pos, 1, memory_order_relaxed);
            if (nonce_pos<params.nonce_sz) {
                nonces[nonce_pos*3 + 0] = nonce;
                nonces[nonce_pos*3 + 1] = hash1;
                nonces[nonce_pos*3 + 2] = hash2;
            }
        }
    }
}


    )";
}

