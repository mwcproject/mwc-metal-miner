//
// Created by Konstantin Bay on 4/20/25.
//

#include "test_st1_build_buckets.h"
#include <set>
#include <unordered_set>
#include "../miner/sip_hash.h"

#include <iostream>
#include <ostream>

void test_st1_build_buckets(const uint64_t * v,
                        MetalContext & context,
                        uint32_t bucket_allocated_size,
                        const std::vector<MemRange> & st1_8B_buckets,
                        const std::vector<MemRange> & st1_4B_buckets,
                        uint EDGE_BITS)
{
    const uint hash_mask = (uint(1) << EDGE_BITS) - 1;
    const uint BUCKETS_NUM = st1_8B_buckets.size() + st1_4B_buckets.size();
    const uint bucket_size = (uint(1) << EDGE_BITS) / BUCKETS_NUM;

    uint32_t found_nonces = 0;
    uint32_t found_zeroes = 0;

    //uint32_t * buf1 = (uint32_t *)bucket_buffers[0]->contents();

    //uint32_t * skipped_data = (uint32_t *) st1_build_buckets_skipped_nonces->contents();

    //std::vector<uint32_t> found_nonce_vals;

    // Checking 8b data
    for (int b8i = 0; b8i < st1_8B_buckets.size(); b8i++) {
        uint64_t * data = (uint64_t *) st1_8B_buckets[b8i].get_data_ptr(context);
        for (int u=0; u < bucket_allocated_size; u++) {
            uint64_t itm = data[u];
            if (itm==0) {
                found_zeroes++;
                continue;
            }

            found_nonces++;

            uint32_t nonce = (uint32_t)(itm>>32);
            //found_nonce_vals.push_back(nonce);

            uint32_t test_hash = (uint32_t)itm;

            uint hash = sip_hash(v, nonce*2) & hash_mask;
            assert( hash == test_hash );
            assert( hash / bucket_size == b8i );
        }
    }

    // Checking 4b data
    for (int b1i = 0; b1i < st1_4B_buckets.size(); b1i++) {
        uint bucket_idx = b1i + st1_8B_buckets.size();

        uint32_t * data = (uint32_t *) st1_4B_buckets[b1i].get_data_ptr(context);
        for (int u=0; u < bucket_allocated_size; u++) {
            uint64_t itm = data[u];
            if (itm==0xFFFFFFFF) {
                found_zeroes++;
                continue;
            }

            found_nonces++;

            uint32_t nonce = itm;

            //found_nonce_vals.push_back(nonce);
            uint hash = sip_hash(v, nonce*2) & hash_mask;
            uint bidx = hash / bucket_size;
            assert( bidx == bucket_idx );
        }
    }

    //std::sort(found_nonce_vals.begin(), found_nonce_vals.end());
    //std::reverse(found_nonce_vals.begin(), found_nonce_vals.end());

    std::cout << "Found " << found_nonces << " from " << (hash_mask+1) << "  Diff: " << (hash_mask-found_nonces+1) << "  found_zeroes: " << found_zeroes  << std::endl;
    assert( hash_mask-found_nonces+1 < hash_mask / 1000 ); // 0.1% we can lost
}
