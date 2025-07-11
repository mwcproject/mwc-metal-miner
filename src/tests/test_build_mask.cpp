#include "test_build_mask.h"
#include <__filesystem/perms.h>
#include <unordered_set>
#include <iostream>

#include "../miner/sip_hash.h"

#define COMPRESS_BUCKETS_SKIPS 0

CollapsedData::CollapsedData(uint64_t edges, uint32_t edge_per_bucket,
                uint32_t ubkt, uint32_t vbkt) {

    uint32_t h1 = uint32_t(edges);
    uint32_t h2 = uint32_t(edges>>32);

    uint32_t ubit = uint32_t(1) << 31;
    uint32_t umask = ubit-1;

    bool u1 = (h1 & ubit)!=0;
    bool u2 = (h2 & ubit)!=0;

    if (u1!=u2) {
        assert(u1);
        assert(!u2);

#ifdef TRACK_COLLAPSED
        counter = (h1 & umask) / edge_per_bucket;
        uint32_t h1mask = (h1 & umask) % edge_per_bucket;
        hash1 = ubit | ( h1mask + edge_per_bucket * ubkt );
        assert(counter % 2==0);
#else
        counter = 1;
        hash1 = h1;
#endif
        hash2 = h2;
    }
    else {
        if (u1) {
            // U-U
#ifdef TRACK_COLLAPSED
            counter = (h1 & umask) / edge_per_bucket;
            uint32_t h1mask = (h1 & umask) % edge_per_bucket;
            hash1 = ubit | ( h1mask + edge_per_bucket * ubkt );
            assert(counter % 2==1);
#else
            counter = 1;
            hash1 = h1;
#endif
            hash2 = h2;
        }
        else {
            // V-V
            assert(u1==false && u2==false);
#ifdef TRACK_COLLAPSED
            counter = h2 / edge_per_bucket;
            uint32_t h2mask = (h2 & umask) % edge_per_bucket;
            hash2 = h2mask + edge_per_bucket * vbkt;
            assert(counter % 2==1);
#else
            counter = 1;
            hash2 = h2;
#endif
            hash1 = h1;
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////

static uint32_t to_scaled_bit(const uint32_t hash, uint32_t mask_scale_k,
            uint32_t mask_scale_mul,
            uint32_t edges_per_bucket,
            uint32_t mask_next_prime) {
    if (mask_scale_k==1 && mask_scale_mul==1)
        return hash % edges_per_bucket;

    assert( mask_next_prime < edges_per_bucket/2 );
    assert( mask_next_prime > edges_per_bucket/2*0.99 );


    uint32_t updated_mask_sz = edges_per_bucket/2 / mask_scale_k;

    uint32_t h = (hash/2) % (edges_per_bucket/2);

    h = ((uint64_t(h) * mask_scale_mul) % mask_next_prime) % updated_mask_sz;

    h = h*2 + hash%2;

    assert(h>=0);
    assert(h < edges_per_bucket/mask_scale_k );

    assert( h % 2 == hash % 2 );

    return h;
}

static void applyHash2Mask(uint32_t hash, uint32_t mask_scale_k,
            uint32_t mask_scale_mul, uint32_t edges_per_bucket, uint32_t mask_next_prime,
            std::vector<bool> & exp_mask, std::vector<bool> & exp_mask_next) {

    uint32_t idx = to_scaled_bit(hash, mask_scale_k, mask_scale_mul, edges_per_bucket, mask_next_prime);
    if (exp_mask[idx]) {
        exp_mask_next[idx/2] = true;
    }
    else {
        exp_mask[idx] = true;
    }
}

void test_build_buckets_st2(void * mask_data, uint32_t mask_size,
            uint32_t in_bucket_idx,
            const std::vector<uint32_t> & bucket_positions,
            const std::vector< std::vector<uint64_t>> & in,
            const std::vector<std::vector<uint64_t>> & out,
            uint32_t EDGE_BITS, uint32_t buckets_num,
            const uint64_t v[4],
            const std::unordered_map<uint64_t, CollapsedData> & found_collapsed_data)
{

    uint32_t edges_per_bucket = ( uint32_t(1) << EDGE_BITS ) / buckets_num;

    const uint hash_mask = (uint(1) << EDGE_BITS) - 1;
    const uint edge_per_bucket = (hash_mask+1)/buckets_num;

    int in_sz = in.size();

    std::vector<bool> exp_mask(edges_per_bucket, false);
    std::vector<bool> exp_mask_next(edges_per_bucket/2, false);

    //int total_in = 0;
    for (int i = 0; i < in_sz; i++) {
        const std::vector<uint64_t> & data = in[i];
        for (int j = 0; j < data.size(); j++) {
            uint64_t dt = data[j];
            assert(dt!=0x6363636363636363ULL);
            if (!dt)
                continue;

            //total_in++;
            uint32_t hash;

            // Expected U-V, U-U or V-V
            hash = uint32_t(dt);
            // nonce to hash
            uint32_t nonce = uint32_t(dt>>32);
            uint32_t exp_hash = sip_hash(v, uint64_t(nonce)*2) & hash_mask;
            assert(exp_hash == hash);

            int hash_bucket_idx = hash / edges_per_bucket;
            assert(hash_bucket_idx == in_bucket_idx);

            applyHash2Mask(hash,  1,
                        1, edges_per_bucket, 1,
                        exp_mask, exp_mask_next);
        }
    }

    std::vector<bool> no_trim_copy(exp_mask);

    // Let's trim the mask
    int edges2plus = 0;
    int edges2 = 0;
    for (int u=0; u<edges_per_bucket; u+=2) {
        if (! (exp_mask[u] && exp_mask[u+1])) {
            exp_mask[u] = false;
            exp_mask[u+1] = false;
        }
        else {
            if (exp_mask_next[u/2]) {
                edges2plus++;
            }
            else {
                edges2++; // can be trimmed
            }
        }
    }

    std::cout << "Edges 2+: " << edges2plus << "  Can be trimmed as twos: " << edges2 << std::endl;

    // Now we can compare what we got to test...
    std::vector<uint32_t> mask_copy(mask_size/4);
    memcpy( mask_copy.data(), mask_data, mask_size );
    uint mask_len = mask_size/4;
    assert(mask_len % 3 == 0);
    assert(mask_len >= edges_per_bucket/32 + edges_per_bucket/2/32);
    uint32_t mask_chunk_sz = edges_per_bucket/32/2;
    assert(mask_chunk_sz*3 <= mask_len);
    assert(mask_chunk_sz*3 <= mask_copy.size());
    //uint32_t * mask_odd_data = mask_copy.data()+mask_chunk_sz;
    uint32_t * mask_secondary_data = mask_copy.data()+mask_chunk_sz*2;
    // Checking trimmed mask
    for (int u=0; u<mask_chunk_sz; u++) {
        assert(u<mask_len/2);
        uint32_t even = mask_copy[u];
        for (int v=0; v<32; v++) {
            uint32_t bit = uint32_t(1) << v;
            int idx1 = u*64 + v*2;
            bool exp1 = (even & bit) != 0;

            assert(exp_mask[idx1] == exp1);
            assert(exp_mask[idx1+1] == exp1);
        }
    }
    // Checking secondary mask
    for (int u=0; u<mask_chunk_sz; u++) {
        uint32_t secondaryMask = mask_secondary_data[u];
        for (int v=0; v<32; v++) {
            uint32_t bit = uint32_t(1) << v;
            int idx = u*32 + v;
            bool exp_sec = (secondaryMask & bit) != 0;

            assert(exp_mask_next[idx] == exp_sec);
        }
    }

    // Checking resulting buckets
    std::vector< std::unordered_set<uint64_t>> exp_buckets(buckets_num);

    std::unordered_map<uint32_t, uint32_t> collapsed_data;
    int found_collapsed_edges = 0;
    for (int bi = 0; bi < in_sz; bi++) {
        const std::vector<uint64_t> & data = in[bi];

        for (int j = 0; j < data.size(); j++) {
            uint64_t dt = data[j];
            if (!dt)
                continue;

            uint32_t hash = uint32_t(dt);
            int hash_bucket_idx = hash / edges_per_bucket;
            assert(hash_bucket_idx == in_bucket_idx);

            uint32_t idx = to_scaled_bit(hash, 1, 1, edges_per_bucket, 1);

            if (!exp_mask[idx])
                continue;

            uint32_t nonce = uint32_t(dt >> 32);
            uint32_t hash2 = sip_hash(v, uint64_t(nonce)*2+1) & hash_mask;
            uint bucket_idx = hash2 / edge_per_bucket;

            if (!exp_mask_next[idx/2]) { // it is collapsed data, processing it separately

                //
                if (collapsed_data.count(idx/2) == 0) {
                    collapsed_data[idx/2] = hash2;
                    continue;
                }
                else {
                    uint32_t other_hash = collapsed_data[idx/2];
                    collapsed_data.erase(idx/2);

                    // At st2, other_hash and hash2 are both Vs
                    if ( other_hash > hash2 )
                        std::swap(other_hash, hash2);

                    assert(other_hash < hash2);

                    CollapsedData col_dt(hash2, other_hash, 1);

                    assert( found_collapsed_data.count(col_dt.toKey()) > 0 );
                    found_collapsed_edges++;
                }
            }
            else {
                // U-V  standard pair
                // Hash must have counter 0
#ifdef TRACK_COLLAPSED
                uint64_t dt = (hash % edge_per_bucket) | (uint32_t(1) << 31) | (uint64_t(hash2) << 32);
#else
                uint64_t dt = hash | (uint32_t(1) << 31) | (uint64_t(hash2) << 32);
#endif
                auto r = exp_buckets[bucket_idx].insert(dt);
                assert(r.second);
           }
        }
    }

    assert(collapsed_data.empty());
    assert(found_collapsed_edges == found_collapsed_data.size());

    uint32_t out_buckets_sz = out.size();
    assert(in.size() == out.size());

    for (int bi = 0; bi < out_buckets_sz; bi++) {
        const std::vector<uint64_t> & data = out[bi];

        uint num = 0;
        uint notfound = 0;
        for (int j = 0; j < data.size(); j++) {
            const uint64_t dt = data[j];
            if (!dt)
                continue;

            if ( exp_buckets[bi].count(dt) > 0 ) {
                num++;
            }
            else {
                // collapsed are not expected at this point, they are not processed yet
                notfound++;
            }
        }

        uint32_t buket_pos = bucket_positions[ in_bucket_idx*buckets_num + bi ];

        assert(buket_pos >= num); // metrics includes some zeroes
        assert(buket_pos <= num + COMPRESS_BUCKETS_SKIPS); // overhead is not expected to be high
        assert(exp_buckets[bi].size() == num);
        assert(notfound < 10);
        // exact_match_expected. If failed, enable metrics and check that there is no data loss
        assert(num <= exp_buckets[bi].size());
        assert(num >= exp_buckets[bi].size()-exp_buckets[bi].size()/100); // duplicates are possible
    }
}

static void process_collapse(std::unordered_map<uint32_t, uint32_t> & collapsed_data,
                std::unordered_map<uint32_t, int> & collapsed_counters,
                uint32_t CYCLE_LEN,
                uint32_t hash_meet_idx, uint32_t COLLAPSE_STASH_COUNTER,
                int32_t idx,
                uint32_t hash,
#ifdef TRACK_COLLAPSED
                int counter,
#endif
                const std::unordered_map<uint64_t, CollapsedData> & in_collapsed_data,
                std::unordered_map<uint64_t, CollapsedData>  &collapsed_res,
                uint32_t collapse_stash_prev_size,
                uint32_t * collapse_stash)
{
    assert(collapsed_data.size() == collapsed_counters.size());

#ifdef TRACK_COLLAPSED
    if (counter >= int(COLLAPSE_STASH_COUNTER)) {
        uint32_t coll_sz = collapse_stash[0];
        assert(coll_sz >= collapse_stash_prev_size);

        uint32_t hash_hash = (hash & 0x07FFFFFF) | ( uint32_t(counter - COLLAPSE_STASH_COUNTER) << (32-5) );
        bool found = false;
        for (uint32_t j=collapse_stash_prev_size; j<coll_sz; j++) {
            if (hash_hash == collapse_stash[1 + j*2] ) {
                if (hash_meet_idx == collapse_stash[1 + j*2+1]) {
                    found = true;
                }
            }
        }
        assert(found);
    }
#endif

    if (collapsed_data.count(idx/2) == 0) {
        collapsed_data[idx/2] = hash;
#ifdef TRACK_COLLAPSED
        collapsed_counters[idx/2] = counter;
#else
        collapsed_counters[idx/2] = 0;
#endif
        assert(collapsed_data.size() == collapsed_counters.size());
    }
    else {
        uint32_t other_hash = collapsed_data[idx/2];
        int other_counter = collapsed_counters[idx/2];
        collapsed_data.erase(idx/2);
        collapsed_counters.erase(idx/2);
        assert(collapsed_data.size() == collapsed_counters.size());

#ifdef TRACK_COLLAPSED
        int res_counter = other_counter+counter+1;
        if (res_counter<0 || res_counter>CYCLE_LEN-1)
            return;
        if (other_hash/2 == hash/2) {
            // assert(other_hash != hash);
            if (res_counter!=CYCLE_LEN-1)
                return;
        }
        assert(res_counter>=0);
        CollapsedData col_dt(other_hash, hash, res_counter);
#else
        CollapsedData col_dt(other_hash, hash, 1);
#endif

        // collapsed_res value, any order of hashes is fine
        assert(collapsed_res.count(col_dt.toKey()) == 0);
        collapsed_res[col_dt.toKey()] = col_dt;

/*        if (in_collapsed_data.count(col_dt.toKey()) == 0) {
            for (const auto & in_c : in_collapsed_data ) {
                if (in_c.second.hash1 == other_hash || in_c.second.hash1 == hash || in_c.second.hash2 == other_hash || in_c.second.hash2 == hash) {
                    int rr=0;
                }
            }
        }*/

        assert( in_collapsed_data.count(col_dt.toKey()) > 0 );
    }
}

void test_build_mask_buckets(void * mask_data, uint32_t mask_size,
            uint32_t in_bucket_idx,
            const std::vector<uint32_t> & bucket_positions,
            const std::vector< std::vector<uint64_t>> & in,
            uint32_t mask_scale_k,
            uint32_t mask_scale_mul,
            const std::vector<std::vector<uint64_t>> & out,
            uint32_t EDGE_BITS, uint32_t buckets_num, uint32_t mask_next_prime,
            bool passU,
            uint32_t CYCLE_LEN,
            const std::unordered_map<uint64_t, CollapsedData> & in_collapsed_data,
            const std::unordered_set<uint64_t> & in_migrated_data,
            std::unordered_map<uint64_t, CollapsedData> & collapsed_res, // expected collapsed data
            uint32_t COLLAPSE_STASH_COUNTER,
            uint32_t collapse_stash_prev_size,
            uint32_t * collapse_stash) {
    //uint32_t * buf1 = (uint32_t *) bucket_buffers[0]->contents();
//    uint32_t * buf2 = (uint32_t *) bucket_buffers[1]->contents();

    assert(collapsed_res.empty());
    uint32_t edges_per_bucket = ( uint32_t(1) << EDGE_BITS ) / buckets_num;

    if (!(mask_scale_k==1 && mask_scale_mul==1)) {
        assert( mask_next_prime < edges_per_bucket/2 );
        assert( mask_next_prime > edges_per_bucket/2*0.99 );
    }

    const uint hash_mask = (uint(1) << EDGE_BITS) - 1;
    const uint32_t UBit = uint32_t(1) << 31;
    const uint32_t UMask = UBit - 1;
    const uint edge_per_bucket = (hash_mask+1)/buckets_num;

    int in_sz = in.size();

    std::vector<bool> exp_mask(edges_per_bucket, false);
    std::vector<bool> exp_mask_next(edges_per_bucket/2, false);

    int total_in = 0;
    for (int i = 0; i < in_sz; i++) {
        const std::vector<uint64_t> & data = in[i];
        for (int j = 0; j < data.size(); j++) {
            uint64_t dt = data[j];
            assert(dt!=0x6363636363636363ULL);
            if (!dt)
                continue;

            total_in++;
            uint32_t hash;

            // Expected U-V, U-U or V-V
            if (passU) {
                hash = uint32_t(dt);

                uint32_t hash2 = uint32_t(dt>>32);
                if (hash & UBit)
                    hash ^= UBit;
                else {
                    assert((hash2 & UBit) == 0);
                    continue; // V-V case
                }
                // U-V or U-U case
            }
            else {
                uint32_t hash2 = uint32_t(dt);
                hash = uint32_t(dt>>32);
                if (hash & UBit) {
                    assert(hash2 & UBit);
                    continue; // U-U case
                }
                // U-V or V-V case
                // in case of V-V
            }

#ifndef TRACK_COLLAPSED
            int hash_bucket_idx = hash / edges_per_bucket;
            assert(hash_bucket_idx == in_bucket_idx);
#endif

            applyHash2Mask(hash,  mask_scale_k,
                        mask_scale_mul, edges_per_bucket, mask_next_prime,
                        exp_mask, exp_mask_next);

            if (passU) {
                hash = uint32_t(dt>>32);
                if (hash & UBit)
                    hash ^= UBit; // U-U
                else
                    continue; // U-V case
            }
            else {
                hash = uint32_t(dt);
                if (hash & UBit) {
                    continue; // U-V case
                }
                // V-V
            }

            int hash_bucket_idx2 = hash / edges_per_bucket;
            if (hash_bucket_idx2!=in_bucket_idx)
                continue;


            applyHash2Mask( hash, mask_scale_k,
                        mask_scale_mul, edges_per_bucket, mask_next_prime,
                        exp_mask, exp_mask_next);
        }
    }

    (void)total_in;
    std::vector<bool> no_trim_copy(exp_mask);

    // Let's trim the mask
    int edges2plus = 0;
    int edges2 = 0;
    for (int u=0; u<edges_per_bucket; u+=2) {
        if (! (exp_mask[u] && exp_mask[u+1])) {
            exp_mask[u] = false;
            exp_mask[u+1] = false;
        }
        else {
            if (exp_mask_next[u/2]) {
                edges2plus++;
            }
            else {
                edges2++; // can be trimmed
            }
        }
    }

    std::cout << "Edges 2+: " << edges2plus << "  Can be trimmed as twos: " << edges2 << std::endl;

    // Now we can compare what we got to test...
    //uint32_t * mask_data = (uint32_t *) mask->contents();
    std::vector<uint32_t> mask_copy(mask_size/4);
    memcpy( mask_copy.data(), mask_data, mask_size );
    uint mask_len = mask_size/4;
    assert(mask_len % 3 == 0);
    assert(mask_len >= edges_per_bucket/32 + edges_per_bucket/2/32);
    uint32_t mask_chunk_sz = edges_per_bucket/32/2;
    assert(mask_chunk_sz*3 <= mask_len);
    assert(mask_chunk_sz*3 <= mask_copy.size());
    //uint32_t * mask_odd_data = mask_copy.data()+mask_chunk_sz;
    uint32_t * mask_secondary_data = mask_copy.data()+mask_chunk_sz*2;
    // Checking trimmed mask
    for (int u=0; u<mask_chunk_sz; u++) {
        assert(u<mask_len/2);
        uint32_t even = mask_copy[u];
        //uint32_t odd = mask_data2[u];
        for (int v=0; v<32; v++) {
            uint32_t bit = uint32_t(1) << v;
            int idx1 = u*64 + v*2;
            //int idx2 = u*64 + v*2 + 1;
            bool exp1 = (even & bit) != 0;
            //bool exp2 = (odd & bit) != 0;

            assert(exp_mask[idx1] == exp1);
            assert(exp_mask[idx1+1] == exp1);
            //assert(exp_mask[idx2] == exp2);
        }
    }
    // Checking secondary mask
    for (int u=0; u<mask_chunk_sz; u++) {
        uint32_t secondaryMask = mask_secondary_data[u];
        for (int v=0; v<32; v++) {
            uint32_t bit = uint32_t(1) << v;
            int idx = u*32 + v;
            bool exp_sec = (secondaryMask & bit) != 0;

            assert(exp_mask_next[idx] == exp_sec);
        }
    }


    // Checking resulting buckets
    std::vector< std::unordered_set<uint64_t>> exp_buckets(buckets_num);

    std::unordered_map<uint32_t, uint32_t> collapsed_data;
    std::unordered_map<uint32_t, int> collapsed_counters;

    int total_by_bucket_pos = 0;
    if (passU) {
        for (int j=0;j<buckets_num;j++)
            total_by_bucket_pos+=bucket_positions[in_bucket_idx*buckets_num + j];
    }
    else {
        for (int j=0;j<buckets_num;j++)
            total_by_bucket_pos+=bucket_positions[in_bucket_idx*j + buckets_num];
    }

    (void)total_by_bucket_pos;

    //std::unordered_set<uint64_t> skipped_edges;
    int migrated_num = 0;
    int trimmed = 0;
    int collapsed = 0;
    int total = 0;
    for (int bi = 0; bi < in_sz; bi++) {
        const std::vector<uint64_t> & data = in[bi];

        for (int j = 0; j < data.size(); j++) {
            const uint64_t dt = data[j];
            if (!dt)
                continue;

            total++;

                // Normal workflow
                uint32_t hash1 = uint32_t(dt);
                uint32_t hash2 = uint32_t(dt>>32);

                bool u1 = (hash1 & UBit)!=0;
                bool u2 = (hash2 & UBit)!=0;

                if (u1 != u2) {
                    assert(u1); // U-V type is expected
                    assert(!u2);

                    uint32_t hash = passU ? (hash1 ^ UBit) : hash2;

#ifdef TRACK_COLLAPSED
                    int counter =  (hash1 & UMask) / edges_per_bucket;

#else
                    int hash_bucket_idx = hash / edges_per_bucket;
                    assert(hash_bucket_idx == in_bucket_idx);
#endif

                    uint32_t idx = to_scaled_bit(hash, mask_scale_k, mask_scale_mul, edges_per_bucket, mask_next_prime);

                    if (!exp_mask[idx]) {
                        trimmed++;
                        continue;
                    }

                    if (!exp_mask_next[idx/2]) {
                        uint32_t hashSecond = passU ? hash2 : hash1;
                        collapsed++;
                        if (!passU) {
                            // hashSecond aka hash1 needs to restore buket data
                            hashSecond = ((hash1 & UMask) % edges_per_bucket + bi * edges_per_bucket) | UBit;
                        }
                        else {
                            hash = ((hash & UMask) % edges_per_bucket + bi * edges_per_bucket) | UBit;
                        }

                        process_collapse(collapsed_data, collapsed_counters, CYCLE_LEN, hash, COLLAPSE_STASH_COUNTER, idx,hashSecond,
#ifdef TRACK_COLLAPSED
                                            counter,
#endif
                                            in_collapsed_data,collapsed_res, collapse_stash_prev_size, collapse_stash);
                        continue;
                    }
                    else {
                        auto r = exp_buckets[bi].insert(dt);
                        assert(r.second);
                    }
                }
                else {
                    assert(u1==u2);
                    if (u1==passU) {
                        if (!passU) { // for U-U first hash is a primary. for V-V the second hash is primary
                            std::swap(hash1,hash2);
                        }
                        // At this point first hash is primary (for current block), second is secondary

                        uint32_t h1 = hash1 & 0x7FFFFFFF;
                        uint32_t h2 = hash2 & 0x7FFFFFFF;

#ifdef TRACK_COLLAPSED
                        int counter = h1 / edges_per_bucket;
#else
                        int hash1_bucket_idx = h1 / edges_per_bucket;
                        assert(hash1_bucket_idx == in_bucket_idx);
#endif
                        int hash2_bucket_idx = h2 / edges_per_bucket;

                        if (hash2_bucket_idx != in_bucket_idx) {
                            // general case
                            uint32_t idx = to_scaled_bit(h1, mask_scale_k, mask_scale_mul, edges_per_bucket, mask_next_prime);

                            if (!exp_mask[idx]) {
                                trimmed++;
                                continue;
                            }

                            if (!exp_mask_next[idx/2]) {
                                hash1 = ((hash1 & UMask) % edges_per_bucket + in_bucket_idx * edges_per_bucket);
                                if (u1)
                                    hash1 |= UBit;

                                process_collapse(collapsed_data, collapsed_counters,CYCLE_LEN,hash1,COLLAPSE_STASH_COUNTER, idx,hash2,
#ifdef TRACK_COLLAPSED
                                            counter,
#endif
                                            in_collapsed_data,collapsed_res, collapse_stash_prev_size, collapse_stash );
                                collapsed++;
                                continue;
                            }

                            // Evething is fine, keeping the edges, but need move to another bucket
                            // Migraiton reqires hash swap
                            uint32_t hash1 = uint32_t(dt);
                            uint32_t hash2 = uint32_t(dt>>32);
                            assert(u1==u2);
#ifdef TRACK_COLLAPSED
                            if (u1) {
                                // U-U
                                uint32_t hash1_mask = (hash1 & UMask) % edge_per_bucket;
                                hash1 = UBit | (hash1_mask + in_bucket_idx * edge_per_bucket);

                                uint32_t hash2_mask = (hash2 & UMask) % edge_per_bucket;
                                hash2 = UBit | (hash2_mask + counter * edge_per_bucket);
                            }
                            else {
                                // V-V
                                uint32_t hash1_mask = hash1 % edge_per_bucket;
                                hash1 = hash1_mask + counter * edge_per_bucket;

                                uint32_t hash2_mask = hash2 % edge_per_bucket;
                                hash2 = hash2_mask + in_bucket_idx * edge_per_bucket;
                            }
#endif
                            uint64_t inv_dt = uint64_t(hash2) | (uint64_t(hash1)<<32);
                            assert(in_migrated_data.count(inv_dt)>0);
                            migrated_num++;
                            continue;
                        }
                        else {
                            assert(hash2_bucket_idx == in_bucket_idx);

                            // both edges in the same THIS bucket
                            uint32_t idx1 = to_scaled_bit(h1, mask_scale_k, mask_scale_mul, edges_per_bucket, mask_next_prime);
                            uint32_t idx2 = to_scaled_bit(h2, mask_scale_k, mask_scale_mul, edges_per_bucket, mask_next_prime);

                            bool del1 = !exp_mask[idx1];
                            bool del2 = !exp_mask[idx2];
                            if (del1 && del2) {
                                trimmed++;
                                continue;
                            }

                            bool collapsed1 = false;
                            if (!del1)
                                collapsed1 = !exp_mask_next[idx1/2];

                            bool collapsed2 = false;
                            if (!del2)
                                collapsed2 = !exp_mask_next[idx2/2];

                            // restoring hash1, it can be U or V. Bucket will be in_bucket_idx in both cases.
                            uint32_t hash1_mask = (hash1 & UMask) % edge_per_bucket;
                            hash1 = (hash1 & UBit) | (hash1_mask + in_bucket_idx * edge_per_bucket);

                            if (del1 || del2) {
                                if (collapsed1) {
                                    // processing as a stab
                                    process_collapse(collapsed_data, collapsed_counters, CYCLE_LEN, hash1, COLLAPSE_STASH_COUNTER, idx1, hash1 ^ 0x01,
#ifdef TRACK_COLLAPSED
                                                -1,
#endif
                                                in_collapsed_data,collapsed_res, collapse_stash_prev_size, collapse_stash );
                                    continue;
                                }
                                else if (collapsed2) {
                                    // processing as a stab
                                    process_collapse(collapsed_data, collapsed_counters, CYCLE_LEN, hash2, COLLAPSE_STASH_COUNTER, idx2,hash2 ^ 0x01,
#ifdef TRACK_COLLAPSED
                                                -1,
#endif
                                                in_collapsed_data,collapsed_res, collapse_stash_prev_size, collapse_stash );
                                    continue;
                                }
                                else {
                                    // just a single del
                                    continue;
                                }
                            }
                            else {
                                assert(del1==false);
                                assert(del2==false);

                                if (collapsed1 && collapsed2) {
                                    if ((hash1 & UMask) % mask_scale_mul > mask_scale_mul/2) {
                                        // first collapse is normal, the second one is a failure
                                        process_collapse(collapsed_data,collapsed_counters, CYCLE_LEN, hash1, COLLAPSE_STASH_COUNTER, idx1,hash2,
#ifdef TRACK_COLLAPSED
                                                counter,
#endif
                                                in_collapsed_data,collapsed_res, collapse_stash_prev_size, collapse_stash);
                                        // Second collapse with stab
                                        process_collapse(collapsed_data,collapsed_counters, CYCLE_LEN, hash2, COLLAPSE_STASH_COUNTER, idx2,hash2 ^ 0x01,
#ifdef TRACK_COLLAPSED
                                                -1,
#endif
                                                in_collapsed_data,collapsed_res, collapse_stash_prev_size, collapse_stash);
                                    }
                                    else {
                                        // First with a stub
                                        process_collapse(collapsed_data,collapsed_counters, CYCLE_LEN, hash1, COLLAPSE_STASH_COUNTER, idx1,hash1 ^ 0x01,
#ifdef TRACK_COLLAPSED
                                                -1,
#endif
                                                in_collapsed_data,collapsed_res, collapse_stash_prev_size, collapse_stash);
                                        // Second with a data
                                        process_collapse(collapsed_data,collapsed_counters, CYCLE_LEN, hash2, COLLAPSE_STASH_COUNTER, idx2,hash1,
#ifdef TRACK_COLLAPSED
                                                counter,
#endif
                                                in_collapsed_data,collapsed_res, collapse_stash_prev_size, collapse_stash);
                                    }
                                    continue;
                                }
                                else if (collapsed1) {
                                    assert(!collapsed2);
                                    process_collapse(collapsed_data,collapsed_counters,CYCLE_LEN, hash1, COLLAPSE_STASH_COUNTER, idx1,hash2,
#ifdef TRACK_COLLAPSED
                                                counter,
#endif
                                                in_collapsed_data,collapsed_res, collapse_stash_prev_size, collapse_stash);
                                    continue;
                                }
                                else if (collapsed2) {
                                    assert(!collapsed1);
                                    process_collapse(collapsed_data,collapsed_counters,CYCLE_LEN, hash2, COLLAPSE_STASH_COUNTER, idx2,hash1,
#ifdef TRACK_COLLAPSED
                                                counter,
#endif
                                                in_collapsed_data,collapsed_res, collapse_stash_prev_size, collapse_stash);
                                    continue;
                                }
                                else {
                                    assert(!collapsed1);
                                    assert(!collapsed2);
                                    // nothing happens, keeping this data into the bucket
                                    auto r = exp_buckets[bi].insert(dt);
                                    assert(r.second);
                                    continue;
                                }
                            }
                        }

                        assert(false);
                    }
                    else {
                        // supposed to be skipped, keeping into the same bucket
                        auto r = exp_buckets[bi].insert(dt);
                        assert(r.second);
                    }
                }
        }
    }

    (void)trimmed;
    (void)collapsed;
    (void)total;

    // Checkiong what exactly is a problem in case of mismatch
    if (collapsed_res.size() != in_collapsed_data.size()) {
        std::unordered_map<uint64_t, CollapsedData> copy11(collapsed_res);
        std::unordered_map<uint64_t, CollapsedData> copy12(in_collapsed_data);
        for (const auto i1 : collapsed_res) {
            copy11.erase(i1.first);
            copy12.erase(i1.first);
        }
        std::unordered_map<uint64_t, CollapsedData> copy21(collapsed_res);
        std::unordered_map<uint64_t, CollapsedData> copy22(in_collapsed_data);
        for (const auto i1 : in_collapsed_data) {
            copy21.erase(i1.first);
            copy22.erase(i1.first);
        }
        assert(false);
    }

    assert(collapsed_res.size() == in_collapsed_data.size());

    // Let's compare the values
    for (const auto & c : collapsed_res) {
        assert(in_collapsed_data.count(c.first) == 1 );
        assert( c.second == in_collapsed_data.at(c.first) );
    }

    assert(migrated_num == in_migrated_data.size());
    assert(collapsed_data.empty());
    assert(collapsed_counters.empty());
    // smaller because both edges can be non pair, it is mean nothing in collapsed_data

    uint32_t out_buckets_sz = out.size();
    assert(in.size() == out.size());

    for (int bi = 0; bi < out_buckets_sz; bi++) {
        const std::vector<uint64_t> & data = out[bi];

        uint num = 0;
        uint notfound = 0;
        for (int j = 0; j < data.size(); j++) {
            const uint64_t dt = data[j];
            if (!dt)
                continue;

            if ( exp_buckets[bi].count(dt) > 0 ) {
                num++;
            }
            else {
                assert(false);
            }
        }

        uint bucket_positions_idx = passU ? in_bucket_idx * in.size() + bi : bi * in.size() + in_bucket_idx;
        uint32_t bucket_pos = bucket_positions[bucket_positions_idx];

        assert(bucket_pos >= num);
        assert(bucket_pos <= num + COMPRESS_BUCKETS_SKIPS);
        assert(exp_buckets[bi].size() == num);
        assert(notfound < 10);

        assert(num <= exp_buckets[bi].size());
        assert(exp_buckets[bi].size() - num < std::max(uint(3), uint(exp_buckets[bi].size()) / 100) );
    }
}

uint32_t validate_bucket_data( bool hasUVbit, uint32_t EDGE_BITS, uint32_t buckets_num,
                    int bucketU, int bucketV, uint32_t bucket_position,
                    const std::vector<uint64_t> & data )
{
    uint32_t edges_per_bucket = ( uint32_t(1) << EDGE_BITS ) / buckets_num;
    const uint32_t UBit = uint32_t(1) << 31;

    uint32_t found_records = 0;

    uint32_t zeroes_before_bucket_position = 0;

    for (int i=0;i<data.size(); i++) {
        const uint64_t dt = data[i];
        assert(dt != 0x6363636363636363ULL);
        if (!dt) {
            if (hasUVbit) {
                if (i<bucket_position) {
                    //const uint64_t * dddt = &data[i];
                    zeroes_before_bucket_position++;
                }
            }
            continue;
        }

        if (hasUVbit) {
            assert(i<bucket_position);
        }

        found_records++;

        uint32_t hash1 = uint32_t(dt);
        uint32_t hash2 = uint32_t(dt>>32);

        bool u1 = (hash1 & UBit) != 0;
        bool u2 = (hash2 & UBit) != 0;

        if (!hasUVbit) {
            assert(bucketV<0);
            int bkt1 = hash1 / edges_per_bucket;
            assert(bkt1 == bucketU);
            continue;
        }

        // General case, both buckets are defined, UV, UU, VV combinations are expected.
        // For UU first should match the bucket.
        // For VV the second must match the bucket.

        if (u1!=u2) {
            // U-V case
            assert(u1);
            assert(!u2);
#ifndef TRACK_COLLAPSED
            int bkt1 = (hash1 & 0x7FFFFFFF) / edges_per_bucket;
            assert(bkt1 == bucketU);
#else
            uint32_t counter = (hash1 & 0x7FFFFFFF) / edges_per_bucket;
            assert(counter % 2 == 0); // it is about collapsed removed, 0 to init
#endif
            int bkt2 = hash2 / edges_per_bucket;
            assert(bkt2 == bucketV);
        }
        else {
            assert(u1==u2);
            if (u1) {
#ifndef TRACK_COLLAPSED
                int bkt1 = (hash1 & 0x7FFFFFFF) / edges_per_bucket;
                assert(bkt1 == bucketU);
#else
                uint32_t counter = (hash1 & 0x7FFFFFFF) / edges_per_bucket;
                assert(counter%2 ==1 );
#endif
                // another bucket can be any, we can't check
            }
            else {
#ifndef TRACK_COLLAPSED
                int bkt2 = hash2 / edges_per_bucket;
                assert(bkt2 == bucketV);
#else
                uint32_t counter = hash2 / edges_per_bucket;
                assert(counter%2 ==1 );
#endif
                // another bucket can be any, we can't check
            }
        }
    }

    assert(zeroes_before_bucket_position<=COMPRESS_BUCKETS_SKIPS); // the real limit is unknown, but we don't want this number to be large
    return found_records;
}

// Simply check if all the data is fill with zeros
void validate_empty_mask_buffers( void * mask_buffer, uint32_t mask_buffer_size ) {
    assert(mask_buffer_size%4==0);
    uint32_t * data = (uint32_t *) mask_buffer;
    for (uint32_t i=0; i<mask_buffer_size/4; i++) {
        assert( data[i] == 0 );
    }
}


void validate_empty_collapse_buffer( void* collapse_buffer, uint32_t collapse_buffer_size,
#ifdef TRACK_COLLAPSED
                void* collapse_edges_counts, uint32_t collapse_edges_counts_size,
#endif
                    uint32_t edge_per_bucket )
{
    assert(collapse_buffer_size == edge_per_bucket*4);
    uint32_t * data = (uint32_t *) collapse_buffer;
    for (int bi = 0; bi < edge_per_bucket; bi++) {
        assert(data[bi] == 0xFFFFFFFF);
    }
#ifdef TRACK_COLLAPSED
    assert(collapse_buffer_size == collapse_edges_counts_size*8);
    uint32_t * data_count = (uint32_t *) collapse_edges_counts;
    for (int bi = 0; bi < edge_per_bucket/8; bi++) {
        assert(data_count[bi] == 0);
    }
#endif
}


void validate_filled_collapse_buffer(
            void * collapse_buffer_data_ptr, uint32_t collapse_buffer_size,
#ifdef TRACK_COLLAPSED
            void * collapse_edges_counts_data_ptr, uint32_t collapse_edges_counts_size,
#endif
            uint32_t edge_per_bucket,
            uint32_t EDGE_BITS, uint32_t CYCLE_LEN, uint32_t buckets_num,
            void* mask_data, uint32_t mask_size,
             std::unordered_map<uint64_t, CollapsedData> & resulting_collapsed_data ) {

    assert(resulting_collapsed_data.empty());

    uint32_t edges_per_bucket = (uint32_t(1) << EDGE_BITS) / buckets_num;

    // size as uint for half mask
    uint32_t half_mask_size = edges_per_bucket / 2 / 32;
    assert(half_mask_size*3*4 == mask_size );
    uint32_t * primary_mask_data = (uint32_t *) mask_data;
    uint32_t * secondary_mask_data = primary_mask_data + half_mask_size*2;

    assert(edge_per_bucket == (uint32_t(1) << EDGE_BITS) / buckets_num);
    assert(edge_per_bucket == half_mask_size * 32 * 2);
    assert(collapse_buffer_size == edge_per_bucket*4);
    std::vector<uint32_t> collapse_data(collapse_buffer_size/4, 0);
    memcpy(collapse_data.data(), collapse_buffer_data_ptr, collapse_buffer_size);

#ifdef TRACK_COLLAPSED
    assert(collapse_edges_counts_size == edge_per_bucket/2);
    std::vector<uint32_t> collapse_data_counts(collapse_edges_counts_size/4, 0);
    memcpy(collapse_data_counts.data(), collapse_edges_counts_data_ptr, collapse_edges_counts_size);
#endif

    //const uint32_t UBit = uint32_t(1) << 31;

    for (int mi = 0; mi < half_mask_size; mi++) {
        uint32_t prm = primary_mask_data[mi];
        uint32_t sec = secondary_mask_data[mi];
        uint32_t secN = ~sec;
        uint32_t collapsed = (prm & secN);
        for (uint b=0; b < 32; b++) {
            uint32_t bit_mask = uint32_t(1) << b;
            uint32_t collapse_index = (mi*32 + b) * 2;
            if (collapsed & bit_mask) {
                uint32_t hash1 = collapse_data[collapse_index];
                uint32_t hash2 = collapse_data[collapse_index+1];

                assert(hash1!=0xFFFFFFFF);
                assert(hash2!=0xFFFFFFFF);

                int counter = 0;
#ifdef TRACK_COLLAPSED
                counter = collapse_data_counts[collapse_index/8];
                counter >>= 8 * ((collapse_index/2) % 4);
                counter &= 0xFF;
                counter = counter+1-2;
                if (counter<0 || counter>=CYCLE_LEN) // Doesn't make sense to keep longer paths.
                    continue; // should be skipped

                if (hash1/2 == hash2/2) {
                    //assert(hash1 != hash2);
                    if (counter!=CYCLE_LEN-1) // 41 is a solution
                        continue; //
                }

#endif

                CollapsedData col_dt(hash1, hash2, counter);

                assert(resulting_collapsed_data.count(col_dt.toKey())==0);
                resulting_collapsed_data[col_dt.toKey()] = col_dt;
            }
            else {
                assert(collapse_data[collapse_index] == 0xFFFFFFFF);
                assert(collapse_data[collapse_index+1] == 0xFFFFFFFF);
#ifdef TRACK_COLLAPSED
                uint32_t counter = collapse_data_counts[collapse_index/8];
                counter >>= 8 * ((collapse_index/2) % 4);
                counter &= 0xFF;
                assert(counter==0);
#endif
            }
        }
    }
}

void validate_collapsed_data(const std::vector<uint32_t> & bucket_threads_positions_prev,
                const std::vector<uint32_t> & bucket_threads_positions_now,
                uint32_t BUCKETS_NUM, uint32_t edge_per_bucket,
                const std::unordered_map<uint64_t, CollapsedData> & collapsed,
                const std::vector<std::vector<MemRange>> & buckets,
                MetalContext & context) {

    assert(bucket_threads_positions_prev.size() == bucket_threads_positions_now.size());
    uint32_t pos_sz = bucket_threads_positions_prev.size();
    uint32_t collapsed_num = 0;
    for (uint32_t i=0; i<pos_sz; i++) {
        assert(bucket_threads_positions_prev[i] <= bucket_threads_positions_now[i]);
        collapsed_num += bucket_threads_positions_now[i] - bucket_threads_positions_prev[i];
    }

    std::unordered_map<uint64_t, CollapsedData> collapsed_data(collapsed.begin(), collapsed.end());

    assert( bucket_threads_positions_prev.size() == BUCKETS_NUM * BUCKETS_NUM );

    for (int ui = 0; ui < BUCKETS_NUM; ui++) {
        uint32_t baseUIdx = ui*BUCKETS_NUM;
        for (int vi = 0; vi < BUCKETS_NUM; vi++) {
            uint32_t pos_idx = baseUIdx + vi;
            uint32_t prevPos = bucket_threads_positions_prev[pos_idx];
            uint32_t nowPos = bucket_threads_positions_now[pos_idx];
            assert(prevPos <= nowPos);
            if (prevPos < nowPos) {
                const uint64_t * data = (const uint64_t * ) buckets[ui][vi].get_data_ptr(context);

                for (int n = prevPos; n<nowPos; n++) {
                    CollapsedData col_dt(data[n], edge_per_bucket,
                                ui, vi);

                    // Order is controlled by the test now
                    if (collapsed_data.count(col_dt.toKey())>0) {
                        auto dt  = collapsed_data.at(col_dt.toKey());
                        assert( dt == col_dt );
                        collapsed_data.erase(col_dt.toKey());
                    }
                    else {
                         assert(false); // not found
                    }
                }
            }
        }
    }
    assert(collapsed_data.empty());
    assert(collapsed_num == collapsed.size());
}

void extract_edge_migrated_data(const std::vector<uint32_t> & bucket_threads_positions_pre,
            const std::vector<uint32_t> & bucket_threads_positions_trim,
            uint32_t EDGE_BITS,
            uint32_t BUCKETS_NUM,
            bool UPass,
            uint32_t bucket_idx,
            const std::vector<std::vector<MemRange>> & buckets,
            MetalContext & context,
            std::unordered_set<uint64_t> & result_migrated_data)
{
    assert(result_migrated_data.empty());
    assert(bucket_threads_positions_pre.size() == bucket_threads_positions_trim.size());
    assert(bucket_threads_positions_pre.size() == BUCKETS_NUM*BUCKETS_NUM);

    const uint32_t UBit = uint32_t(1) << 31;

    for (int ui = 0; ui < BUCKETS_NUM; ui++) {
        if (UPass && bucket_idx==ui)
            continue;
        for (int vi = 0; vi < BUCKETS_NUM; vi++) {
            if (!UPass && bucket_idx==vi)
                continue;
            uint32_t pos_idx = ui*BUCKETS_NUM + vi;
            uint32_t prevPos = bucket_threads_positions_pre[pos_idx];
            uint32_t nowPos = bucket_threads_positions_trim[pos_idx];
            assert(prevPos <= nowPos);
            if (prevPos < nowPos) {
                    const uint64_t * data = (const uint64_t * ) buckets[ui][vi].get_data_ptr(context);

                    for (int n = prevPos; n<nowPos; n++) {
                        uint64_t migrated_edge = data[n];
                        assert(result_migrated_data.count(migrated_edge)==0);
                        result_migrated_data.insert(migrated_edge);

                        // Let's check that migrated data is really qialify...
                        uint32_t hash1 = uint32_t(migrated_edge);
                        uint32_t hash2 = uint32_t(migrated_edge>>32);
                        bool u1 = (hash1 & UBit) != 0;
                        bool u2 = (hash2 & UBit) != 0;
                        assert(u1==u2);
                        assert(u1 == UPass);
#ifndef TRACK_COLLAPSED
                        const uint32_t edge_per_bucket = (uint32_t(1) << EDGE_BITS) / BUCKETS_NUM;
                        if (u1) {
                            uint32_t bkt1 = (hash1 ^ UBit) / edge_per_bucket;
                            uint32_t bkt2 = (hash2 ^ UBit) / edge_per_bucket;
                            assert(bkt1>=0 && bkt1<BUCKETS_NUM);
                            assert(bkt2>=0 && bkt2<BUCKETS_NUM);
                            if (bkt1 < bucket_idx) {
                                assert(bkt1<bkt2);
                                assert(bkt2 == bucket_idx);
                            }
                            else if (bkt1>bucket_idx) {
                                assert(bkt1>bkt2);
                                assert(bkt2 == bucket_idx);
                            }
                            else {
                                assert(false);
                            }
                        }
                        else {
                            // V-V first processing has2, next hash1
                            uint32_t bkt1 = hash2 / edge_per_bucket;
                            uint32_t bkt2 = hash1 / edge_per_bucket;
                            assert(bkt1>=0 && bkt1<BUCKETS_NUM);
                            assert(bkt2>=0 && bkt2<BUCKETS_NUM);
                            if (bkt1 < bucket_idx) {
                                assert(bkt1<bkt2);
                                assert(bkt2 == bucket_idx);
                            }
                            else if (bkt1>bucket_idx) {
                                assert(bkt1>bkt2);
                                assert(bkt2 == bucket_idx);
                            }
                            else {
                                assert(false);
                            }
                        }
#endif
                    }
            }
        }
    }

}

void validate_expected_collapsed(const std::unordered_map<uint64_t, CollapsedData> & collapsed_data,  // Got from  collapse_edges buffer
            const std::unordered_map<uint64_t, CollapsedData> & collapsed_res) // Got from input data (hash order might be invalid)
{
    assert(collapsed_data.size() == collapsed_res.size());
    std::unordered_map<uint64_t, CollapsedData> expected_collapsed(collapsed_res);

    for (const auto & c : collapsed_data) {
        assert(expected_collapsed.count(c.first));
        expected_collapsed.erase(c.first);
    }

    assert(expected_collapsed.empty());
}

void check_data_duplication(MetalContext & context,
            const std::vector<std::vector<MemRange>> & buckets) {

     std::unordered_set<uint64_t> test_data;
     int BUCKETS_NUM = buckets.size();

     for ( int ui=0; ui<BUCKETS_NUM; ui++) {
          const std::vector<MemRange> & urow = buckets[ui];
          assert(urow.size() == BUCKETS_NUM);
          for (int vi=0; vi<BUCKETS_NUM; vi++) {
               std::vector<uint64_t> data;
               urow[vi].copy_data_into(context, data);

               for (uint64_t d : data) {
                   if (!d)
                        continue;

                    auto r = test_data.insert(d);
                    if (!r.second) {
                        std::cout << "Data duplication issue: " << std::hex << d << std::dec << std::endl;
                    }
               }
          }
     }
}
