#ifndef TEST_BUILD_MASK_H
#define TEST_BUILD_MASK_H

#include "../miner/metal.h"
#include <unordered_set>

struct CollapsedData {
    uint32_t hash1 = 0xFFFFFFFF;
    uint32_t hash2 = 0xFFFFFFFF;
    uint32_t counter = 0xFFFFFFFF;

    CollapsedData() = default;
    CollapsedData(const CollapsedData &) = default;
    CollapsedData & operator = (const CollapsedData &) = default;
    inline bool operator == (const CollapsedData & other) const {
        return toKey()==other.toKey() && counter==other.counter;
    }

    CollapsedData(uint32_t _hash1,  uint32_t _hash2, uint32_t _counter) :
        hash1(_hash1), hash2(_hash2), counter(_counter)
    {
        uint32_t ubit = uint32_t(1) << 31;
        bool u1 = (hash1 & ubit)!=0;
        bool u2 = (hash2 & ubit)!=0;

        if (u1==u2) {
            assert(counter%2==1);
        }
        else {
            assert(counter%2==0);
        }
    }

    CollapsedData(uint64_t edges, uint32_t edge_per_bucket,
                uint32_t ubkt, uint32_t vbkt);

    inline uint64_t toKey() const {
        uint32_t h1 = std::min(hash1, hash2);
        uint32_t h2 = std::max(hash1, hash2);
        return uint64_t(h1) | (uint64_t(h2) << 32);
    }
};

void test_build_buckets_st2(void * mask_data, uint32_t mask_size,
            uint32_t in_bucket_idx,
            const std::vector<uint32_t> & bucket_positions,
            const std::vector< std::vector<uint64_t>> & in,
            const std::vector<std::vector<uint64_t>> & out,
            uint32_t EDGE_BITS, uint32_t buckets_num,
            const uint64_t v[4],
            const std::unordered_map<uint64_t, CollapsedData> & collapsed_data);

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
            std::unordered_map<uint64_t, CollapsedData>  &collapsed_res,
            uint32_t COLLAPSE_STASH_COUNTER,
            uint32_t collapse_stash_prev_size,
            uint32_t * collapse_stash);

uint32_t validate_bucket_data( bool hasUVbit, uint32_t EDGE_BITS, uint32_t buckets_num,
                    int bucketU, int bucketV, uint32_t bucket_position,
                    const std::vector<uint64_t> & data );

// Simply check if all the data is fill with zeros
void validate_empty_mask_buffers( void * mask_buffer, uint32_t mask_buffer_size );

void validate_empty_collapse_buffer( void* collapse_buffer, uint32_t collapse_buffer_size,
#ifdef TRACK_COLLAPSED
                void* collapse_edges_counts, uint32_t collapse_edges_counts_size,
#endif
                uint32_t edge_per_bucket );

void validate_filled_collapse_buffer( void * collapse_buffer_data, uint32_t collapse_buffer_size,
#ifdef TRACK_COLLAPSED
            void * collapse_edges_counts_data, uint32_t collapse_edges_counts_size,
#endif
            uint32_t edge_per_bucket,
            uint32_t EDGE_BITS, uint32_t CYCLE_LEN, uint32_t buckets_num,
            void* mask_data, uint32_t mask_size,
            std::unordered_map<uint64_t, CollapsedData> & resulting_collapsed_data);

void validate_collapsed_data(const std::vector<uint32_t> & bucket_threads_positions_prev,
                const std::vector<uint32_t> & bucket_threads_positions_now,
                uint32_t BUCKETS_NUM, uint32_t edge_per_bucket,
                const std::unordered_map<uint64_t, CollapsedData> & collapsed,
                const std::vector<std::vector<MemRange>> & buckets,
                MetalContext & context);

void validate_expected_collapsed(const std::unordered_map<uint64_t, CollapsedData> & collapsed_data, // Got from  collapse_edges buffer
            const std::unordered_map<uint64_t, CollapsedData> & collapsed_res);

void extract_edge_migrated_data(const std::vector<uint32_t> & bucket_threads_positions_pre,
            const std::vector<uint32_t> & bucket_threads_positions_trim,
            uint32_t EDGE_BITS,
            uint32_t BUCKETS_NUM,
            bool UPass,
            uint32_t bucket_idx,
            const std::vector<std::vector<MemRange>> & buckets,
            MetalContext & context,
            std::unordered_set<uint64_t> & result_migrated_data);

void check_data_duplication(MetalContext & context,
            const std::vector<std::vector<MemRange>> & buckets);

#endif //TEST_BUILD_MASK_H
