#ifndef TEST_ST1_BUILD_BUCKETS_U_H
#define TEST_ST1_BUILD_BUCKETS_U_H

#include "../miner/metal.h"

void test_st1_build_buckets(
            const uint64_t * v,
            MetalContext & context,
            uint32_t bucket_allocated_size,
            const std::vector<MemRange> & st1_8B_buckets,
            const std::vector<MemRange> & st1_1B_buckets,
            uint EDGE_BITS);


#endif //TEST_ST1_BUILD_BUCKETS_U_H
