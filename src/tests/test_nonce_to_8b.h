#ifndef TEST_NONCE_TO_8B_H
#define TEST_NONCE_TO_8B_H

#include "../miner/metal.h"

void test_nonce_to_8b_check(const uint64_t * v,
            MetalContext & context,
            const MemRange & in,
            const MemRange & out,
            uint EDGE_BITS);

#endif //TEST_NONCE_TO_8B_H
