//
// Created by Konstantin Bay on 1/31/25.
//

#ifndef SIP_HASH_SIMD8_H
#define SIP_HASH_SIMD8_H


#include <arm_neon.h>
#include <stdint.h>

void sip_hash_8x(const uint64_t * v, uint64x2_t nonces[4], uint32_t * results);


#endif //SIP_HASH_SIMD8_H
