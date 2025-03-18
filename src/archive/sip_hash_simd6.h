//
// Created by Konstantin Bay on 1/31/25.
//

#ifndef SIP_HASH_SIMD6_H
#define SIP_HASH_SIMD6_H



#include <arm_neon.h>
#include <stdint.h>

void sip_hash_6x(const uint64_t * v, uint32_t * nonces, uint32_t * results4);


#endif //SIP_HASH_SIMD6_H
