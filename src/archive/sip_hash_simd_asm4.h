//
// Created by Konstantin Bay on 1/31/25.
//

#ifndef SIP_HASH_SIMD_ASM4_H
#define SIP_HASH_SIMD_ASM4_H


#include <arm_neon.h>
#include <stdint.h>
#include <iostream>

void sip_hash_4x_asm(const uint64_t * v, uint32_t * nonces, uint32_t *results);


#endif //SIP_HASH_SIMD_ASM4_H
