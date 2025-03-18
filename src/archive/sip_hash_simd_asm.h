//
// Created by Konstantin Bay on 1/30/25.
//

#ifndef SIP_HASH_SIMD_ASM_H
#define SIP_HASH_SIMD_ASM_H

#include <arm_neon.h>
#include <stdint.h>
#include <iostream>

void sip_hash_2x_asm(const uint64_t * v, uint32_t nonce1, uint32_t nonce2, uint32_t *results);

#endif //SIP_HASH_SIMD_ASM_H
