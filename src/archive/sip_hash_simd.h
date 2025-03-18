#ifndef SIP_HASH_SIMD_H
#define SIP_HASH_SIMD_H

#include <arm_neon.h>
#include <stdint.h>


void sip_hash_2x(const uint64_t * v, uint32_t nonces[2], uint32_t * results2);


#endif //SIP_HASH_SIMD_H
