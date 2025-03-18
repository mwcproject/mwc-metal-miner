//
// Created by Konstantin Bay on 1/29/25.
//

#include "sip_hash_simd.h"

// SipHash round macro (vectorized for NEON)
// uint64x2_t vaddq_u64( uint64x2_t a,uint64x2_t b ), Addition            ADD Vd.2D,Vn.2D,Vm.2D
// uint64x2_t vorrq_u64(uint64x2_t a, uint64x2_t b);  Bitwise OR          VORR q0,q0,q0
// uint64x2_t vshlq_n_u64(uint64x2_t a,const int n);  Vector shift left   SHL Vd.2D,Vn.2D,#n
// uint64x2_t vshrq_n_u64(uint64x2_t a,const int n);  Vector shift right  USHR Vd.2D,Vn.2D,#n
// uint64x2_t veorq_u64(uint64x2_t a,uint64x2_t b);   XOR                 EOR Vd.16B,Vn.16B,Vm.16B

#define SIPROUND_NEON(v0, v1, v2, v3)          \
    v0 = vaddq_u64(v0, v1);                    \
    v2 = vaddq_u64(v2, v3);                    \
    v1 = vorrq_u64(vshlq_n_u64(v1, 13), vshrq_n_u64(v1, 64 - 13)); \
    v3 = vorrq_u64(vshlq_n_u64(v3, 16), vshrq_n_u64(v3, 64 - 16)); \
    v1 = veorq_u64(v1, v0);                    \
    v3 = veorq_u64(v3, v2);                    \
    v0 = vorrq_u64(vshlq_n_u64(v0, 32), vshrq_n_u64(v0, 64 - 32)); \
    v2 = vaddq_u64(v2, v1);                    \
    v0 = vaddq_u64(v0, v3);                    \
    v1 = vorrq_u64(vshlq_n_u64(v1, 17), vshrq_n_u64(v1, 64 - 17)); \
    v3 = vorrq_u64(vshlq_n_u64(v3, 21), vshrq_n_u64(v3, 64 - 21)); \
    v1 = veorq_u64(v1, v2);                    \
    v3 = veorq_u64(v3, v0);                    \
    v2 = vorrq_u64(vshlq_n_u64(v2, 32), vshrq_n_u64(v2, 64 - 32));

// Parallel SipHash for 2 nonces
void sip_hash_2x(const uint64_t * v, uint32_t nonces[2], uint32_t * results) {
    // Load constants into NEON registers
    uint64x2_t v0 = vdupq_n_u64(v[0]); // Load all lanes of vector to the same literal value, // VMOV d0,r0,r0
    uint64x2_t v1 = vdupq_n_u64(v[1]);
    uint64x2_t v2 = vdupq_n_u64(v[2]);
    uint64x2_t v3 = vdupq_n_u64(v[3]);

    // Load nonces into a NEON register and zero-extend to 64 bits
    // vld1_u32:  Load a single vector from memory, VLD1.32 {d0}, [r0]
    // vmovl_u32  Vector long move,   // VMOVL.U32 q0,d0
    uint64x2_t nonce = vmovl_u32(vld1_u32(nonces));

    // XOR v3 with nonces
    v3 = veorq_u64(v3, nonce); // Bitwise exclusive or (EOR or XOR), // VEOR q0,q0,q0

    // Perform SipHash rounds
    SIPROUND_NEON(v0, v1, v2, v3);
    SIPROUND_NEON(v0, v1, v2, v3);

    // XOR v0 with nonces
    v0 = veorq_u64(v0, nonce); // XOR

    // XOR v2 with 0xff
    v2 = veorq_u64(v2, vdupq_n_u64(0xff));

    // Final rounds
    SIPROUND_NEON(v0, v1, v2, v3);
    SIPROUND_NEON(v0, v1, v2, v3);
    SIPROUND_NEON(v0, v1, v2, v3);
    SIPROUND_NEON(v0, v1, v2, v3);

    // Final XOR and extract results
    uint64x2_t final = veorq_u64(veorq_u64(v0, v1), veorq_u64(v2, v3));

    // Store results
    uint32_t *final_ptr = (uint32_t *)&final;
    results[0] = final_ptr[0];
    results[1] = final_ptr[1];
}