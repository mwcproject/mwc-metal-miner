//
// Created by Konstantin Bay on 1/31/25.
//

#include "sip_hash_simd4.h"


// SipHash round macro (vectorized for NEON)
// uint64x2_t vaddq_u64( uint64x2_t a,uint64x2_t b ), Addition            ADD Vd.2D,Vn.2D,Vm.2D
// uint64x2_t vorrq_u64(uint64x2_t a, uint64x2_t b);  Bitwise OR          VORR q0,q0,q0
// uint64x2_t vshlq_n_u64(uint64x2_t a,const int n);  Vector shift left   SHL Vd.2D,Vn.2D,#n
// uint64x2_t vshrq_n_u64(uint64x2_t a,const int n);  Vector shift right  USHR Vd.2D,Vn.2D,#n
// uint64x2_t veorq_u64(uint64x2_t a,uint64x2_t b);   XOR                 EOR Vd.16B,Vn.16B,Vm.16B

#define SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v12, v13)          \
    v0 = vaddq_u64(v0, v1);                    \
    v2 = vaddq_u64(v2, v3);                    \
    v10 = vaddq_u64(v10, v11);                    \
    v12 = vaddq_u64(v12, v13);                    \
    v1 = vorrq_u64(vshlq_n_u64(v1, 13), vshrq_n_u64(v1, 64 - 13)); \
    v3 = vorrq_u64(vshlq_n_u64(v3, 16), vshrq_n_u64(v3, 64 - 16)); \
    v11 = vorrq_u64(vshlq_n_u64(v11, 13), vshrq_n_u64(v11, 64 - 13)); \
    v13 = vorrq_u64(vshlq_n_u64(v13, 16), vshrq_n_u64(v13, 64 - 16)); \
    v1 = veorq_u64(v1, v0);                    \
    v3 = veorq_u64(v3, v2);                    \
    v11 = veorq_u64(v11, v10);                    \
    v13 = veorq_u64(v13, v12);                    \
    v0 = vorrq_u64(vshlq_n_u64(v0, 32), vshrq_n_u64(v0, 64 - 32)); \
    v10 = vorrq_u64(vshlq_n_u64(v10, 32), vshrq_n_u64(v10, 64 - 32)); \
    v2 = vaddq_u64(v2, v1);                    \
    v0 = vaddq_u64(v0, v3);                    \
    v12 = vaddq_u64(v12, v11);                    \
    v10 = vaddq_u64(v10, v13);                    \
    v1 = vorrq_u64(vshlq_n_u64(v1, 17), vshrq_n_u64(v1, 64 - 17)); \
    v3 = vorrq_u64(vshlq_n_u64(v3, 21), vshrq_n_u64(v3, 64 - 21)); \
    v11 = vorrq_u64(vshlq_n_u64(v11, 17), vshrq_n_u64(v11, 64 - 17)); \
    v13 = vorrq_u64(vshlq_n_u64(v13, 21), vshrq_n_u64(v13, 64 - 21)); \
    v1 = veorq_u64(v1, v2);                    \
    v3 = veorq_u64(v3, v0);                    \
    v11 = veorq_u64(v11, v12);                    \
    v13 = veorq_u64(v13, v10);                    \
    v2 = vorrq_u64(vshlq_n_u64(v2, 32), vshrq_n_u64(v2, 64 - 32)); \
    v12 = vorrq_u64(vshlq_n_u64(v12, 32), vshrq_n_u64(v12, 64 - 32));

// Parallel SipHash for 2 nonces
void sip_hash_4x(const uint64_t * v, uint32_t nonces[4], uint32_t * results) {
    // Load constants into NEON registers
    uint64x2_t v0 = vdupq_n_u64(v[0]); // Load all lanes of vector to the same literal value, // VMOV d0,r0,r0
    uint64x2_t v1 = vdupq_n_u64(v[1]);
    uint64x2_t v2 = vdupq_n_u64(v[2]);
    uint64x2_t v3 = vdupq_n_u64(v[3]);

    uint64x2_t v10 = v0;
    uint64x2_t v11 = v1;
    uint64x2_t v12 = v2;
    uint64x2_t v13 = v3;

    // Load nonces into a NEON register and zero-extend to 64 bits
    // vld1_u32:  Load a single vector from memory, VLD1.32 {d0}, [r0]
    // vmovl_u32  Vector long move,   // VMOVL.U32 q0,d0
    uint64x2_t nonce1 = vmovl_u32(vld1_u32(nonces));
    uint64x2_t nonce2 = vmovl_u32(vld1_u32(&nonces[2]));

    // XOR v3 with nonces
    v3 = veorq_u64(v3, nonce1); // Bitwise exclusive or (EOR or XOR), // VEOR q0,q0,q0
    v13 = veorq_u64(v13, nonce2);

    // Perform SipHash rounds
    SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v11, v12);
    SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v11, v12);

    // XOR v0 with nonces
    v0 = veorq_u64(v0, nonce1); // XOR
    v10 = veorq_u64(v10, nonce2);

    // XOR v2 with 0xff
    v2 = veorq_u64(v2, vdupq_n_u64(0xff));
    v12 = veorq_u64(v12, vdupq_n_u64(0xff));

    // Final rounds
    SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v11, v12);
    SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v11, v12);
    SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v11, v12);
    SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v11, v12);

    // Final XOR and extract results
    uint64x2_t final1 = veorq_u64(veorq_u64(v0, v1), veorq_u64(v2, v3));
    uint64x2_t final2 = veorq_u64(veorq_u64(v10, v11), veorq_u64(v12, v13));

    // Store results
    uint32_t *final_ptr1 = (uint32_t *)&final1;
    results[0] = final_ptr1[0];
    results[1] = final_ptr1[1];
    uint32_t *final_ptr2 = (uint32_t *)&final2;
    results[2] = final_ptr2[0];
    results[3] = final_ptr2[1];
}