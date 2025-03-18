//
// Created by Konstantin Bay on 1/31/25.
//

#include "sip_hash_simd8.h"


#define SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v12, v13, v20, v21, v22, v23, v30, v31, v32, v33)          \
    v0 = vaddq_u64(v0, v1);                    \
    v2 = vaddq_u64(v2, v3);                    \
    v10 = vaddq_u64(v10, v11);                    \
    v12 = vaddq_u64(v12, v13);                    \
    v20 = vaddq_u64(v20, v21);                    \
    v22 = vaddq_u64(v22, v23);                    \
    v30 = vaddq_u64(v30, v31);                    \
    v32 = vaddq_u64(v32, v33);                    \
    v1 = vorrq_u64(vshlq_n_u64(v1, 13), vshrq_n_u64(v1, 64 - 13)); \
    v3 = vorrq_u64(vshlq_n_u64(v3, 16), vshrq_n_u64(v3, 64 - 16)); \
    v11 = vorrq_u64(vshlq_n_u64(v11, 13), vshrq_n_u64(v11, 64 - 13)); \
    v13 = vorrq_u64(vshlq_n_u64(v13, 16), vshrq_n_u64(v13, 64 - 16)); \
    v21 = vorrq_u64(vshlq_n_u64(v21, 13), vshrq_n_u64(v21, 64 - 13)); \
    v23 = vorrq_u64(vshlq_n_u64(v23, 16), vshrq_n_u64(v23, 64 - 16)); \
    v31 = vorrq_u64(vshlq_n_u64(v31, 13), vshrq_n_u64(v31, 64 - 13)); \
    v33 = vorrq_u64(vshlq_n_u64(v33, 16), vshrq_n_u64(v33, 64 - 16)); \
    v1 = veorq_u64(v1, v0);                    \
    v3 = veorq_u64(v3, v2);                    \
    v11 = veorq_u64(v11, v10);                    \
    v13 = veorq_u64(v13, v12);                    \
    v21 = veorq_u64(v21, v20);                    \
    v23 = veorq_u64(v23, v22);                    \
    v31 = veorq_u64(v31, v30);                    \
    v33 = veorq_u64(v33, v32);                    \
    v0 = vorrq_u64(vshlq_n_u64(v0, 32), vshrq_n_u64(v0, 64 - 32)); \
    v10 = vorrq_u64(vshlq_n_u64(v10, 32), vshrq_n_u64(v10, 64 - 32)); \
    v20 = vorrq_u64(vshlq_n_u64(v20, 32), vshrq_n_u64(v20, 64 - 32)); \
    v30 = vorrq_u64(vshlq_n_u64(v30, 32), vshrq_n_u64(v30, 64 - 32)); \
    v2 = vaddq_u64(v2, v1);                    \
    v0 = vaddq_u64(v0, v3);                    \
    v12 = vaddq_u64(v12, v11);                    \
    v10 = vaddq_u64(v10, v13);                    \
    v22 = vaddq_u64(v22, v21);                    \
    v20 = vaddq_u64(v20, v23);                    \
    v32 = vaddq_u64(v32, v31);                    \
    v30 = vaddq_u64(v30, v33);                    \
    v1 = vorrq_u64(vshlq_n_u64(v1, 17), vshrq_n_u64(v1, 64 - 17)); \
    v3 = vorrq_u64(vshlq_n_u64(v3, 21), vshrq_n_u64(v3, 64 - 21)); \
    v11 = vorrq_u64(vshlq_n_u64(v11, 17), vshrq_n_u64(v11, 64 - 17)); \
    v13 = vorrq_u64(vshlq_n_u64(v13, 21), vshrq_n_u64(v13, 64 - 21)); \
    v21 = vorrq_u64(vshlq_n_u64(v21, 17), vshrq_n_u64(v21, 64 - 17)); \
    v23 = vorrq_u64(vshlq_n_u64(v23, 21), vshrq_n_u64(v23, 64 - 21)); \
    v31 = vorrq_u64(vshlq_n_u64(v31, 17), vshrq_n_u64(v31, 64 - 17)); \
    v33 = vorrq_u64(vshlq_n_u64(v33, 21), vshrq_n_u64(v33, 64 - 21)); \
    v1 = veorq_u64(v1, v2);                    \
    v3 = veorq_u64(v3, v0);                    \
    v11 = veorq_u64(v11, v12);                    \
    v13 = veorq_u64(v13, v10);                    \
    v21 = veorq_u64(v21, v22);                    \
    v23 = veorq_u64(v23, v20);                    \
    v31 = veorq_u64(v31, v32);                    \
    v33 = veorq_u64(v33, v30);                    \
    v2 = vorrq_u64(vshlq_n_u64(v2, 32), vshrq_n_u64(v2, 64 - 32)); \
    v12 = vorrq_u64(vshlq_n_u64(v12, 32), vshrq_n_u64(v12, 64 - 32)); \
    v22 = vorrq_u64(vshlq_n_u64(v22, 32), vshrq_n_u64(v22, 64 - 32)); \
    v32 = vorrq_u64(vshlq_n_u64(v32, 32), vshrq_n_u64(v32, 64 - 32));

// Parallel SipHash for 2 nonces
void sip_hash_8x(const uint64_t * v, uint64x2_t nonces[4], uint32_t * results) {
    // Load constants into NEON registers
    uint64x2_t v0 = vdupq_n_u64(v[0]); // Load all lanes of vector to the same literal value, // VMOV d0,r0,r0
    uint64x2_t v1 = vdupq_n_u64(v[1]);
    uint64x2_t v2 = vdupq_n_u64(v[2]);
    uint64x2_t v3 = vdupq_n_u64(v[3]);

    uint64x2_t v10 = v0;
    uint64x2_t v11 = v1;
    uint64x2_t v12 = v2;
    uint64x2_t v13 = v3;

    uint64x2_t v20 = v0;
    uint64x2_t v21 = v1;
    uint64x2_t v22 = v2;
    uint64x2_t v23 = v3;

    uint64x2_t v30 = v0;
    uint64x2_t v31 = v1;
    uint64x2_t v32 = v2;
    uint64x2_t v33 = v3;

    // XOR v3 with nonces
    v3 = veorq_u64(v3, nonces[0]); // Bitwise exclusive or (EOR or XOR), // VEOR q0,q0,q0
    v13 = veorq_u64(v13, nonces[1]);
    v23 = veorq_u64(v23, nonces[2]);
    v33 = veorq_u64(v33, nonces[3]);

    // Perform SipHash rounds
    SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v11, v12, v20, v22, v22, v23, v30, v31, v32, v33);
    SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v11, v12, v20, v22, v22, v23, v30, v31, v32, v33);

    // XOR v0 with nonces
    v0 = veorq_u64(v0, nonces[0]); // XOR
    v10 = veorq_u64(v10, nonces[1]);
    v20 = veorq_u64(v20, nonces[2]);

    // XOR v2 with 0xff
    uint64x2_t ff = vdupq_n_u64(0xff);
    v2 = veorq_u64(v2, ff);
    v12 = veorq_u64(v12, ff);
    v22 = veorq_u64(v22, ff);
    v32 = veorq_u64(v32, ff);

    // Final rounds
    SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v11, v12, v20, v22, v22, v23, v30, v31, v32, v33);
    SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v11, v12, v20, v22, v22, v23, v30, v31, v32, v33);
    SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v11, v12, v20, v22, v22, v23, v30, v31, v32, v33);
    SIPROUND_NEON(v0, v1, v2, v3, v10, v11, v11, v12, v20, v22, v22, v23, v30, v31, v32, v33);

    // Final XOR and extract results
    uint32x4_t final1 = vreinterpretq_u32_u64(veorq_u64(veorq_u64(v0, v1), veorq_u64(v2, v3)));
    uint32x4_t final2 = vreinterpretq_u32_u64(veorq_u64(veorq_u64(v10, v11), veorq_u64(v12, v13)));
    uint32x4_t final3 = vreinterpretq_u32_u64(veorq_u64(veorq_u64(v20, v21), veorq_u64(v22, v23)));
    uint32x4_t final4 = vreinterpretq_u32_u64(veorq_u64(veorq_u64(v30, v31), veorq_u64(v32, v33)));

    // Store results
    results[0] = vgetq_lane_u32(final1, 0);
    results[1] = vgetq_lane_u32(final1, 2);
    results[2] = vgetq_lane_u32(final2, 0);
    results[3] = vgetq_lane_u32(final2, 2);
    results[4] = vgetq_lane_u32(final3, 0);
    results[5] = vgetq_lane_u32(final3, 2);
    results[6] = vgetq_lane_u32(final4, 0);
    results[7] = vgetq_lane_u32(final4, 2);
}