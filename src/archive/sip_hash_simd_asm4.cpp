
#include "sip_hash_simd_asm4.h"

// v0 += v1
#define ADD_VV(v0,v1)  "add " v0 ".2d, " v0 ".2d, " v1 ".2d\n"

#define ADD_VVVV(v0,v1, v2,v4) \
            ADD_VV(v0,v1) \
            ADD_VV(v2,v4)

// v1 ^= v0
#define XOR_VV(v1, v0) "eor " v1 ".16b, " v1 ".16b, " v0 ".16b\n"

#define XOR_VVVV(v1, v0, v3,v4) \
            XOR_VV(v1,v0) \
            XOR_VV(v3,v4)

#define STR(x) #x

// v1 = rotl(v1, X)
// v2 = rotl(v2, X)
// t1, t2, t3, t4 - tmp registers
#define ROTL_VVVV(v1, v2, X, X_INV, t1, t2, t3 ,t4 ) \
        "shl " t1 ".2d, " v1 ".2d, #" STR(X) "\n"           \
        "shl " t3 ".2d, " v2 ".2d, #" STR(X) "\n"           \
        "sri " t2 ".2d, " v1 ".2d, #" STR(X_INV) "\n"           \
        "sri " t4 ".2d, " v2 ".2d, #" STR(X_INV) "\n"           \
        "orr " v1 ".16b, " t1 ".16b, " t2 ".16b\n" \
        "orr " v2 ".16b, " t3 ".16b, " t4 ".16b\n"

// SipHash round macro in assembly
// Using v4,v5 tmp register
// v0,v1,v2,v3  - v values, series 1
// v4,v5 - tmp registers 1
// v6,v7,v8,v9 - v values, series 2
// v10,v11 - tmp registers 2
#define SIPROUND_ASM_VVVV(v0,v1,v2,v3, v4,v5, v6,v7,v8,v9, v10,v11)      \
        ADD_VVVV(v0,v1, v6,v7) \
        ADD_VVVV(v2,v3, v8,v9) \
        ROTL_VVVV(v1, v6, 13, 64-13, v4,v5,v10,v11) \
        ROTL_VVVV(v3, v9, 16, 64-16, v4,v5,v10,v11) \
        XOR_VVVV(v1,v0, v7,v6) \
        XOR_VVVV(v3,v2, v9,v8) \
        ROTL_VVVV(v0, v6, 32, 64-62, v4,v5,v10,v11) \
        ADD_VVVV(v2,v1, v8,v7) \
        ADD_VVVV(v0,v3, v6,v9) \
        ROTL_VVVV(v1, v7, 17, 64-17, v4,v5,v10,v11) \
        ROTL_VVVV(v3, v9, 21, 64-21, v4,v5,v10,v11) \
        XOR_VVVV(v1,v2, v7,v8) \
        XOR_VVVV(v3,v0, v9,v6) \
        ROTL_VVVV(v2, v8, 32, 64-32, v4,v5,v10,v11)

// Parallel SipHash for 2 nonces using inline assembly
void sip_hash_4x_asm(const uint64_t * v, uint32_t * nonces, uint32_t *results) {
    asm (
        // we have v0, v1,v2,v3 - hash 1,2 constants
        // v4,v5 - hash 1,2, ROUND tmp values
        // v6 - nonces 1,2
        // v7 - nonces 3,4
        // v8, v9,v10,v11 - hash 3,4 constants
        // v12,v13 - hash 3,4, ROUND tmp values
        // v15 - 0xFF values

        "mov x0, %[v]\n"
        "mov x5, %[nonces]\n"
        "mov x10, %[results]\n"     ///

        // Load nonces into a NEON register and zero-extend to 64 bits
        //"movi v6.4s, #0\n"              // Set all in zeroes
        //"movi v7.4s, #0\n"

        "ldp w6, w7, [x5]\n" // load all 4 nonces
        "ldp w8, w9, [x5, #8]\n"

        "ldr x1, [x0]\n" // load all 4 v constants
        "ldr x2, [x0, #8]\n"
        "ldr x3, [x0, #16]\n"
        "ldr x4, [x0, #24]\n"

        // v15 using for 0xFF masks
        "mov x5, 0xff\n"

        // Load nonces into register, as well as data. Mixed because of latency
        "mov v6.d[0], x6\n"            // Load nonce1 into v6
        "mov v7.d[0], x8\n"
        "mov v6.d[1], x7\n"
        "mov v7.d[1], x9\n"

        "mov v15.d[0], x5\n"       // Set first byte to 0xff

        // Load constants into NEON registers (128-bit)
        "mov v0.d[0], x1\n"                    // Load lower 64 bits into v0
        "mov v1.d[0], x2\n"                    // Load lower 64 bits into v1
        "mov v2.d[0], x3\n"                    // Load lower 64 bits into v2
        "mov v3.d[0], x4\n"                    // Load lower 64 bits into v3
        "mov v0.d[1], x1\n"                    // Load upper 64 bits into v0
        "mov v1.d[1], x2\n"                    // Load upper 64 bits into v1
        "mov v2.d[1], x3\n"                    // Load upper 64 bits into v2
        "mov v3.d[1], x4\n"                    // Load upper 64 bits into v3

        "mov v15.d[1], x5\n"       // Set 9th byte to 0xff

        // Copy v0-v3 to v8-v11
        "mov v8.16b, v0.16b\n"
        "mov v9.16b, v1.16b\n"
        "mov v10.16b,v2.16b\n"
        "mov v11.16b,v3.16b\n"

        // XOR v3 with nonces
        "eor v3.16b, v3.16b, v6.16b\n"
        "eor v11.16b, v11.16b, v7.16b\n"

        // Perform SipHash rounds
        SIPROUND_ASM_VVVV("v0","v1","v2","v3", "v4","v5", "v8","v9","v10","v11", "v12","v13")
        SIPROUND_ASM_VVVV("v0","v1","v2","v3", "v4","v5", "v8", "v9","v10","v11", "v12","v13")


        // XOR v0 with nonces
        "eor v0.16b, v0.16b, v6.16b\n"
        "eor v8.16b, v8.16b, v7.16b\n"


        // XOR v2 with 0xff
        "eor v2.16b, v2.16b, v15.16b\n"
        "eor v9.16b, v2.16b, v15.16b\n"

        // Final rounds
        SIPROUND_ASM_VVVV("v0","v1","v2","v3", "v4","v5", "v8", "v9","v10","v11", "v12","v13")
        SIPROUND_ASM_VVVV("v0","v1","v2","v3", "v4","v5", "v8", "v9","v10","v11", "v12","v13")
        SIPROUND_ASM_VVVV("v0","v1","v2","v3", "v4","v5", "v8", "v9","v10","v11", "v12","v13")
        SIPROUND_ASM_VVVV("v0","v1","v2","v3", "v4","v5", "v8", "v9","v10","v11", "v12","v13")

        // Final XOR and store results
        "eor v4.16b, v0.16b, v1.16b\n"
        "eor v5.16b, v2.16b, v3.16b\n"

        "eor v12.16b, v8.16b, v9.16b\n"
        "eor v13.16b, v10.16b, v11.16b\n"

        "eor v1.16b, v4.16b, v5.16b\n"
        "eor v8.16b, v12.16b, v13.16b\n"

        // Store results
        "umov w1, v1.s[0]\n"
        "umov w2, v1.s[3]\n"
        "umov w3, v8.s[0]\n"
        "umov w4, v8.s[3]\n"
        "str w1, [x10]\n"        // Store w3 at results[0]
        "str w2, [x10, #4]\n"
        "str w3, [x10, #8]\n"
        "str w4, [x10, #12]\n"
        :
        : [v] "r"(v), [nonces] "r"(nonces), [results] "r"(results)
        : "v0", "v1", "v2", "v3", "v4", "v5", "v6", "x0", "x1", "x3", "x4", "x5", "x6","x7", "x8", "x9", "x10", "memory"
    );
}