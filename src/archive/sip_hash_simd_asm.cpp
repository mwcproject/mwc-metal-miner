//
// Created by Konstantin Bay on 1/30/25.
//

#include "sip_hash_simd_asm.h"


// SipHash round macro in assembly
// Using v4,v5 tmp register
#define SIPROUND_ASM      \
    "add v0.2d, v0.2d, v1.2d\n"       \
    "add v2.2d, v2.2d, v3.2d\n"       \
    "shl v4.2d, v1.2d, #13\n"           \
    "sri v5.2d, v1.2d, #51\n"           \
    "orr v1.16b, v4.16b, v5.16b\n"      \
    "shl v4.2d, v3.2d, #16\n"           \
    "sri v5.2d, v3.2d, #48\n"           \
    "orr v3.16b, v4.16b, v5.16b\n"      \
    "eor v1.16b, v1.16b, v0.16b\n"    \
    "eor v3.16b, v3.16b, v2.16b\n"    \
    "shl v4.2d, v0.2d, #32\n"           \
    "sri v5.2d, v0.2d, #32\n"           \
    "orr v0.16b, v4.16b, v5.16b\n"      \
    "add v2.2d, v2.2d, v1.2d\n"       \
    "add v0.2d, v0.2d, v3.2d\n"       \
    "shl v4.2d, v1.2d, #17\n"           \
    "sri v5.2d, v1.2d, #47\n"           \
    "orr v1.16b, v4.16b, v5.16b\n"      \
    "shl v4.2d, v3.2d, #21\n"           \
    "sri v5.2d, v3.2d, #43\n"           \
    "orr v3.16b, v4.16b, v5.16b\n"      \
    "eor v1.16b, v1.16b, v2.16b\n"    \
    "eor v3.16b, v3.16b, v0.16b\n"    \
    "shl v4.2d, v2.2d, #32\n"           \
    "sri v5.2d, v2.2d, #32\n"           \
    "orr v2.16b, v4.16b, v5.16b\n"

// Parallel SipHash for 2 nonces using inline assembly
void sip_hash_2x_asm(const uint64_t * v, uint32_t nonce1, uint32_t nonce2, uint32_t *results) {
    asm volatile(
        // we have v0, v1,v2,v3 - hash constants
        // v4,v5 - ROUND tmp values
        // v6 - nonces (tmp value)

        "mov x0, %[results]\n"     ///
        "ldr x1, [x0]\n" // load all 4 v constants
        "ldr x2, [x0, #8]\n"
        "ldr x3, [x0, #16]\n"
        "ldr x4, [x0, #24]\n"

        // Load nonces into a NEON register and zero-extend to 64 bits
        "movi v6.4s, #0\n"              // Set all in zeroes
        "ins v6.s[0], %w[nonce1]\n"            // Load nonce1 into v6
        "mov v6.s[2], %w[nonce2]\n"            // Load nonce2 into v6

        // Load constants into NEON registers (128-bit)
        "mov v0.d[0], x1\n"                    // Load lower 64 bits into v0
        "mov v0.d[1], x1\n"                    // Load upper 64 bits into v0

        "mov v1.d[0], x2\n"                    // Load lower 64 bits into v1
        "mov v1.d[1], x2\n"                    // Load upper 64 bits into v1

        "mov v2.d[0], x3\n"                    // Load lower 64 bits into v2
        "mov v2.d[1], x3\n"                    // Load upper 64 bits into v2

        "mov v3.d[0], x4\n"                    // Load lower 64 bits into v3
        "mov v3.d[1], x4\n"                    // Load upper 64 bits into v3

        // XOR v3 with nonces
        "eor v3.16b, v3.16b, v6.16b\n"

        // Perform SipHash rounds
        SIPROUND_ASM
        SIPROUND_ASM

        // XOR v0 with nonces
        "eor v0.16b, v0.16b, v4.16b\n"

        // XOR v2 with 0xff
        "movi v5.2d, #0xff\n"
        "eor v2.16b, v2.16b, v5.16b\n"

        "mov x3, %[results]\n"     ///

        // Final rounds
        SIPROUND_ASM
        SIPROUND_ASM
        SIPROUND_ASM
        SIPROUND_ASM

        // Final XOR and store results
        "eor v4.16b, v0.16b, v1.16b\n"
        "eor v5.16b, v2.16b, v3.16b\n"
        "eor v1.16b, v4.16b, v5.16b\n"

        // Store results
        "umov w0, v1.s[0]\n"
        "umov w1, v1.s[3]\n"
        "str w0, [x3]\n"        // Store w3 at results[0]
        "str w1, [x3, #4]\n"    // Store w4 at results[1]
        :
        : [v] "r"(v), [nonce1] "r"(nonce1), [nonce2] "r"(nonce2), [results] "r"(results)
        : "v0", "v1", "v2", "v3", "v4", "v5", "v6", "x0", "x1", "x3", "x4", "memory"
    );
}