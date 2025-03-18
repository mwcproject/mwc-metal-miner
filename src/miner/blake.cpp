// Copyright 2025 The MWC Developers
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "blake.h"
#include <vector>
#include <assert.h>
#include "cuckatoo.h"


// BLAKE2b constants
constexpr uint64_t IV[8] = {
    0x6a09e667f3bcc908, 0xbb67ae8584caa73b,
    0x3c6ef372fe94f82b, 0xa54ff53a5f1d36f1,
    0x510e527fade682d1, 0x9b05688c2b3e6c1f,
    0x1f83d9abfb41bd6b, 0x5be0cd19137e2179
};

constexpr uint8_t SIGMA[12][16] = {
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
    {14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3},
    {11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4},
    {7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8},
    {9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13},
    {2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9},
    {12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11},
    {13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10},
    {6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5},
    {10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0},
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
    {14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3}
};

// Rotate right operation
inline uint64_t rotr64(uint64_t x, uint8_t n) {
    return (x >> n) | (x << (64 - n));
}

// BLAKE2b mixing function
inline void G(uint64_t& a, uint64_t& b, uint64_t& c, uint64_t& d, uint64_t x, uint64_t y) {
    a = a + b + x;
    d = rotr64(d ^ a, 32);
    c = c + d;
    b = rotr64(b ^ c, 24);
    a = a + b + y;
    d = rotr64(d ^ a, 16);
    c = c + d;
    b = rotr64(b ^ c, 63);
}

// BLAKE2b compression function
void compress(uint64_t* h, const uint8_t* block, uint64_t t, bool last) {
    uint64_t v[16];
    for (int i = 0; i < 8; ++i) v[i] = h[i];
    for (int i = 0; i < 8; ++i) v[i + 8] = IV[i];

    v[12] ^= t; // Low 64 bits of the offset
    //v[13] ^= 0; // High 64 bits of the offset (always 0 for BLAKE2b)
    if (last) v[14] ^= ~0ULL;

    uint64_t m[16];
    for (int i = 0; i < 16; ++i) {
        m[i] = *reinterpret_cast<const uint64_t*>(block + i * 8);
    }

    for (int i = 0; i < 12; ++i) {
        G(v[0], v[4], v[8], v[12], m[SIGMA[i][0]], m[SIGMA[i][1]]);
        G(v[1], v[5], v[9], v[13], m[SIGMA[i][2]], m[SIGMA[i][3]]);
        G(v[2], v[6], v[10], v[14], m[SIGMA[i][4]], m[SIGMA[i][5]]);
        G(v[3], v[7], v[11], v[15], m[SIGMA[i][6]], m[SIGMA[i][7]]);
        G(v[0], v[5], v[10], v[15], m[SIGMA[i][8]], m[SIGMA[i][9]]);
        G(v[1], v[6], v[11], v[12], m[SIGMA[i][10]], m[SIGMA[i][11]]);
        G(v[2], v[7], v[8], v[13], m[SIGMA[i][12]], m[SIGMA[i][13]]);
        G(v[3], v[4], v[9], v[14], m[SIGMA[i][14]], m[SIGMA[i][15]]);
    }

    for (int i = 0; i < 8; ++i) h[i] ^= v[i] ^ v[i + 8];
}

// BLAKE2b hash function
void blake2b(uint8_t* out, const uint8_t* in, size_t inlen) {
    uint64_t h[8];
    for (int i = 0; i < 8; ++i) h[i] = IV[i];

    // Parameter block: digest length = 32, key length = 0, fanout = 1, depth = 1
    h[0] ^= 0x01010000 ^ 32;

    uint64_t t = 0;
    const uint8_t* block = in;
    size_t remaining = inlen;

    while (remaining > 128) {
        t += 128;
        compress(h, block, t, false);
        block += 128;
        remaining -= 128;
    }

    uint8_t last_block[128] = {0};
    memcpy(last_block, block, remaining);
    compress(h, last_block, t + remaining, true);

    memcpy(out, h, 32);
}

