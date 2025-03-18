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

#ifndef DATA_STORAGE_H
#define DATA_STORAGE_H

#include <assert.h>
#include <cstdint>
#include <iostream>
#include <span>
#include <vector>
#include "nonce_hash.h"

template <uint32_t ELEMENT_SIZE, uint8_t HASH_BITS>
inline HashNonce extract_hash_nonce(uint64_t * data, uint32_t idx) {
    size_t byte_pos = idx * ELEMENT_SIZE;       // Compute byte position
    size_t ulong_index = byte_pos / 8;        // Find which ulong to store in
    size_t bit_offset  = (byte_pos % 8) * 8;  // Bit offset inside ulong

    // Fetch existing value
    uint64_t lower = data[ulong_index];

    // Store in one or two ulong slots
    lower >>= bit_offset;  // Set new value

    // If hash crosses into the next ulong, store remaining bits
    if (bit_offset > (64 - ELEMENT_SIZE*8)) {
        uint64_t upper = data[ulong_index + 1];
        upper <<= (64 - bit_offset);  // Clear existing bits
        lower |= upper;
    }

    // Applly full mask, so res.nonce will be clear
    lower &= (uint64_t(1) << (ELEMENT_SIZE*8)) - 1;

    HashNonce res;
    res.hash = lower & ((uint64_t(1)<<HASH_BITS)-1);
    res.nonce = (lower >> HASH_BITS);
    return res;
}

template <uint32_t ELEMENT_SIZE, uint8_t HASH_BITS>
inline void store_hash_nonce(uint64_t * data, uint32_t idx, const HashNonce & hn) {
    uint64_t dt = uint64_t(hn.hash) | (uint64_t(hn.nonce) << HASH_BITS);

    size_t byte_pos = idx * ELEMENT_SIZE;       // Compute byte position
    size_t ulong_index = byte_pos / 8;        // Find which ulong to store in
    size_t bit_offset  = (byte_pos % 8) * 8;  // Bit offset inside ulong

    uint64_t mask = (uint64_t(1) << (ELEMENT_SIZE*8)) - 1;

    // Fetch existing value
    uint64_t lower = data[ulong_index];

    // Store in one or two ulong slots
    lower &= ~(mask << bit_offset);  // Clear existing bits
    lower |= (dt & mask) << bit_offset;  // Set new value
    data[ulong_index] = lower;

    // If hash crosses into the next ulong, store remaining bits
    if (bit_offset > (64 - ELEMENT_SIZE*8)) {
        uint64_t upper = data[ulong_index + 1];
        upper &= ~(mask >> (64 - bit_offset));  // Clear existing bits
        upper |= (dt & mask) >> (64 - bit_offset);  // Set new value
        data[ulong_index + 1] = upper;
    }

}


#endif //DATA_STORAGE_H
