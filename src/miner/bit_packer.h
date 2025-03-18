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

#ifndef BIT_PACKER_H
#define BIT_PACKER_H


#include <iostream>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <assert.h>

// Note bitpacker is expected to work in write only, then read only mode.
// read/write operations can't be mixed
// Write must be done first, than writeFlush -> reset and we can read
// It is not expected to store large number. Currently max value is 16, same as in the mask
struct BitPacker {
public:
    // Constructor
    BitPacker() : bit_position(0), buf_item(0) {}

    // Add a number to the encoder
    inline void encode_number(uint32_t num) {
        // Encode small numbers [0, COMPACT_NUM_LIMIT) in compact format
        encode_small_number(num);
    }

    // Valid for writer
    inline uint32_t get_writer_bit_pos() const {
        return buffer.size()*64+bit_position;
    }

    // Decode the next number from the buffer
    inline uint32_t decode_number() const {
        uint32_t num = decode_small_number();
        return num;
    }

    // Finish with a writing. add whatever we already have
    void write_flush() {
        if (bit_position>0) {
            buffer.push_back(buf_item);
            buf_item = 0;
            bit_position = 0;
        }
    }

    // Reset the buffer and bit position
    void reset_read() {
        buffer_index = 0;
        read_mask = 0;
        read_pos = 0;
    }

    inline uint32_t get_reader_bit_pos() const {
        return read_pos;
    }

    // Expected to called for reading
    inline void set_bit_pos(uint32_t bit_pos) const {
        if (bit_pos>0) {
            read_pos = bit_pos;
            bit_pos--;
            read_mask = uint64_t(1) << (bit_pos % 64);
            buffer_index = bit_pos / 64;
            buf_item = buffer[buffer_index];
            buffer_index++;
        }
        else {
            buffer_index = 0;
            read_mask = 0;
            read_pos = 0;
        }
    }

    inline size_t buf_size() const {return buffer.size();}
private:

    // Encode a small number (0–COMPACT_NUM_LIMIT) in compact format
    inline void encode_small_number(uint32_t num) {
        // Encode the number as a sequence of `1`s followed by a `0`
        uint32_t encodedValue = (1 << num) - 1; // Create a sequence of `num` `1`s
        write_bits(encodedValue, num + 1);
    }

    // Write a sequence of bits to the buffer
    inline  void write_bits(uint64_t value, int bits) {
        while (bits > 0) {
            assert(bit_position<64);

            // Calculate how many bits can be written to the current buffer element
            int bitsAvailable = 64 - bit_position;
            int bitsToWrite = std::min(bits, bitsAvailable);

            // Write the bits to the current buffer element
            buf_item |= (value << bit_position);

            // Update the value and bit position
            value >>= bitsToWrite;
            bit_position += bitsToWrite;
            bits -= bitsToWrite;

            // If the current buffer element is full, add a new one
            if (bit_position >= 64) {
                buffer.push_back(buf_item);
                buf_item = 0;
                bit_position = 0;
            }
        }
    }

    // Decode a small number in compact format
    inline uint32_t decode_small_number() const {
        uint32_t num = 0;

        // Count the number of consecutive `1`s
        while (read_next_bit()) {
            ++num;
        }

        return num;
    }

    // Read the next bit and advance the bit position
    inline bool read_next_bit() const {
        read_mask <<= 1;
        read_pos++;
        if (read_mask==0) {
            read_mask = 1;
            buf_item = buffer[buffer_index];
            buffer_index++;
        }

        return (buf_item & read_mask) != 0;
    }

private:
    std::vector<uint64_t> buffer; // The packed buffer
    int bit_position;              // Current bit position in the current buffer element
    mutable uint64_t buf_item;

    mutable int buffer_index;
    mutable uint64_t read_mask;
    mutable int read_pos;
};

#endif //BIT_PACKER_H
