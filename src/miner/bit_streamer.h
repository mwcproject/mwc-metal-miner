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

#ifndef BITSTREAMER_H
#define BITSTREAMER_H

#include <vector>
#include <cstdint>
#include <assert.h>
#include <iostream>
#include <ostream>
#include <iomanip>

// We need to be able to write/read the numbers below 32 bits as a stream
template <uint8_t BITS_NUM>
class BitStreamer {
public:
    BitStreamer() = default;
    BitStreamer(const BitStreamer & ) = delete;
    BitStreamer(BitStreamer && ) = default;
    BitStreamer & operator = (const BitStreamer & ) = delete;

    // Add a number to the encoder
    inline void encode_number(uint32_t num) {
        assert(uint32_t(num) < (uint32_t(1) << BITS_NUM));
        assert(bit_position<64);

        // Calculate how many bits can be written to the current buffer element
        uint8_t bitsAvailable = 64 - bit_position;
        if (BITS_NUM < bitsAvailable) {
            buf_item |= (uint64_t(num) << bit_position);
            bit_position += BITS_NUM;
            assert(bit_position < 64);
        }
        else {
            // Write the bits to the current buffer element
            buf_item |= (uint64_t(num) << bit_position);

            // Update the value and bit position
            num >>= bitsAvailable;
            assert(bit_position+bitsAvailable == 64);

            assert(bitsAvailable<=BITS_NUM);
            buffer.push_back(buf_item);
            buf_item = num;
            bit_position = BITS_NUM-bitsAvailable;
        }
        size++;
    }

    void write_flush() {
        if (size>0) {
            buffer.push_back(buf_item);
            buf_item = 0;
            bit_position = 0;
        }
    }

    inline uint32_t get_size() {return size;}

    void reset_read() {
        bit_position = 0;
        assert(!buffer.empty());
        buf_item = buffer[0];
        read_buf_idx = 1;
    }

    inline uint32_t decode_number() {
        uint32_t result = 0;
        assert(bit_position<64);

        // Calculate how many bits can be written to the current buffer element
        uint8_t bitsAvailable = 64 - bit_position;
        if (BITS_NUM < bitsAvailable) {
            result = (buf_item >> bit_position);
            bit_position += BITS_NUM;
            assert(bit_position < 64);
        }
        else {
            result = (buf_item >> bit_position);

            // Update the value and bit position
            assert(bit_position+bitsAvailable == 64);

            assert(bitsAvailable<=BITS_NUM);
            buf_item = buffer[read_buf_idx++];
            bit_position = BITS_NUM-bitsAvailable;
            result = result  | ((buf_item & ((1<<bit_position)-1) )  << bitsAvailable);
        }

        return  result & ((1 << BITS_NUM)-1);
    }
private:
    std::vector<uint64_t> buffer;
    int bit_position = 0;              // Current bit position in the current buffer element
    uint64_t buf_item = 0;
    uint32_t read_buf_idx = 0;

    uint32_t size = 0;

};

// Note, dump data is destructive for write
template<uint8_t BITS_NUM>
void dump_bit_steam_data( BitStreamer<BITS_NUM> & bs ) {
    uint sz = bs.get_size();

    bs.reset_read();

    std::cout << std::hex << std::uppercase << std::setw((BITS_NUM+3)/4) << std::setfill('0');
    for( uint i = 0; i < sz; i++ ) {
        uint32_t n = bs.decode_number();
        std::cout << n << "  ";
        if ( i % 16 == 15 )
            std::cout << std::endl;
    }
    std::cout << std::endl  << std::endl;
}

#endif //BITSTREAMER_H
