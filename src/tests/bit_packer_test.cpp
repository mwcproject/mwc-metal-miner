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


#include "../miner/bit_packer.h"
#include "bit_packer_test.h"

void test_bit_packer_simple() {
    BitPacker bit_packer;

    assert(bit_packer.get_writer_bit_pos()==0);
    bit_packer.encode_number(0);
    assert(bit_packer.get_writer_bit_pos()==1);
    bit_packer.encode_number(1);
    assert(bit_packer.get_writer_bit_pos()==1+2);
    bit_packer.encode_number(5);
    assert(bit_packer.get_writer_bit_pos()==1+2+6);

    assert(bit_packer.buf_size() == 0);
    bit_packer.write_flush();
    assert(bit_packer.buf_size() == 1);

    bit_packer.reset_read();
    uint n1 = bit_packer.decode_number();
    assert(n1==0);
    uint n2 = bit_packer.decode_number();
    assert(n2==1);
    uint n3 = bit_packer.decode_number();
    assert(n3==5);

    bit_packer.reset_read();
    n1 = bit_packer.decode_number();
    assert(n1==0);

    bit_packer.set_bit_pos(1+2);
    n3 = bit_packer.decode_number();
    assert(n3==5);
}

void test_bit_packer_heavy() {
    BitPacker bit_packer;

    std::vector<uint> bit_pos;
    std::vector<uint> added_val;

    for (uint i=0; i<256; i++) {
        uint v = i % 16;
        bit_pos.push_back(bit_packer.get_writer_bit_pos());
        added_val.push_back(v);
        bit_packer.encode_number(v);
    }
    bit_packer.write_flush();
    assert(bit_packer.buf_size() <= bit_pos.back()/64 + 1);

    // Let's check if bit_pos match expected added_val
    for (int i=1; i<bit_pos.size(); i++ ) {
        int bits = bit_pos[i] - bit_pos[i-1];
        int val = added_val[i-1];
        assert(bits == val+1);
    }

    // Sequentual reading
    bit_packer.reset_read();
    for (uint i=0; i<256; i++) {
        uint v = i % 16;
        uint d = bit_packer.decode_number();
        assert(v==d);
    }

    bit_packer.set_bit_pos( bit_pos[10] );
    auto expected = added_val[10];
    uint d = bit_packer.decode_number();
    assert(d==expected);

    // Sequentual reading with set bit pos
    for (uint i=0; i<256; i++) {
        uint v = i % 16;
        bit_packer.set_bit_pos( bit_pos[i] );
        uint d = bit_packer.decode_number();
        assert(v==d);
    }

    // Random access reading
    for (uint i=0; i<256; i++) {
        uint idx = (i * 17) / bit_pos.size();
        bit_packer.set_bit_pos( bit_pos[idx] );
        uint d = bit_packer.decode_number();
        assert(d==added_val[idx]);
    }
}

