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

#include "../miner/blake.h"
#include "blake_test.h"
#include "../miner/cuckatoo.h"

void test_blake2() {
    // basic test
    uint8_t h1[32]; // resulting hash
    blake2b(h1, (uint8_t*)"Hello, BLAKE2b!", strlen("Hello, BLAKE2b!"));  //  (uint8_t *) nonces.data(), nonces.size()*sizeof(uint64_t));
    // Check with: https://toolkitbay.com/tkb/tool/BLAKE2b_256
    // Output: 5d5415132e08e43caed7cd44413114c4167fd6719656c4cd0003eab79a96adb8
    assert(h1[0] == 0x5d);
    assert(h1[1] == 0x54);
    assert(h1[2] == 0x15);


    // Data comes from mwc_node test   fn test_proof_rw()
    std::vector<uint64_t> nonces = {0x1128e07, 0xc181131, 0x110fad36, 0x1135ddee, 0x1669c7d3, 0x1931e6ea, 0x1c0005f3,
        0x1dd6ecca, 0x1e29ce7e, 0x209736fc, 0x2692bf1a, 0x27b85aa9, 0x29bb7693, 0x2dc2a047,
        0x2e28650a, 0x2f381195, 0x350eb3f9, 0x3beed728, 0x3e861cbc, 0x41448cc1, 0x41f08f6d,
        0x42fbc48a, 0x4383ab31, 0x4389c61f, 0x4540a5ce, 0x49a17405, 0x50372ded, 0x512f0db0,
        0x588b6288, 0x5a36aa46, 0x5c29e1fe, 0x6118ab16, 0x634705b5, 0x6633d190, 0x6683782f,
        0x6728b6e1, 0x67adfb45, 0x68ae2306, 0x6d60f5e1, 0x78af3c4f, 0x7dde51ab, 0x7faced21};
    assert(nonces.size() == 42);
    // must be 128 byte aligned
    //nonces.resize( (nonces.size() + 15)/16*16, 0ULL );
    //assert(nonces.size() == 48);

    std::vector<uint8_t> nonces_compressed = compress_nonces<31>(nonces);
    //nonces_compressed.resize(129);

    uint8_t h[32]; // resulting hash
    blake2b(h, nonces_compressed.data(), nonces_compressed.size());
    assert(h[0] == 0x5c);
    assert(h[1] == 0x29);
    assert(h[2] == 0x17);
    assert(h[3] == 0x4b);
    assert(h[4] == 0xce);
    assert(h[5] == 0x6e);
    assert(h[6] == 0xa6);
    assert(h[7] == 0x5a);
    assert(h[8] == 0xe9);
    assert(h[9] == 0xf8);
    assert(h[10] == 0x87);
    assert(h[11] == 0xe6);
    assert(h[12] == 0x1e);
    assert(h[13] == 0xdf);
    assert(h[14] == 0x0c);
    assert(h[15] == 0xad);
    assert(h[16] == 0x38);
    assert(h[17] == 0x02);
    assert(h[18] == 0x6b);
    assert(h[19] == 0x62);
    assert(h[20] == 0x1f);
    assert(h[21] == 0xf4);
    assert(h[22] == 0x78);
    assert(h[23] == 0x7e);
    assert(h[24] == 0xc9);
    assert(h[25] == 0x8d);
    assert(h[26] == 0xa2);
    assert(h[27] == 0x26);
    assert(h[28] == 0x8d);
    assert(h[29] == 0x88);
    assert(h[30] == 0xe6);
    assert(h[31] == 0x24);
}
