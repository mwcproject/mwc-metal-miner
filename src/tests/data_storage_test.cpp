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

#include "data_storage_test.h"
#include "../miner/data_storage.h"

template <uint N>
static void compare_data(uint64_t d1[N], uint64_t d2[N]) {
    for (int i = 0; i < N; i++) {
        assert(d1[i] == d2[i]);
    }
}

// It is easier to debug case
void test_extract_store_4b() {
    uint64_t data[4] = {0,0,0,0};

    store_hash_nonce<4,16>(data, 0, HashNonce(3,5) );
    compare_data<4>(data, (uint64_t[4]){0x00050003, 0,0 ,0});

    store_hash_nonce<4,16>(data, 2, HashNonce(7,1) );
    compare_data<4>(data, (uint64_t[4]){0x00050003, 0x00010007,0 ,0});

    store_hash_nonce<4,16>(data, 1, HashNonce(1,2) );
    compare_data<4>(data, (uint64_t[4]){0x0002000100050003, 0x00010007,0 ,0});

    store_hash_nonce<4,16>(data, 7, HashNonce(3,1) );
    compare_data<4>(data, (uint64_t[4]){0x0002000100050003, 0x00010007,0 ,0x0001000300000000});

    store_hash_nonce<4,16>(data, 2, HashNonce(3,2) );
    compare_data<4>(data, (uint64_t[4]){0x0002000100050003, 0x00020003,0 ,0x0001000300000000});

    store_hash_nonce<4,24>(data, 3, HashNonce(7,4) );
    compare_data<4>(data, (uint64_t[4]){0x0002000100050003, 0x0400000700020003,0 ,0x0001000300000000});

    HashNonce hn = extract_hash_nonce<4,16>(data, 0);
    assert(hn.hash == 3 && hn.nonce==5);
    hn = extract_hash_nonce<4,16>(data, 1);
    assert(hn.hash == 1 && hn.nonce==2);
    hn = extract_hash_nonce<4,16>(data, 2);
    assert(hn.hash == 3 && hn.nonce==2);
    hn = extract_hash_nonce<4,24>(data, 3);
    assert(hn.hash == 7 && hn.nonce==4);
    hn = extract_hash_nonce<4,16>(data, 4);
    assert(hn.is_zero());
    hn = extract_hash_nonce<4,16>(data, 5);
    assert(hn.is_zero());
    hn = extract_hash_nonce<4,16>(data, 6);
    assert(hn.is_zero());
    hn = extract_hash_nonce<4,16>(data, 7);
    assert(hn.hash == 3 && hn.nonce==1);
}

void test_extract_store_C31() {
    //sasdasda
    uint64_t data[4] = {0,0,0,0};

    store_hash_nonce<5,31>(data, 0, HashNonce(0x53212345, 0x21) );
    compare_data<4>(data, (uint64_t[4]){0x0010D3212345, 0,0 ,0});

    store_hash_nonce<5,31>(data, 1, HashNonce(0x12345678, 0x1) );
    compare_data<4>(data, (uint64_t[4]){0x34567810D3212345, 0x0092,0 ,0});

    store_hash_nonce<5,31>(data, 2, HashNonce(0x33221100, 0xff) );
    compare_data<4>(data, (uint64_t[4]){0x34567810D3212345, 0x007fB32211000092,0 ,0});

    store_hash_nonce<5,31>(data, 3, HashNonce(0x11221122, 0x11) );
    compare_data<4>(data, (uint64_t[4]){0x34567810D3212345, 0x227fB32211000092, 0x08912211 ,0});


    HashNonce hn = extract_hash_nonce<5,31>(data, 2);
    assert(hn.hash == 0x33221100 && hn.nonce==0xff);

    hn = extract_hash_nonce<5,31>(data, 4);
    assert(hn.is_zero());

    hn = extract_hash_nonce<5,31>(data, 0);
    assert(hn.hash == 0x53212345 && hn.nonce==0x21);

    hn = extract_hash_nonce<5,31>(data, 1);
    assert(hn.hash == 0x12345678 && hn.nonce==0x1);

}

void test_extract_store_C32() {
    //sasdasda
    uint64_t data[4] = {0,0,0,0};

    store_hash_nonce<5,32>(data, 0, HashNonce(0x53212345, 0x21) );
    compare_data<4>(data, (uint64_t[4]){0x002153212345, 0,0 ,0});

    store_hash_nonce<5,32>(data, 1, HashNonce(0x12345678, 0x1) );
    compare_data<4>(data, (uint64_t[4]){0x3456782153212345, 0x0112,0 ,0});

    store_hash_nonce<5,32>(data, 2, HashNonce(0x33221100, 0xff) );
    compare_data<4>(data, (uint64_t[4]){0x3456782153212345, 0xff332211000112,0 ,0});

    store_hash_nonce<5,32>(data, 3, HashNonce(0x11221122, 0x11) );
    compare_data<4>(data, (uint64_t[4]){0x3456782153212345, 0x22ff332211000112, 0x11112211 ,0});


    HashNonce hn = extract_hash_nonce<5,32>(data, 2);
    assert(hn.hash == 0x33221100 && hn.nonce==0xff);

    hn = extract_hash_nonce<5,32>(data, 4);
    assert(hn.is_zero());

    hn = extract_hash_nonce<5,32>(data, 0);
    assert(hn.hash == 0x53212345 && hn.nonce==0x21);

    hn = extract_hash_nonce<5,32>(data, 1);
    assert(hn.hash == 0x12345678 && hn.nonce==0x1);

}
