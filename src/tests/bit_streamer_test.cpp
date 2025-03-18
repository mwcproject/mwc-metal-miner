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

#include "../miner/bit_streamer.h"
#include "bit_streamer_test.h"


template<uint8_t BITS_NUM>
static void test_bit_streamer_impl2(const uint sz) {
    BitStreamer<BITS_NUM> streamer;

    for (uint i=0;i<sz; i++) {
        streamer.encode_number(1);
    }
    streamer.write_flush();
    assert(streamer.get_size() == sz);
    streamer.reset_read();
    for (int i=0;i<sz; i++) {
        uint n = streamer.decode_number();
        assert(n == 1);
    }
}

template<uint8_t BITS_NUM>
static void test_bit_streamer_impl() {
    for (int i=90; i<123;i++) {
        test_bit_streamer_impl2<BITS_NUM>(i);
    }
}

void test_bit_streamer() {
    test_bit_streamer_impl<7>();
    test_bit_streamer_impl<8>();
    test_bit_streamer_impl<9>();
    test_bit_streamer_impl<10>();
    test_bit_streamer_impl<11>();
    test_bit_streamer_impl<12>();
    test_bit_streamer_impl<13>();
    test_bit_streamer_impl<14>();
    test_bit_streamer_impl<15>();
    test_bit_streamer_impl<16>();
    test_bit_streamer_impl<17>();
    test_bit_streamer_impl<18>();
    test_bit_streamer_impl<19>();
    test_bit_streamer_impl<20>();
    test_bit_streamer_impl<21>();
    test_bit_streamer_impl<22>();
    test_bit_streamer_impl<23>();
    test_bit_streamer_impl<24>();
    test_bit_streamer_impl<25>();
    test_bit_streamer_impl<26>();
    test_bit_streamer_impl<27>();
    test_bit_streamer_impl<28>();
    test_bit_streamer_impl<29>();
    test_bit_streamer_impl<30>();
}

