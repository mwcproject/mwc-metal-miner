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


#include <iostream>
#include <sstream>
#include <cstdint>
#include <chrono>

#include "bit_packer_test.h"
#include "data_storage_test.h"
#include "inplace_bucket_sorting_test.h"
#include "metal_test.h"
#include "row_builder_test.h"
#include "trimmer_test.h"
#include "cuckatoo_test.h"
#include "blake_test.h"
#include "bit_streamer_test.h"
#include "miner_network_test.h"

int main(int argc, char* argv[]) {

    test_bucket_read_write();

    test_bit_packer_simple();
    test_bit_packer_heavy();

    test_extract_store_4b();
    test_extract_store_C31();
    test_extract_store_C32();

    test_inplace_sorting1();
    test_inplace_sorting2();

    test_metal_4b();
    test_metal_5b();

    test_row_builder_4b();
    test_row_builder_5b();

    test_row_builder_C25();

    test_row_builder_C27();
    test_row_builder_C28();

    test_trim_step1();

    test_nonces_to_hash();

    test_bit_streamer();

    test_cuckatoo_12_4b();
    test_cuckatoo_12_5b();

    test_cuckatoo_25();
    test_cuckatoo_29();

    test_blake2();

    test_cuckatoo_15_with_solution();

    test_utils();

    test_pow_nonce_hash();

    std::cout << "Tests: DONE" << std::endl;

    return 1;
}