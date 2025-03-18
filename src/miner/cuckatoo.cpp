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

#include "cuckatoo.h"

static bool is_prime(uint32_t num) {
    if (num < 2) return false;
    if (num % 2 == 0) return num == 2;
    for (uint32_t i = 3; i * i <= num; i += 2) {
        if (num % i == 0) return false;
    }
    return true;
}

// Function to get the next prime >= start
uint32_t next_prime(uint32_t start) {
    while (!is_prime(start)) {
        start++;
    }
    return start;
}
