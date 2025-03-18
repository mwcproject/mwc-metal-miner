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

#include "sip_hash.h"

#define rotl(x,b) ((x << b) | ( x >> (64 - b)))

#define SIPROUND \
v0 += v1; v2 += v3; v1 = rotl(v1,13); 	\
v3 = rotl(v3,16); v1 ^= v0; v3 ^= v2;	\
v0 = rotl(v0,32); v2 += v1; v0 += v3;	\
v1 = rotl(v1,17);   v3 = rotl(v3,21);	\
v1 ^= v2; v3 ^= v0; v2 = rotl(v2,32);


uint32_t sip_hash(const uint64_t * v, uint64_t nonce) {
    uint64_t v0 = v[0];
    uint64_t v1 = v[1];
    uint64_t v2 = v[2];
    uint64_t v3 = v[3] ^ nonce;

    SIPROUND; SIPROUND;
    v0 ^= nonce;
    v2 ^= 0xff;
    SIPROUND; SIPROUND; SIPROUND; SIPROUND;
    return uint32_t(v0 ^ v1 ^ v2  ^ v3);
}

