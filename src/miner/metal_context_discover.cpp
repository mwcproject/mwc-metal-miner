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

#include "metal_context_discover.h"
#include "metal.h"
#include "cycle.h"

void MetalContextDiscover::init(MetalOps * metal, const SmallHashTable<bool> & keys, bool upass, uint64_t v[4]) {
    assert(metal);

    if (params==nullptr) {
        params = metal->device->newBuffer( sizeof(DiscoverParams), MTL::ResourceStorageModeShared );
    }

    const auto & keys_data = keys.get_keys();
    if (hh_keys==nullptr || hh_keys->length() < keys_data.size()*4) {
        RELEASE(hh_keys);
        hh_keys = metal->device->newBuffer( keys_data.size()*4, MTL::ResourceStorageModeShared );
    }

    if (nonces==nullptr || nonces->length() < keys_data.size()*4 ) {
        RELEASE(nonces);
        nonces = metal->device->newBuffer( keys_data.size()*4, MTL::ResourceStorageModeShared );
    }

    // Init params
    DiscoverParams * prms = (DiscoverParams *) params->contents();
    prms->key0 = v[0];
    prms->key1 = v[1];
    prms->key2 = v[2];
    prms->key3 = v[3];

    prms->nonce_pos = 0;
    prms->nonce_sz = keys_data.size();

    prms->hh_len = keys.get_hasher_sz();
    prms->hh_k = keys.get_hasher_mul();
    assert(keys.get_hasher_sz() == keys_data.size());
    prms->upass = upass;

    // init keys
    memcpy( hh_keys->contents(), keys_data.data(), keys_data.size()*4 );

    // nonces - no need to init, noise is fine
#ifndef NDEBUG
    memset( nonces->contents(), 0x69, keys_data.size()*4 );
#endif
}

void MetalContextDiscover::read_results(std::vector<NonceInfo> & discovered_nonces ) {
    DiscoverParams * prms = (DiscoverParams *) params->contents();
    discovered_nonces.resize( std::min(prms->nonce_pos, prms->nonce_sz) );
    memcpy( discovered_nonces.data(), nonces->contents(), discovered_nonces.size()*sizeof(NonceInfo) );
}


void MetalContextDiscover::release() {
    RELEASE(params);
    RELEASE(hh_keys);
    RELEASE(nonces);
}
