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

#ifndef TRIM_NONCES_H
#define TRIM_NONCES_H

#include "metal.h"

// Simple class that store and manage trimmed nonces data into the Metal buffers.
// Since we are using it in two routines, lets encapsulate it here
class TrimNonces {
private:
    MetalOps *metal;
    uint32_t max_nonces;
    std::vector<MTL::Buffer *> trimmed_nonces;
    uint32_t * cur_buf;
    uint32_t cur_idx;
    uint32_t total_added_nonces = 0;
public:
    TrimNonces(MetalOps * _metal, uint32_t _max_nonces) : metal(_metal), max_nonces(_max_nonces) {
        assert(max_nonces >= 64);
        MTL::Buffer * buf = metal->create_buffer(max_nonces * sizeof(uint32_t));
        trimmed_nonces.push_back(buf);
        cur_buf = (uint32_t *) buf->contents();
        cur_idx = 0;
    }

    // Init call is expected:
    TrimNonces() {}

    ~TrimNonces() {
        assert(trimmed_nonces.empty()); // expected that caller takes case about those buffers
    }

    void init(MetalOps * _metal, uint32_t _max_nonces) {
        metal = _metal;
        max_nonces = _max_nonces;
        assert(max_nonces >= 64);
        MTL::Buffer * buf = metal->create_buffer(max_nonces * sizeof(uint32_t));
        trimmed_nonces.push_back(buf);
        cur_buf = (uint32_t *) buf->contents();
        cur_idx = 0;
    }

    inline void add_nonce(uint32_t nonce) {
        if (cur_idx>=max_nonces) {
            MTL::Buffer * buf = metal->create_buffer(max_nonces * sizeof(uint32_t));
            trimmed_nonces.push_back(buf);
            cur_buf = (uint32_t *) buf->contents();
            cur_idx = 0;
        }
        cur_buf[cur_idx++] = nonce;
        total_added_nonces++;
    }


    // get buffer to write into.
    inline uint32_t * get_nonce_buf(uint32_t reserved_cpacity) {
        if (cur_idx+reserved_cpacity<=max_nonces) {
            return &cur_buf[cur_idx];
        }
        else {
            for (; cur_idx<max_nonces; cur_idx++) {
                cur_buf[cur_idx] = UINT_MAX;
            }
            MTL::Buffer * buf = metal->create_buffer(max_nonces * sizeof(uint32_t));
            trimmed_nonces.push_back(buf);
            cur_buf = (uint32_t *) buf->contents();
            cur_idx = 0;
            return cur_buf;
        }
    }


    inline void update_written_num(uint32_t num) {
        cur_idx+=num;
        assert(cur_idx <= max_nonces);
        total_added_nonces+=num;
    }

    void finish();

    std::vector<MTL::Buffer *> & get_trimmed_nonces() {return trimmed_nonces;}

    uint32_t get_total_added_nonces() const {return total_added_nonces;}
};



#endif //TRIM_NONCES_H
