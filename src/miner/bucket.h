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

#ifndef MASKED_BUCKET_H
#define MASKED_BUCKET_H

#include "bit_packer.h"
#include "sip_hash.h"

struct IntervalsIndex {
    uint32_t bitN = 0;
    uint32_t intervalsPos = 0;
    uint32_t nonceIdx = 0;

    inline void reset() {
        bitN = intervalsPos = nonceIdx = 0;
    }

    inline void set(uint32_t _bitN, uint32_t _intervalsPos, uint32_t _nonceIdx) {
        bitN = _bitN;
        intervalsPos = _intervalsPos;
        nonceIdx = _nonceIdx;
    }
};

// Here we trying to solve problem of storage multiple integers related to the mask.
// We could have another mask, but it has larger footprint. Alternative is a BitPacker
// that can provide siquentual access.
// This bucket size is expected to be around 4096 for C31 and 8192 for C32. or tests in can be much larger
// interval_indexes are pointed at bits at every around
template <uint8_t EDGE_BITS, uint8_t BUCKET_BITS, uint8_t MASK_BITS, uint32_t BIT_PACKER_INDEX_INTERVAL>
struct Bucket {
    Bucket() = default;
    Bucket( const Bucket &) = delete;
    Bucket & operator = ( const Bucket & ) = delete;

    // Need to build, but still never should be called
    Bucket( const Bucket && other) {
        assert(false);
    }

    BitPacker intervals; // hash related intervals for nonces

    // get<0> - bit number at intervals, get<1> - index at interval Indexes
    IntervalsIndex interval_indexes[ (1<<(EDGE_BITS-BUCKET_BITS*2))/BIT_PACKER_INDEX_INTERVAL-1]; // it is indexes for the nonces at offsets 128, 128*2, 128*3 ...
    std::vector<uint32_t> nonces;

    // idx - index inside this bucket
    uint get_nonces(uint idx, IntervalsIndex & index, uint32_t resulting_nonces[1<<MASK_BITS] ) const {
        uint32_t ii = idx / BIT_PACKER_INDEX_INTERVAL;
        size_t interval_indexes_length = sizeof(interval_indexes) / sizeof(interval_indexes[0]);
        assert(ii<interval_indexes_length+1);

        if ( index.bitN > idx || (ii>0 && index.bitN < interval_indexes[ii-1].bitN)  ) {
            if (ii>0)
                index = interval_indexes[ii-1];
            else
                index.reset();
        }

        intervals.set_bit_pos(index.intervalsPos);
        while ( index.bitN < idx ) {
            index.bitN++;
            index.nonceIdx += intervals.decode_number();
        }
        int num = intervals.decode_number();
        assert(index.nonceIdx + num <= nonces.size());

        // read the data from nonces into resulting_nonces
        for ( uint i = 0; i < num; i++ ) {
            resulting_nonces[i] = nonces[index.nonceIdx + i];
        }

        index.nonceIdx += num;
        index.bitN++;
        index.intervalsPos = intervals.get_reader_bit_pos();

        return num;
    }

    // idx - index inside this bucket
    uint get_nonces(uint idx, uint32_t resulting_nonces[1<<MASK_BITS] ) const {
        uint32_t ii = idx / BIT_PACKER_INDEX_INTERVAL;
        size_t interval_indexes_length = sizeof(interval_indexes) / sizeof(interval_indexes[0]);
        assert(ii<interval_indexes_length+1);

        IntervalsIndex index = ii>0 ? interval_indexes[ii-1] : IntervalsIndex();

        intervals.set_bit_pos(index.intervalsPos);
        while ( index.bitN < idx ) {
            index.bitN++;
            index.nonceIdx += intervals.decode_number();
        }
        int num = intervals.decode_number();
        assert(index.nonceIdx + num <= nonces.size());

        // read the data from nonces into resulting_nonces
        for ( uint i = 0; i < num; i++ ) {
            resulting_nonces[i] = nonces[index.nonceIdx + i];
        }

        index.nonceIdx += num;
        index.bitN++;
        index.intervalsPos = intervals.get_reader_bit_pos();

        return num;
    }

};



#endif //MASKED_BUCKET_H
