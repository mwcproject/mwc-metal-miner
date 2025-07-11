#include "test_nonce_to_8b.h"
#include <map>
#include <assert.h>
#include "../miner/sip_hash.h"

void test_nonce_to_8b_check(const uint64_t * v,
            MetalContext & context,
            const MemRange & in,
            const MemRange & out,
            uint EDGE_BITS)
{
    const uint hash_mask = (uint(1) << EDGE_BITS) - 1;

    uint32_t * data_in = (uint32_t *) in.get_data_ptr(context);
    uint64_t * data_out = (uint64_t *) out.get_data_ptr(context);

    uint32_t sz = in.get_length_bytes()/4;
    assert(sz == out.get_length_bytes()/8);

    for (uint i = 0; i < sz; i++) {
        uint32_t nonce = data_in[i];
        uint64_t hashes = data_out[i];

        if (nonce == 0xFFFFFFFF) {
            assert(hashes == 0);
        }
        else {
            uint32_t hash11 = sip_hash(v, nonce*2) & hash_mask;

            uint32_t hash21 = uint32_t(hashes);
            uint32_t nonce22 = uint32_t(hashes>>32);

            assert(hash11 == hash21);
            assert(nonce22 == nonce);
        }
    }

}
