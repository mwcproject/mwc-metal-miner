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

#ifndef SMALL_HASH_TABLE_H
#define SMALL_HASH_TABLE_H

#include <vector>
#include <assert.h>
#include "utils.h"

const uint32_t PRIMES_NUM = 18;
const uint32_t PRIMES[PRIMES_NUM] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};

// It is a small hash table with predifined keys. Idea is to find the perfect has that will satisfy those int keys (no collisiopns)
// Because of birthday theorem, we limit the table size and don't include all the keys (accept our best)
template <typename VAL>
class SmallHashTable {
private:
    uint32_t hasher_sz = 0;
    uint32_t hasher_mul = 0;
    std::vector<uint32_t> keys;
    std::vector<VAL> values;
public:
    SmallHashTable() {
    }

    // return rejected keys
    std::vector<uint32_t> reset(const std::vector<uint32_t> & init_keys, uint32_t size_limit, const uint32_t key_empty, const VAL & init_value) {
        assert(!init_keys.empty());

        uint32_t table_prime = next_prime(std::max(32U, size_limit));
        uint32_t best_pr = PRIMES[0];
        uint32_t best_accepted = 0;
        for (uint32_t pi=0; pi<PRIMES_NUM; pi++) {
            uint32_t accepted = validate_hasher(table_prime, PRIMES[pi], init_keys);
            if (accepted>best_accepted) {
                best_accepted = accepted;
                best_pr = PRIMES[pi];
            }
        }

        hasher_sz = table_prime;
        hasher_mul = best_pr;
        keys.resize(0);
        keys.resize(hasher_sz, key_empty);
        values.resize(0);
        values.resize(hasher_sz, init_value);

        std::vector<uint32_t> rejected;

        for (auto k : init_keys) {
            if ( keys[hash(k)] != key_empty ) {
                rejected.push_back(k);
                continue;
            }
            keys[hash(k)] = k;
        }
        return rejected;
    }

    inline bool contains(uint32_t k) const {
        return keys[hash(k)] == k;
    }

    inline VAL& get(uint32_t k) {
        assert(contains(k));
        return values[hash(k)];
    }

    inline const VAL& get(uint32_t k) const {
        assert(contains(k));
        return values[hash(k)];
    }

    const std::vector<uint32_t> & get_keys() const { return keys; }

    uint32_t get_hasher_sz() const {return hasher_sz;}
    uint32_t get_hasher_mul() const {return hasher_mul;}
private:
    // return number of keys that fit without collision.
    uint32_t validate_hasher(uint32_t hasher, uint32_t mul, const std::vector<uint32_t> & ks) const {
        assert(hasher>ks.size());
        std::vector<bool> tracker(hasher, false);
        uint32_t accepted = 0;
        for (auto k : ks) {
            uint32_t idx = (k*mul) % hasher;
            if (tracker[idx])
                continue;
            tracker[idx] = true;
            accepted++;
        }
        return accepted;
    }

    inline uint32_t hash(uint32_t key) const {
        return (key*hasher_mul) % hasher_sz;
    }
};



#endif //SMALL_HASH_TABLE_H
