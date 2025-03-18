//
// Created by Konstantin Bay on 3/16/25.
//

#include "utils.h"
#include  <assert.h>

static inline uint8_t char2u8(const char ch) {
    if (ch >= '0' && ch <= '9') {
        return ch - '0';
    }
    if (ch >= 'a' && ch <= 'f') {
        return 10 + (ch - 'a');
    }
    if (ch >= 'A' && ch <= 'F') {
        return 10 + (ch - 'A');
    }

    assert(false);
    return 0;
}

std::vector<uint8_t> hexstr2bin(const std::string& hex_string) {
    // Check if the string length is even
    std::vector<uint8_t> result;
    int len = hex_string.length();
    assert(len % 2 == 0); // invalid data, will be converted with error
    len -= len % 2; // just in case of currupted data
    const char* data = hex_string.c_str();
    result.reserve(len / 2); // 8 for the nonce (future usage)

    for (size_t i = 0; i < len; i += 2) {
        // Convert two hex characters to a uint8_t value
        const char ch1 = data[0];
        const char ch2 = data[1];
        data += 2;
        const uint8_t hi = char2u8(ch1);
        const uint8_t lo = char2u8(ch2);
        result.push_back( lo | (hi << 4) );
    }
    return result;
}

static const char * HEX_SYMBOLS = "0123456789ABCDEF";

std::string bin2hexstr(std::vector<uint8_t> & data) {
    std::string result;
    for (uint8_t d : data) {
        result.push_back(HEX_SYMBOLS[(d>>4) & 0xf]);
        result.push_back(HEX_SYMBOLS[d & 0xf]);
    }
    return result;
}

