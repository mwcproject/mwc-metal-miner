//
// Created by Konstantin Bay on 3/16/25.
//

#ifndef UTILS_H
#define UTILS_H

#include <inttypes.h>
#include <vector>

std::vector<uint8_t> hexstr2bin(const std::string& hex_string);

std::string bin2hexstr(std::vector<uint8_t> & result);

// Calculating next prime NOT efficient way. Ise this method for a debug or a single time calculations during init process
uint32_t next_prime(uint32_t start);


#endif //UTILS_H
