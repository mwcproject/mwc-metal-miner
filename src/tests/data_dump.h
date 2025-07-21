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

#ifndef DATA_DUMP_H
#define DATA_DUMP_H

#include <Metal/Metal.hpp>
#include <fstream>
#include <iostream>
#include "../miner/features.h"

template<typename T>
void write_dump( const std::vector<T> & data, const std::string & file_name ) {
    // Open file for binary output with truncation (overwrite if it exists)
    std::ofstream out( std::string(DUMP_DATA_DIR) + file_name, std::ios::binary | std::ios::trunc);

    if (!out) {
        assert(false);
        std::cout << "Unable to open the output dump file " << file_name << std::endl;
        exit(1);
    }

    // Write the contents of the vector to the file as raw binary
    out.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(T));
}

template<typename T>
void read_dump( const std::string & file_name, std::vector<T> & data ) {
    // Open file for binary output with truncation (overwrite if it exists)
    const std::string filename = std::string(DUMP_DATA_DIR) + file_name;
    std::ifstream file(filename, std::ios::binary);

    // Check if the file is open
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        assert(false);
        return;
    }

    // Determine the size of the file
    file.seekg(0, std::ios::end);
    std::streamsize fileSize = file.tellg();
    file.seekg(0, std::ios::beg);

    assert( fileSize % sizeof(T) == 0 );
    data.resize( fileSize / sizeof(T) );

    // Read the file into the buffer
    file.read((char*) data.data(), data.size() * sizeof(T));

    // Check if the read operation was successful
    if (!file) {
        std::cerr << "Error reading file: " << filename << std::endl;
        assert(false);
        return;
    }

    // Close the file
    file.close();
}


#endif //DATA_DUMP_H
