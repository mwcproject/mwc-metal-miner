# Copyright 2025 The MWC Developers
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.30)
project(mwc_cpu_miner)

find_package(jsoncpp REQUIRED)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_COMPILER /usr/bin/clang++)
set(CMAKE_C_COMPILER /usr/bin/clang)

if(NOT CMAKE_BUILD_TYPE)
    message(FATAL_ERROR "CMAKE_BUILD_TYPE is not set. Please specify a build type: Debug, Release).")
endif()

if(NOT BUILD_TARGET)
    set(BUILD_TARGET "CPU_MINER")
    #message(FATAL_ERROR "BUILD_TARGET is not set. Please specify a build target: CPU_MINER,TESTS).")
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
    # Enable ARMv9 optimizations (Apple M4)
    add_compile_options(-march=armv8.5-a+simd)
else()
    message(FATAL_ERROR "Cpu miner is optimized for ARMv9 (Apple M4). On other platfroms it will not run efficeitly.")
endif()

# Disable RTTI (Run-Time Type Information)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti -fno-exceptions  -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_DEBUG")
# Disable exceptions
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions")

# Compiler flags for Debug
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0")
# Compiler flags for Release
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

include_directories(../metal-cpp)

set(SOURCES
        src/miner/bit_packer.h
        src/miner/bit_streamer.h
        src/miner/blake.cpp
        src/miner/blake.h
        src/miner/bucket.h
        src/miner/cuckatoo.cpp
        src/miner/cuckatoo.h
        src/miner/data_storage.h
        src/miner/inplace_bucket_sorting.h
        src/miner/metal.cpp
        src/miner/metal.h
        src/miner/miner_network.cpp
        src/miner/miner_network.h
        src/miner/nonce_hash.h
        src/miner/row_builder.h
        src/miner/sip_hash.cpp
        src/miner/sip_hash.h
        src/miner/trim_nonces.cpp
        src/miner/trim_nonces.h
        src/miner/trimmer.h
        src/miner/utils.cpp
        src/miner/utils.h
        src/tests/miner_network_test.cpp
        src/tests/miner_network_test.h
)

if(BUILD_TARGET STREQUAL "CPU_MINER")
    list(APPEND SOURCES
            src/miner/main_cpu_miner.cpp
    )
elseif(BUILD_TARGET STREQUAL "TESTS")
    list(APPEND SOURCES
            src/tests/main_test.cpp
            src/tests/bit_packer_test.cpp
            src/tests/bit_packer_test.h
            src/tests/bit_streamer_test.cpp
            src/tests/bit_streamer_test.h
            src/tests/blake_test.cpp
            src/tests/blake_test.h
            src/tests/cuckatoo_test.cpp
            src/tests/cuckatoo_test.h
            src/tests/data_storage_test.cpp
            src/tests/data_storage_test.h
            src/tests/inplace_bucket_sorting_test.cpp
            src/tests/inplace_bucket_sorting_test.h
            src/tests/metal_test.cpp
            src/tests/metal_test.h
            src/tests/row_builder_test.cpp
            src/tests/row_builder_test.h
            src/tests/test_lean_miner.h
            src/tests/trimmer_test.cpp
            src/tests/trimmer_test.h
    )
else()
    message(FATAL_ERROR "Invalid BUILD_TARGET value. Expected CPU_MINER,TESTS.")
endif()

add_executable(mwc_cpu_miner ${SOURCES})

# Link Metal Frameworks and Metal-Cpp
target_link_libraries(mwc_cpu_miner PRIVATE
        jsoncpp_lib
        "-framework Metal"
        "-framework Foundation"
        "-framework QuartzCore"
        "-framework MetalKit"
        "-framework CoreGraphics"
        "-framework CoreFoundation"
        "-lobjc"
)