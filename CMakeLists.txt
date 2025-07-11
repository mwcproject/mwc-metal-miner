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
project(mwc-metal-miner)

find_package(jsoncpp REQUIRED)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_COMPILER /usr/bin/clang++)
set(CMAKE_C_COMPILER /usr/bin/clang)

if(NOT CMAKE_BUILD_TYPE)
    message(FATAL_ERROR "CMAKE_BUILD_TYPE is not set. Please specify a build type: Debug, Release).")
endif()

if(NOT BUILD_TARGET)
    #set(BUILD_TARGET "PERF_TEST")
    set(BUILD_TARGET "CPU_MINER")
    #message(FATAL_ERROR "BUILD_TARGET is not set. Please specify a build target: CPU_MINER, PERF_TEST).")
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
    # Enable ARMv9 optimizations (Apple M4)
    add_compile_options(-march=armv8.5-a+simd)
else()
    message(FATAL_ERROR "Cpu miner is optimized for ARMv9 (Apple M4). On other platforms it will not run well.")
endif()

# Disable RTTI(Run-Time Type Information)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-rtti -fno-exceptions  ")
# Disable exceptions
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-exceptions")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_DEBUG")

# Compiler flags for Debug
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -O0")
# Compiler flags for Release
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG")

# Add debug info to release
# set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g")

include_directories(../metal-cpp)

set(SOURCES
        src/miner/blake.cpp
        src/miner/blake.h
        src/miner/cuckatoo.cpp
        src/miner/cuckatoo.h
        src/miner/cuckatoo_pipeline.cpp
        src/miner/cuckatoo_pipeline.h
        src/miner/features.h
        src/miner/metal_structs.h
        src/miner/metal.cpp
        src/miner/metal.h
        src/miner/metal_code.cpp
        src/miner/metal_code.h
        src/miner/metal_context.cpp
        src/miner/metal_context.h
        src/miner/metal_context_phase1.cpp
        src/miner/metal_context_phase1.h
        src/miner/metal_context_phase2.cpp
        src/miner/metal_context_phase2.h
        src/miner/mem_pool.cpp
        src/miner/mem_pool.h
        src/miner/miner_network.cpp
        src/miner/miner_network.h
        src/miner/sip_hash.cpp
        src/miner/sip_hash.h
        src/miner/utils.cpp
        src/miner/utils.h
        src/miner/events_tracker.h
        src/miner/thread_pool.h
        src/miner/cycle_finder.h
        src/miner/cycle_finder.cpp
        src/miner/small_hash_table.h
        src/miner/cycle.h
        src/miner/cycle.cpp
        src/miner/metal_context_discover.h
        src/miner/metal_context_discover.cpp
        src/tests/test_st1_build_buckets.cpp
        src/tests/test_st1_build_buckets.h
        src/tests/test_mem_pool.cpp
        src/tests/test_mem_pool.h
        src/tests/test_nonce_to_8b.cpp
        src/tests/test_nonce_to_8b.h
        src/tests/test_build_mask.h
        src/tests/test_build_mask.cpp
)

if(BUILD_TARGET STREQUAL "CPU_MINER")
    list(APPEND SOURCES
            src/miner/main_cpu_miner.cpp
    )
elseif(BUILD_TARGET STREQUAL "PERF_TEST")

    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DTRACK_EVENTS")

    list(APPEND SOURCES
            src/perf/events_tracker.cpp
            src/perf/main_perf.cpp
    )
else()
    message(FATAL_ERROR "Invalid BUILD_TARGET value. Expected CPU_MINER,PERF_TEST.")
endif()

add_executable(mwc-metal-miner ${SOURCES})

# Link Metal Frameworks and Metal-Cpp
target_link_libraries(mwc-metal-miner PRIVATE
        jsoncpp_lib
        "-framework Metal"
        "-framework Foundation"
        "-framework QuartzCore"
        "-framework MetalKit"
        "-framework CoreGraphics"
        "-framework CoreFoundation"
        "-lobjc"
)