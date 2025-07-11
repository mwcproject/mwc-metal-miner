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


#include <arpa/inet.h>
#include <assert.h>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstring>
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <json/json.h>
#include <mutex>
#include <queue>
#include <random>
#include <string>
#include <sys/socket.h>
#include <thread>
#include <unistd.h>
#include <csignal>
#include <pthread.h>

#include "blake.h"
#include "cuckatoo.h"
#include "cuckatoo_pipeline.h"
#include "miner_network.h"
#include "utils.h"

#define VERSION "1.0.0"

// Miner

std::atomic<bool> exiting(false);

// Signal handler function
void signal_handler(int signal) {
    if (signal == SIGINT) {
        std::cout << "Exiting the miner, please wait..." << std::endl;
        exiting.store(true, std::memory_order_relaxed);
    }
}

int main(int argc, char* argv[]) {
    // Parse command-line arguments
    if (argc < 5) {
        std::cerr << "mwc-metal-miner Version: " << VERSION << std::endl;
        std::cerr << "Usage: " << argv[0] << " -node <host:port> -login <user_name> [-pass <password>]" << std::endl;
        return 1;
    }

    std::string nodeHost;
    int nodePort = -1;
    std::string loginName;
    std::string password;

    for (int i = 2; i < argc; i += 2) {
        std::string key = argv[i-1];
        std::string value = argv[i];

        if (key == "-node") {
            size_t colonPos = value.find(':');
            if (colonPos == std::string::npos) {
                std::cerr << "Invalid node address. Use <host:port> format." << std::endl;
                return 1;
            }
            nodeHost = value.substr(0, colonPos);
            nodePort = std::stoi(value.substr(colonPos + 1));
        } else if (key == "-login") {
            loginName = value;
        } else if (key == "-pass") {
            password = value;
        }
        else {
            std::cerr << "Unknown arguments: " << key << std::endl;
            return 1;
        }
    }

    if (nodePort == -1 || nodeHost.empty()) {
        std::cerr << "Please define -node host:port to connect to." << std::endl;
        return 1;
    }

    std::cout << "mwc-metal-miner Version: " << VERSION << std::endl;
    std::cout << "Connecting to node: " << nodeHost << ":" << nodePort << std::endl;
    std::cout << "Login: " << loginName << std::endl;
    std::cout << "Algorithm: C31" << std::endl;

    std::signal(SIGINT, signal_handler);

    const int edge_bits = 31;

    MinerNetwork network(edge_bits);

    std::random_device rd;
    std::mt19937_64 nonce_gen(rd()); // Fixed seed
    std::uniform_int_distribution<uint64_t> uint64_dist(0, UINT64_MAX);

    int iteration = 0;

    while(!exiting.load(std::memory_order_relaxed)) {
        std::cout << "Connecting to the node..." << std::endl;
        if (!network.connect(nodeHost, nodePort)) {
            std::cout << "Unable connect to the node. Waiting some time to reconnect." << std::endl;
            std::this_thread::sleep_for(std::chrono::seconds(30));
            continue;
        }

        // Start threads
        std::thread readerThread(&MinerNetwork::networkReaderThread, &network);
        std::thread writerThread(&MinerNetwork::networkWriterThread, &network);

        network.sendLoginMessage(loginName, password, true);

        std::chrono::time_point last_get_job_request_time = std::chrono::steady_clock::now();
        std::chrono::time_point last_keep_alive_request_time = last_get_job_request_time;

        CuckatooPipeline cuckatoo_pipeline(&network);

        while (network.is_running() && !exiting.load(std::memory_order_relaxed)) {
            auto now = std::chrono::steady_clock::now();
            if (now - last_keep_alive_request_time > std::chrono::seconds(60)) {
                network.sendKeepAliveRequest();
                last_keep_alive_request_time = now;
            }

            // Wait for a job
            CuckaJob currentTask = network.getActiveJob();

            if (!currentTask.is_valid()) {
                // waiting for the next job
                if (now - last_get_job_request_time > std::chrono::seconds(5)) {
                    std::cout << "Job pool is empty, requesting a new job from the node" << std::endl;
                    network.sendGetJobRequest();
                    last_get_job_request_time = now;
                }

                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                continue;
            }

            // Calculating nonces and V hashes
            assert(currentTask.is_valid());
            uint64_t nonce = uint64_dist(nonce_gen);

            auto sys_now = std::chrono::system_clock::now();
            std::time_t now_time_t = std::chrono::system_clock::to_time_t(sys_now);
            std::tm now_tm = *std::localtime(&now_time_t);

            std::cout << std::put_time(&now_tm, "%Y-%m-%d %H:%M:%S") <<  " Starting " << (++iteration) << " C31 job: " << currentTask.jobId << " for height: " << currentTask.height <<
                ", difficulty: " << currentTask.difficulty << ", nonce: " << std::hex << nonce << std::dec << std::endl;

            uint64_t v[4];
            currentTask.calculate_seed_hash(nonce, v);

            // Known solution
            //uint64_t v[4] = {0xcdc22e3228ae6ce4,0x9b702732c5917f3a,0x987dab75b205767,0x358d6c79d4355934};

            // Starting cucckatoo calculations...

            cuckatoo_pipeline.submit_task( currentTask.height, currentTask.jobId, nonce, v );

            //std::this_thread::sleep_for(std::chrono::seconds(1));
        }

        cuckatoo_pipeline.release();
        network.stop_running();

        // Wait for threads to finish
        readerThread.join();
        writerThread.join();
    }

    return 0;
}