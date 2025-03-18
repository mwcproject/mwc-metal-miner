//
// Created by Konstantin Bay on 3/15/25.
//

#ifndef MINER_NETWORK_H
#define MINER_NETWORK_H

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

struct CuckaJob {
    int jobId = -1;
    std::vector<uint8_t> prePow;
    uint64_t difficulty = 0;
    int height = -1;

    void reset() {
        jobId = -1;
        prePow.clear();
        difficulty = 0;
        height = -1;
    }

    bool is_valid() {
        return jobId>=0 && !prePow.empty() && difficulty>0 && height>0;
    }

    void calculate_seed_hash( uint64_t nonce, uint64_t res_v[4]);
};

struct MinerNetwork {
public:
    MinerNetwork();
    ~MinerNetwork();

    bool connect(const std::string & nodeHost, int nodePort);

    bool is_running() {return running.load();}
    void stop_running() {running.store(false);}

    void networkReaderThread();
    void networkWriterThread();

    void sendLoginMessage(const std::string & loginName, const std::string & password, bool sendGetJobRequest);
    void sendKeepAliveRequest();
    void sendGetJobRequest();
    void sendResponseRequest(int edge_bits, const CuckaJob & job, uint64_t nonce, const std::vector<uint64_t> & res_nonces);

    CuckaJob getActiveJob() const;

private:
    void sendJsonMessage(const Json::Value& jsonMessage);
    std::string readJsonMessage();
private:
    std::atomic<bool> running{true};
    int sockfd = -1;

    int request_id = 1;

    // Job and response pools
    CuckaJob activeJob;
    mutable std::mutex jobMutex;
    std::queue<Json::Value> requestsPool;
    std::mutex requestsMutex;
    std::condition_variable requestsCV;

    // Persistent buffer for incomplete messages
    std::string remainingData;
};


#endif //MINER_NETWORK_H
