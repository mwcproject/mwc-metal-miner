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

/**
 * Cuckatoo Job from the mining node
 */
struct CuckaJob {
    // Job Id. Needed for the response.
    int jobId = -1;
    // Header pre POW.
    std::vector<uint8_t> prePow;
    // Difficulty threshold
    uint64_t difficulty = 0;
    // Blockchain height
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

    // calculate the hash foe this Job. Hash is calculated based on nonce and prePow
    void calculate_seed_hash( uint64_t nonce, uint64_t res_v[4]);
};

/**
 * Interface to submit the response back to the mining node
 */
class CuckatooResponser {
public:
    virtual void submit_response(int height, int jobId, uint64_t nonce, const std::vector<uint32_t> & res_nonces) = 0;
};

/**
 * Network connection to the mining node.
 */
struct MinerNetwork : public CuckatooResponser {
public:
    // Expected Edge bits for this miber. Expected only C31
    MinerNetwork(int edge_bits);
    ~MinerNetwork();

    // Connect to the moining node
    bool connect(const std::string & nodeHost, int nodePort);

    bool is_running() {
        return running.load();
    }
    void stop_running() {
        running.store(false);
        requestsCV.notify_all();
    }

    // Read/write threads
    void networkReaderThread();
    void networkWriterThread();

    // Stratum messages:
    // login
    void sendLoginMessage(const std::string & loginName, const std::string & password, bool sendGetJobRequest);
    // keep alive
    void sendKeepAliveRequest();
    // Get a new job request
    void sendGetJobRequest();

    // Submit solution back to miner node
    virtual void submit_response(int height, int jobId, uint64_t nonce, const std::vector<uint32_t> & res_nonces) override;

    // Get current active job that everybody is working on
    CuckaJob getActiveJob() const;

private:
    // send message to the node
    void sendJsonMessage(const Json::Value& jsonMessage);
    // read message fom the node
    std::string readJsonMessage();
private:
    int edge_bits; // for C31 must be 31
    std::atomic<bool> running{true};
    int sockfd = -1;

    // Request Id is requred by the stratum protocol
    std::atomic<int> request_id{1};

    // Job and response pools
    CuckaJob activeJob;
    mutable std::mutex jobMutex;
    // Requests messages weiting to be sent.
    std::queue<Json::Value> requestsPool;
    std::mutex requestsMutex;
    std::condition_variable requestsCV;

    // Persistent buffer for incomplete messages
    std::string remainingData;
};


#endif //MINER_NETWORK_H
