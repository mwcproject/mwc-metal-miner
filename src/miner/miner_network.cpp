//
// Created by Konstantin Bay on 3/15/25.
//

#include "miner_network.h"
#include "blake.h"
#include "utils.h"
#include "features.h"

const std::string AGENT = "mwc-metal-miner";

void CuckaJob::calculate_seed_hash( uint64_t nonce, uint64_t res_v[4]) {
    for (int i = 0; i < 8; i++) {
        prePow.push_back( nonce >> ((7-i)*8) & 0xff );
    }
    blake2b((uint8_t*) res_v, prePow.data(), prePow.size() );
    // revert back inserted nonce
    prePow.resize(prePow.size() - sizeof(uint64_t));
}


////////////////////////////////////////////////////////////////////////////////////////////
// MinerNetwork

MinerNetwork::MinerNetwork(int _edge_bits) : edge_bits(_edge_bits) {
}

MinerNetwork::~MinerNetwork() {
    if (sockfd>=0) {
        close(sockfd);
    }
}


bool MinerNetwork::connect(const std::string & nodeHost, int nodePort) {
    if (sockfd>=0) {
        close(sockfd);
        sockfd = -1;
    }

    // Create TCP socket
    sockfd = socket(AF_INET, SOCK_STREAM, 0);
    if (sockfd < 0) {
        std::cerr << "Failed to create socket!" << std::endl;
        return false;
    }

    // Set 5 seconds connect timeout
    struct timeval timeout;
    timeout.tv_sec = 5;  // 5 seconds timeout
    timeout.tv_usec = 0; // No microseconds

    if (setsockopt(sockfd, SOL_SOCKET, SO_SNDTIMEO, &timeout, sizeof(timeout)) < 0) {
        std::cerr << "Failed to set timeout settings to socket. Error: " << errno << std::endl;
        return false;
    }

    // Connect to the server
    struct sockaddr_in serverAddr;
    serverAddr.sin_family = AF_INET;
    serverAddr.sin_port = htons(nodePort);
    inet_pton(AF_INET, nodeHost.c_str(), &serverAddr.sin_addr);

    if (::connect(sockfd, (struct sockaddr*)&serverAddr, sizeof(serverAddr)) < 0) {
        std::cerr << "Failed to connect to the server! Error: " << errno << std::endl;
        close(sockfd);
        sockfd = -1;
        return false;
    }

    // Set socket to non-blocking mode
    //fcntl(sockfd, F_SETFL, O_NONBLOCK);

    return sockfd >=0;
}


std::string json2str(const Json::Value & val) {
    if (val.isString())
        return val.asString();

    Json::StreamWriterBuilder writer;
    std::string jsonString = Json::writeString(writer, val);
    return jsonString;
}



// Send login message first
static Json::Value generateLoginRequest(int request_id, const std::string & loginName, const std::string & password) {
    Json::Value request;
    request["id"] = request_id;
    request["jsonrpc"] = "2.0";
    request["method"] = "login";
    Json::Value & params = request["params"];

    params["login"] = loginName;
    params["pass"] = password;
    params["agent"] = AGENT;
    return request;
}

// Get job template request
static Json::Value generateGetJobRequest(int request_id) {
    Json::Value request;
    request["id"] = request_id;
    request["jsonrpc"] = "2.0";
    request["method"] = "getjobtemplate";
    request["params"] =  Json::nullValue;
    return request;
}

// Get keepalive request
static Json::Value generateKeepAliveRequest(int request_id) {
    Json::Value request;
    request["id"] = request_id;
    request["jsonrpc"] = "2.0";
    request["method"] = "keepalive";
    request["params"] =  Json::nullValue;
    return request;
}

static Json::Value generateSubmitRequest(int request_id, int edge_bits, int height, int jobId,
        uint64_t nonce, const std::vector<uint32_t> & res_nonces )
{
    Json::Value request;
    request["id"] = request_id;
    request["jsonrpc"] = "2.0";
    request["method"] = "submit";
    Json::Value & params = request["params"];
    params["edge_bits"] = edge_bits;
    params["height"] = height;
    params["job_id"] = jobId;
    params["nonce"] = nonce;
    Json::Value pow(Json::arrayValue); // Create a JSON array
    for (uint64_t nonce : res_nonces) {
        pow.append(uint64_t(nonce));
    }
    params["pow"] = pow;
    return request;
}


// Function to send a JSON message over the TCP socket
void MinerNetwork::sendJsonMessage(const Json::Value& jsonMessage) {
    Json::StreamWriterBuilder writer;
    writer["indentation"] = "";
    std::string message = Json::writeString(writer, jsonMessage) + "\n";
    send(sockfd, message.c_str(), message.size(), 0);
}

// Function to read a single JSON message from the TCP socket
std::string MinerNetwork::readJsonMessage() {
    char buffer[1024];
    std::string message;

    while (true) {
        // Check if we already have a complete message in remainingData
        size_t newlinePos = remainingData.find('\n');
        if (newlinePos != std::string::npos) {
            // Extract the complete message
            message = remainingData.substr(0, newlinePos);
            // Remove the processed message from remainingData
            remainingData.erase(0, newlinePos + 1);
            return message;
        }

        // Read more data from the socket
        ssize_t bytesRead = recv(sockfd, buffer, sizeof(buffer) - 1, 0);
        if (bytesRead <= 0) {
            if (bytesRead == 0) {
                std::cerr << "Server closed the connection." << std::endl;
            } else if (errno != EAGAIN && errno != EWOULDBLOCK) {
                std::cerr << "Failed to read from socket: " << strerror(errno) << std::endl;
            }
            //int err = errno;
            running = false;
            requestsCV.notify_one(); // triggering exit for Writer thread
            return "";
        }

        // Append the new data to remainingData
        buffer[bytesRead] = '\0';
        remainingData += buffer;
    }
}

// Function to parse JSON response
Json::Value parseJson(const std::string& jsonString) {
    Json::Value root;
    Json::CharReaderBuilder reader;
    std::string errs;

    std::istringstream s(jsonString);
    if (!Json::parseFromStream(reader, s, &root, &errs)) {
        std::cerr << "Failed to parse JSON: " << jsonString << "  Error:" << errs << std::endl;
    }

    return root;
}

// Network reader thread
void MinerNetwork::networkReaderThread() {
    pthread_setname_np("networkReaderThread");
    while (running.load(std::memory_order_relaxed)) {
        std::string message = readJsonMessage();
        if (!message.empty()) {
            Json::Value jsonMessage = parseJson(message);

            if (!jsonMessage.isMember("method")) {
                std::cout << "Got unexpected message from the node: " << message << std::endl;
            }
            else {
                std::string method = jsonMessage["method"].asString();
                if (method == "job" || method == "getjobtemplate" ) {
                    // calculating prepow here because we don't want do free main solver thread

                    Json::Value jobParams;
                    if (jsonMessage.isMember("params")) {
                        jobParams = jsonMessage["params"];
                    }
                    else if (jsonMessage.isMember("result")) {
                        jobParams = jsonMessage["result"];
                    }

                    if (jobParams.isNull() || !jobParams.isMember("job_id")  || !jobParams.isMember("pre_pow")
                            || !jobParams.isMember("difficulty")  || !jobParams.isMember("height") ) {

                        std::cout << "Get invalid job params from the node: " << message << std::endl;
                        continue;
                    }

                    CuckaJob task;
                    task.jobId = jobParams["job_id"].asInt();
                    task.prePow = hexstr2bin( jobParams["pre_pow"].asString() );
                    task.difficulty = jobParams["difficulty"].asUInt64();
                    task.height = jobParams["height"].asInt();

                    // Evaluating prePOW
                    {
                        std::lock_guard<std::mutex> lock(jobMutex);
                        activeJob = task;
                    }
                }
                else {
                    if (jsonMessage.isMember("result")) {
#ifdef SHOW_NETWORK
                        std::cout << "Response for " << method << ": " << json2str(jsonMessage["result"]) << std::endl;
#endif
                    }
                    if (jsonMessage.isMember("error")) {
                        auto err_val = jsonMessage["error"];
                        if (!err_val.isNull()) {
                            std::cout << "Error response for " << method << ": " << json2str(err_val) << std::endl;
                            if (err_val.isMember("code")) {
                                if (err_val["code"].asInt() == -32000) { // Node is syncing - please wait
                                    std::lock_guard<std::mutex> lock(jobMutex);
                                    activeJob.reset();
                                }
                            }
                        }
                    }

                }
            }
        }
    }
#ifdef SHOW_NETWORK
    std::cout << "Network Reader is exited..." << std::endl;
#endif
}

// Network writer thread
void MinerNetwork::networkWriterThread() {
    pthread_setname_np("networkWriterThread");
    while (running.load(std::memory_order_relaxed)) {
        std::unique_lock<std::mutex> lock(requestsMutex);
        requestsCV.wait(lock, [&] {
            return !requestsPool.empty() || !running.load(std::memory_order_relaxed);
        });

        if (!running.load(std::memory_order_relaxed))
            break;

        while (!requestsPool.empty()) {
            Json::Value response = requestsPool.front();
            requestsPool.pop();
            lock.unlock();

            sendJsonMessage(response);

            lock.lock();
        }
    }
}

void MinerNetwork::sendLoginMessage(const std::string & loginName, const std::string & password, bool sendGetJobRequest) {
    std::lock_guard<std::mutex> lock(requestsMutex);
    requestsPool.push(generateLoginRequest(request_id++, loginName, password));
    if (sendGetJobRequest) {
        requestsPool.push(generateGetJobRequest(request_id.fetch_add(1)));
    }
    requestsCV.notify_one();
}

void MinerNetwork::sendKeepAliveRequest() {
    std::lock_guard<std::mutex> lock(requestsMutex);
    requestsPool.push(generateKeepAliveRequest(request_id.fetch_add(1)));
    requestsCV.notify_one();
}

void MinerNetwork::sendGetJobRequest() {
    std::lock_guard<std::mutex> lock(requestsMutex);
    requestsPool.push(generateGetJobRequest(request_id.fetch_add(1)));
    requestsCV.notify_one();
}

void MinerNetwork::submit_response(int height, int jobId, uint64_t nonce, const std::vector<uint32_t> & res_nonces) {
    {
        std::lock_guard<std::mutex> lock(jobMutex);
        if (activeJob.height > height) {
            std::cout << "Skipping outdated solution..." << std::endl;
            return;
        }
    }

    std::cout << "Submitting found solution..." << std::endl;

    Json::Value submit_request = generateSubmitRequest(request_id.fetch_add(1), edge_bits, height, jobId, nonce, res_nonces);
#ifdef SHOW_NETWORK
    std::cout << "Submitted response: " << json2str(submit_request) << std::endl;
#endif
    std::lock_guard<std::mutex> lock(requestsMutex);
    requestsPool.push(submit_request);
    requestsCV.notify_one();
}

CuckaJob MinerNetwork::getActiveJob() const {
    std::lock_guard<std::mutex> lock(jobMutex);
    CuckaJob job = activeJob;
    return job;
}

