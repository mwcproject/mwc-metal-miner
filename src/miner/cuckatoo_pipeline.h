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

#ifndef CUCKATOO_PIPELINE_H
#define CUCKATOO_PIPELINE_H

#include "cuckatoo.h"

#define CUCKATOO_STAGE_NUM 2

class CuckatooResponser;

/**
 * Pipeline to run cuckatoo cycle. Using two threads to run.
 * Probably should switch to the tread pool, might be more efficient.
 */
class CuckatooPipeline {
public:
    /**
     * Init the class
     * @param responser - interfave to send found solution
     */
    CuckatooPipeline(CuckatooResponser * responser);
    ~CuckatooPipeline() {release();}

    /**
     * Process phase1 of a new task sync, then trigger further processing
     * @param height - need for response. Current height
     * @param jobId - need for response. Job Id
     * @param nonce - need for response. Nonce that was used to generate the hash_v4
     * @param hash_v4  - hash data, input for task
     */
    void submit_task( int height, uint32_t jobId, uint64_t nonce, uint64_t hash_v4[4] );

    /**
     * Release the memory
     */
    void release();

public:
    // Request task to process
    CuckatooSolution * get_task_for_stage(int stage_id) const;
    // Delete task from the Q, resease the task instance
    void delete_task_from_q(CuckatooSolution * task);

    // Conditional varibale to sync when a new task available
    std::condition_variable cond_vars[CUCKATOO_STAGE_NUM];
    std::mutex mutexes[CUCKATOO_STAGE_NUM];
    std::atomic<bool> stop{false};

    // Solver to process the task.
    CuckatooSolver<31, 42> solver; // C31 solver
private:
    // Size of the Q
    int get_queue_size() const;
private:
    // Q to hold tasks that are in progress
    std::deque<CuckatooSolution*> queue;
    mutable std::mutex queue_mutex; // Protects the queue

    // Stage2 thread
    std::thread * phase2_thread = nullptr;
    // Enrich found solution thread
    std::thread * enrich_cycles_thread = nullptr;
};


#endif //CUCKATOO_PIPELINE_H
