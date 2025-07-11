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

#include "cuckatoo_pipeline.h"
#include "miner_network.h"
#include "features.h"

// Stage2 processing therad
void process_phase2(CuckatooPipeline * pipeline, int stage_id) {
    pthread_setname_np("Phase2_thread");
    while (true) {
        std::unique_lock<std::mutex> lock(pipeline->mutexes[stage_id]);
        pipeline->cond_vars[stage_id].wait(lock, [&] {
            return pipeline->stop.load(std::memory_order_relaxed) || pipeline->get_task_for_stage(stage_id)!=nullptr;
        });

        if (pipeline->stop.load(std::memory_order_relaxed))
            break;

        CuckatooSolution * task = pipeline->get_task_for_stage(stage_id);
        if (!task)
            continue;

        if ( pipeline->solver.phase2(*task) ) {
#ifdef SHOW_NETWORK
            std::cout << "Solution candidates are found." << std::endl;
#endif
            task->current_stage = stage_id+1;
            pipeline->cond_vars[stage_id+1].notify_one();
        }
        else {
#ifdef SHOW_NETWORK
            std::cout << "Empty solution." << std::endl;
#endif
            // deleting the task - no solution is found
            pipeline->delete_task_from_q(task);
        }
    }
}

// Enrich and response processing
void enrich_and_response(CuckatooPipeline * pipeline, CuckatooResponser * responser, int stage_id) {
    pthread_setname_np("Enrich_thread");
    while (true) {
        std::unique_lock<std::mutex> lock(pipeline->mutexes[stage_id]);
        pipeline->cond_vars[stage_id].wait(lock, [&] {
            return pipeline->stop.load(std::memory_order_relaxed) || pipeline->get_task_for_stage(stage_id)!=nullptr;
        });

        if (pipeline->stop.load(std::memory_order_relaxed))
            break;

        CuckatooSolution * task = pipeline->get_task_for_stage(stage_id);
        if (!task)
            continue;

        pipeline->solver.enrich_cycles(*task);
        //  Note, it is possible that solition still not found.
        //  If there are several soluitons, we have to find the best
#ifdef SHOW_NETWORK
        std::cout << "Found cycles: " << task->cycles.size() << std::endl;
#endif

        for (Cycle & cycle : task->cycles) {
            assert(cycle.cycle.size() == 42 || cycle.cycle.size() == 0);
            if (cycle.cycle.size() == 42) {
                std::sort(cycle.cycle.begin(), cycle.cycle.end());
                responser->submit_response(task->height, task->jobId, task->nonce, cycle.cycle);
            }
            else {
#ifdef SHOW_NETWORK
                std::cout << "Skipping cycle because ot size: " << cycle.cycle.size() << std::endl;
#endif
            }
        }

        // done with a task, we can release it
        pipeline->delete_task_from_q(task);
    }
}

/**
 * Init the class
 * @param responser - interfave to send found solution
 */
CuckatooPipeline::CuckatooPipeline(CuckatooResponser * responser)
{
    phase2_thread = new std::thread( process_phase2, this, 0);
    enrich_cycles_thread = new std::thread( enrich_and_response, this, responser, 1);
}

/**
 * Release the memory
 */
void CuckatooPipeline::release() {
    stop.store( true, std::memory_order_relaxed);
    // Trigger threads to finish processing
    for (int i=0; i<CUCKATOO_STAGE_NUM; i++) {
        cond_vars[i].notify_all();
    }

    if (phase2_thread) {
        phase2_thread->join();
        delete phase2_thread;
        phase2_thread = nullptr;
    }
    if (enrich_cycles_thread) {
        enrich_cycles_thread->join();
        delete enrich_cycles_thread;
        enrich_cycles_thread = nullptr;
    }

    // Q is own the content, need to delete it
    for (CuckatooSolution* q : queue) {
        delete q;
    }
    queue.clear();
}

// Request task to process
CuckatooSolution* CuckatooPipeline::get_task_for_stage(int stage_id) const {
    std::lock_guard<std::mutex> lock(queue_mutex);
    for (CuckatooSolution * task : queue) {
        if (task->current_stage == stage_id)
            return task;
    }
    return nullptr;
}

// Delete task from the Q, resease the task instance
void CuckatooPipeline::delete_task_from_q(CuckatooSolution * task) {
    std::lock_guard<std::mutex> lock(queue_mutex);
    queue.erase( std::remove_if( queue.begin(), queue.end(), [task](CuckatooSolution * t) { return t==task; } ),
        queue.end());
    delete task;
}

int CuckatooPipeline::get_queue_size() const {
    std::lock_guard<std::mutex> lock(queue_mutex);
    return queue.size();
}

/**
 * Process phase1 of a new task sync, then trigger further processing
 * @param height - need for response. Current height
 * @param jobId - need for response. Job Id
 * @param nonce - need for response. Nonce that was used to generate the hash_v4
 * @param hash_v4  - hash data, input for task
 */
void CuckatooPipeline::submit_task( int height, uint32_t jobId, uint64_t nonce, uint64_t hash_v4[4] ) {
    while( get_queue_size() > 2 ) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    FindCycleStartResult ph1_res = solver.phase1(hash_v4);
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        CuckatooSolution * solution = new CuckatooSolution(height, jobId, nonce, hash_v4, 0);
        solution->phase1_result = std::move(ph1_res);
        queue.push_back( solution );
    }
    cond_vars[0].notify_one(); // Wake up stage 0
}

