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

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>
#include <assert.h>

// Simple but efficient thread pool
class ThreadPool {
public:
    ThreadPool(size_t);

    ThreadPool() = delete;
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool(ThreadPool&&) = delete;

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
        -> std::shared_future<std::invoke_result_t<F, Args...>>;

    ~ThreadPool();

    inline int getThreadsNum() const {return threads_num;}
private:
    int threads_num;
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

 inline ThreadPool::ThreadPool(size_t threads)
    :   threads_num(threads), stop(false)
{
    for(size_t i = 0;i<threads;++i)
        workers.emplace_back(
            [this]
            {
                // worker thread executing tasks form the Queue until it will be stopped
                for(;;)
                {
                    std::function<void()> task;

                    { // Waiting for a next task.
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock,
                                             [this]{ return this->stop || !this->tasks.empty(); });
                        if(this->stop && this->tasks.empty())
                            return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }

                    // Executing the task
                    task();
                }
            }
            );
}

template<class F, class... Args>
auto ThreadPool::enqueue(F&& f, Args&&... args)
    -> std::shared_future<std::invoke_result_t<F, Args...>>
{
    using return_type = std::invoke_result_t<F, Args...>;

    // Creating a future object for task result.
    auto task = std::make_shared< std::packaged_task<return_type()> >(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // don't allow enqueueing after stopping the pool
        if(stop) {
            assert(false);
        }

        tasks.emplace([task](){ (*task)(); });
    }
    // Send notification that a new task is available.
    condition.notify_one();
    return res;
}

inline ThreadPool::~ThreadPool()
{
    { // Stopping all threads with a flag
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    // Notify threads that waiting for a new task
    condition.notify_all();
    // Waiting until threads finish processing current task
    for(std::thread &worker: workers)
        worker.join();
}

#endif //THREAD_POOL_H
