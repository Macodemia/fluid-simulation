#pragma once
#include <cstdint>
#include <vector>
#include <memory>
#include <thread>
#include <iostream>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <Eigen/Dense>
#include "ThreadPool.hpp"

typedef std::int8_t  int8;
typedef std::int16_t int16;
typedef std::int32_t int32;
typedef std::int64_t int64;

typedef std::uint8_t  uint8;
typedef std::uint16_t uint16;
typedef std::uint32_t uint32;
typedef std::uint64_t uint64;

namespace learnSPH {

#if 0
// Threadpool from https://github.com/f-squirrel/thread_pool
static thread_pool::ThreadPool threadPool;

template<typename Function>
void forEachParticleDo(int numberOfThreads, std::vector<Eigen::Vector3d>& particles, Function callback) {
    const unsigned int particlesPerThread = particles.size() / numberOfThreads;

    std::vector<std::shared_ptr<std::future<void>>> tasks;
    tasks.reserve(numberOfThreads);
    for (int threadIndex = 0; threadIndex < numberOfThreads; threadIndex++) {
        const unsigned int startIndex = particlesPerThread * threadIndex;
        const unsigned int endIndex = threadIndex == numberOfThreads - 1 ? particles.size() : startIndex + particlesPerThread;

        std::shared_ptr<std::future<void>> task = std::make_shared<std::future<void>>(threadPool.submit([startIndex, endIndex, &callback]() -> void {
            for (unsigned int i = startIndex; i < endIndex; i++) {
                callback(i);
            }
        }));

        tasks.push_back(task);
    }

    for (auto task : tasks) {
        task->wait();
    }
}
#endif

#if 1
// Without Threadpool
template<typename Function>
void forEachParticleDo(int numberOfThreads, std::vector<Eigen::Vector3d>& particles, Function callback) {
    const unsigned int particlesPerThread = particles.size() / numberOfThreads;
    std::vector<std::shared_ptr<std::thread>> threads;

    for (int threadIndex = 0; threadIndex < numberOfThreads; threadIndex++) {
        const unsigned int startIndex = particlesPerThread * threadIndex;
        const unsigned int endIndex = threadIndex == numberOfThreads - 1 ? particles.size() : startIndex + particlesPerThread;

        std::shared_ptr<std::thread> thread = std::make_shared<std::thread>([startIndex, endIndex, &callback]()->void{
            for(unsigned int i = startIndex; i < endIndex; i++){
                callback(i);
            }
        });
        threads.emplace_back(thread);
    }

    for (auto& thread : threads) {
        thread->join();
    }
}
#endif

#if 0
// Weird hacky selfmade half threadpool
struct ParticleLoopTasks {
    std::condition_variable conditionVariable;
    bool readyToRun = false;
    unsigned int startIndex;
    unsigned int endIndex;
    std::function<void(unsigned int)> callback;
};

static std::mutex taskMutex;
// We need to construct with a size as it cannot be made larger later (condition_variable cannot be copied)
// NOTE: This means we cannot have more than 100 threads!
static std::vector<ParticleLoopTasks> tasks(100); 
static std::vector<std::shared_ptr<std::thread>> threads;

template<typename Function>
void forEachParticleDo(int numberOfThreads, std::vector<Eigen::Vector3d>& particles, Function callback) {
    // NOTE: Assuming that numberOfThreads will stay the same during program execution (which it currently does)!
    if (threads.empty()) {
        threads.reserve(numberOfThreads);

        for (int i = 0; i < numberOfThreads; i++) {
            threads.push_back(std::make_shared<std::thread>([i]() -> void {
                while(true) {
                    //std::unique_lock<std::mutex> lock(taskMutex);
                    //tasks[i].conditionVariable.wait(lock);
                    if (tasks[i].readyToRun) {
                        for(unsigned int particleIndex = tasks[i].startIndex; particleIndex < tasks[i].endIndex; particleIndex++) {
                            tasks[i].callback(particleIndex);
                        }
                        tasks[i].readyToRun = false;
                    }
                    std::this_thread::sleep_for(std::chrono::nanoseconds(200));
                }
            }));
        }
    }

    const unsigned int particlesPerThread = particles.size() / numberOfThreads;

    for (int threadIndex = 0; threadIndex < numberOfThreads; threadIndex++) {
        const unsigned int startIndex = particlesPerThread * threadIndex;
        const unsigned int endIndex = threadIndex == numberOfThreads - 1 ? particles.size() : startIndex + particlesPerThread;

        tasks[threadIndex].startIndex = startIndex;
        tasks[threadIndex].endIndex = endIndex;
        tasks[threadIndex].callback = callback;
        tasks[threadIndex].readyToRun = true;
        //tasks[threadIndex].conditionVariable.notify_all();
    }

    bool allFinished = false;
    while (!allFinished) {
        allFinished = true;
        for (int i = 0; i < threads.size(); i++) {
            allFinished = allFinished & !tasks[i].readyToRun;
        }
    }
}
#endif

};
