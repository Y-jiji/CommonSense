#pragma once


#include <chrono>

namespace ONIAK {

class StopWatch{
private:
    bool running_;
    std::chrono::time_point<std::chrono::steady_clock> start_time_;
    double total_time_;
public:
    explicit StopWatch(bool with_start = true):
       running_(false), total_time_(0)
    {
        if (with_start){
            start();
        }
    }

    bool is_running() const {
        return running_;
    }

    bool start(){
        if (running_){
            // clock is already running
            return false;
        } else {
            running_ = true;
            start_time_ = std::chrono::steady_clock::now();
            return true;
        }
    }

    bool stop(){
        if (!running_){
            // clock is already stopped
            return false;
        } else {
            running_ = false;
            std::chrono::time_point<std::chrono::steady_clock> end_time = std::chrono::steady_clock::now();
            std::chrono::duration<double> diff = end_time - start_time_;
            total_time_ += diff.count();
            return true;
        }
    }

    // All time units are in seconds.
    double peek(bool stop_it = false){
        if (running_){
            stop();
        }
        if (!stop_it) start();
        return total_time_;
    }

    void reset_and_start(){
        running_ = true;
        total_time_ = 0;
        start_time_ = std::chrono::steady_clock::now();
    }
};

}

