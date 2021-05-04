#ifndef _TIMER_H
#define _TIMER_H
#include <chrono>

namespace custom_timer {
class Timer {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_time;
public:
    Timer() {}
    void start();
    void stop();
    long double getDuration();
};
}
#endif