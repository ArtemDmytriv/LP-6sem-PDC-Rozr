#include <timer.h>

namespace custom_timer {
void Timer::start() {
    start_time = std::chrono::high_resolution_clock::now();
}

void Timer::stop() {
    end_time = std::chrono::high_resolution_clock::now();
}

long double Timer::getDuration() {
    std::chrono::duration<long double> duration = end_time - start_time;
    return std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
}
}
