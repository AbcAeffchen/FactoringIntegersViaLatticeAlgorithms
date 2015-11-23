
#include "Timer.h"

Timer::Timer(unsigned int threads) : startTime(std::vector<double>(threads)) {}

void Timer::startTimer(int thread)
{
    this->startTime[thread] = clock();
    return;
}

double Timer::stopTimer(int thread)
{
    double time = (clock() - this->startTime[thread])/CLOCKS_PER_SEC;
    this->totalTime += time;
    return time;
}

double Timer::step(int thread)
{
    return ((clock() - this->startTime[thread])/CLOCKS_PER_SEC);
}

double Timer::getTotalTime()
{
    return this->totalTime;
}