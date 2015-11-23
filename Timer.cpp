
#include "Timer.h"

Timer::Timer() {}

void Timer::startTimer()
{
    this->startTime = clock();
    return;
}

double Timer::stopTimer()
{
    double time = (clock() - this->startTime)/CLOCKS_PER_SEC;
    this->totalTime += time;
    return time;
}

double Timer::step()
{
    return ((clock() - this->startTime)/CLOCKS_PER_SEC);
}

double Timer::getTotalTime()
{
    return this->totalTime;
}