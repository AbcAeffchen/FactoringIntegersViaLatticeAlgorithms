
#include "Timer.h"

Timer::Timer() {}
/**
 * starts a timer. Notice: If you start the timer again without stopping
 * it first, the first timer gets overwritten and is lost.
 */
void Timer::startTimer()
{
    this->startTime = clock();
    return;
}
/**
 * stops the timer, computes the time difference in sec and adds the time to
 * an overall timer.
 * @return the time in sec
 */
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
/**
 * @return the overall time in sec
 */
double Timer::getTotalTime()
{
    return this->totalTime;
}
