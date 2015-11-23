
#ifndef TIMER_H
#define	TIMER_H

#include <vector>
#include <time.h>

class Timer
{
private:

    double totalTime = 0;
    std::vector<double> startTime;

public:
    Timer(unsigned int threads = 1);
    /**
     * starts a timer. Notice: If you start the timer again without stopping
     * it first, the first timer gets overwritten and is lost.
     */
    void startTimer(int thread = 0);
    /**
     * stops the timer, computes the time difference in sec and adds the time to
     * an overall timer.
     * @return the time in sec
     */
    double stopTimer(int thread = 0);

    double step(int thread = 0);

    /**
     * @return the overall time in sec
     */
    double getTotalTime();
};

#endif	/* TIMER_H */