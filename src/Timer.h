
#ifndef TIMER_H
#define	TIMER_H

#include <ctime>

class Timer
{
private:

    double totalTime = 0;
    double startTime;

public:
    Timer() = default;
    /**
     * starts a timer. Notice: If you start the timer again without stopping
     * it first, the first timer gets overwritten and is lost.
     */
    void startTimer()
    {
        startTime = clock();
    }

    /**
     * stops the timer, computes the time difference in sec and adds the time to
     * an overall timer.
     * @return the time in sec
     */
    double stopTimer()
    {
        double time = (clock() - startTime) / CLOCKS_PER_SEC;
        totalTime += time;
        return time;
    }

    double step()
    {
        return ((clock() - startTime) / CLOCKS_PER_SEC);
    }

    /**
     * @return the overall time in sec
     */
    double getTotalTime()
    {
        return totalTime;
    }
};

#endif	/* TIMER_H */