
#ifndef TIMER_H
#define	TIMER_H

#include <vector>
#include <time.h>

class Timer
{
private:

    double totalTime = 0;
    double startTime;

public:
    Timer();
    /**
     * starts a timer. Notice: If you start the timer again without stopping
     * it first, the first timer gets overwritten and is lost.
     */
    void startTimer();
    /**
     * stops the timer, computes the time difference in sec and adds the time to
     * an overall timer.
     * @return the time in sec
     */
    double stopTimer();

    double step();

    /**
     * @return the overall time in sec
     */
    double getTotalTime();
};

#endif	/* TIMER_H */