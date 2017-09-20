
#ifndef STATISTICS_H
#define	STATISTICS_H

#include <NTL/RR.h>

using namespace NTL;

struct Statistics
{
private:
    double sumSlightBkz = 0;
    double sumNewEnum = 0;
    double sumNewEnumWE = 0;
    double sumDistanceReduction = 0;

public:
    long roundsTotal = 0;
    long roundsWithEquations = 0;           // rounds that found at least one equations
    long roundsWithoutDelayedStages = 0;
    long roundsWithoutReduction = 0;

    unsigned long long totalStagesCheckedForEquationsWithEquations = 0;
    unsigned long long minStagesCheckedForEquationsWithEquations = 100000000;   // this will be reduced
    unsigned long long maxStagesCheckedForEquationsWithEquations = 0;
    double avgStagesCheckedForEquationsWithEquations = 0;
    unsigned long long totalStagesCheckedForEquationsWithoutEquations = 0;
    unsigned long long minStagesCheckedForEquationsWithoutEquations = 100000000;   // this will be reduced
    unsigned long long maxStagesCheckedForEquationsWithoutEquations = 0;
    double avgStagesCheckedForEquationsWithoutEquations = 0;

    long eqnUniqueTotal = 0;
    long eqnDuplicates = 0;
    double avgNumEqnPerRoundWithEqn = 0;
    double avgNumUniqEqnPerRoundWithEqn = 0;

    double minSlightBkz = INFINITY;
    double maxSlightBkz = 0;
    double avgSlightBkz = 0;

    double minNewEnum = INFINITY;
    double maxNewEnum = 0;
    double avgNewEnum = 0;

    double minNewEnumWE = INFINITY;
    double maxNewEnumWE = 0;
    double avgNewEnumWE = 0;

    double maxDistanceReduction = INFINITY;
    double minDistanceReduction = 0;
    double avgDistanceReduction = 0;

    Statistics() = default;

    void newSlightBkzTime(const double time)
    {
        this->minSlightBkz = std::min(time, this->minSlightBkz);
        this->maxSlightBkz = std::max(time, this->maxSlightBkz);
        this->sumSlightBkz += time;
    }

    void newNewEnumTime(const double time, const bool eqn)
    {
        this->minNewEnum = std::min(time, this->minNewEnum);
        this->maxNewEnum = std::max(time, this->maxNewEnum);
        this->sumNewEnum += time;

        if(!eqn) return;

        this->minNewEnumWE = std::min(time, this->minNewEnumWE);
        this->maxNewEnumWE = std::max(time, this->maxNewEnumWE);
        this->sumNewEnumWE += time;
    }

    void updateDistanceStats(const RR& theoretical, const RR&heuristic, const RR& reduced)
    {
        auto reduce_ratio = conv<double>(reduced / theoretical);

        this->minDistanceReduction = std::max(reduce_ratio, this->minDistanceReduction);
        this->maxDistanceReduction = std::min(reduce_ratio, this->maxDistanceReduction);

        if(heuristic == reduced)
            this->roundsWithoutReduction++;
        else
            this->sumDistanceReduction += reduce_ratio;
    }


    void updateRoundStats(const bool delayed, const bool eqn)
    {
        if(!delayed)
            this->roundsWithoutDelayedStages++;

        if(eqn)
            this->roundsWithEquations++;
    }

    void updateStagesCheckedForEquations(const unsigned long long stagesCheckedForEquations, const bool equationsFound)
    {
        if(equationsFound)
        {
            this->totalStagesCheckedForEquationsWithEquations += stagesCheckedForEquations;
            if(this->minStagesCheckedForEquationsWithEquations > stagesCheckedForEquations)
                this->minStagesCheckedForEquationsWithEquations = stagesCheckedForEquations;
            if(this->maxStagesCheckedForEquationsWithEquations < stagesCheckedForEquations)
                this->maxStagesCheckedForEquationsWithEquations = stagesCheckedForEquations;
        }
        else
        {
            this->totalStagesCheckedForEquationsWithoutEquations += stagesCheckedForEquations;
            if(this->minStagesCheckedForEquationsWithoutEquations > stagesCheckedForEquations)
                this->minStagesCheckedForEquationsWithoutEquations = stagesCheckedForEquations;
            if(this->maxStagesCheckedForEquationsWithoutEquations < stagesCheckedForEquations)
                this->maxStagesCheckedForEquationsWithoutEquations = stagesCheckedForEquations;
        }
    }

    void closeStatistics(const long totalRounds, const long uniqueEquations, const long duplicateEquations)
    {
        this->roundsTotal = totalRounds;
        this->eqnUniqueTotal = uniqueEquations;
        this->eqnDuplicates = duplicateEquations;

        this->avgStagesCheckedForEquationsWithEquations = 1.0 * this->totalStagesCheckedForEquationsWithEquations / totalRounds;
        this->avgStagesCheckedForEquationsWithoutEquations = 1.0 * this->totalStagesCheckedForEquationsWithoutEquations / totalRounds;

        this->avgSlightBkz = this->sumSlightBkz / this->roundsTotal;

        this->avgNewEnum = this->sumNewEnum / this->roundsTotal;
        this->avgNewEnumWE = this->sumNewEnumWE / this->roundsWithEquations;

        this->avgNumEqnPerRoundWithEqn = 1.0 * (this->eqnUniqueTotal + this->eqnDuplicates) / this->roundsWithEquations;
        this->avgNumUniqEqnPerRoundWithEqn =  1.0 * this->eqnUniqueTotal / this->roundsWithEquations;

        this->avgDistanceReduction = this->sumDistanceReduction / (this->roundsTotal - this->roundsWithoutReduction);
    }
};

#endif	/* STATISTICS_H */