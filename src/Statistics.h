
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
        minSlightBkz = std::min(time, minSlightBkz);
        maxSlightBkz = std::max(time, maxSlightBkz);
        sumSlightBkz += time;
    }

    void newNewEnumTime(const double time, const bool eqn)
    {
        minNewEnum = std::min(time, minNewEnum);
        maxNewEnum = std::max(time, maxNewEnum);
        sumNewEnum += time;

        if(!eqn) return;

        minNewEnumWE = std::min(time, minNewEnumWE);
        maxNewEnumWE = std::max(time, maxNewEnumWE);
        sumNewEnumWE += time;
    }

    void updateDistanceStats(const RR& theoretical, const RR&heuristic, const RR& reduced)
    {
        auto reduce_ratio = conv<double>(reduced / theoretical);

        minDistanceReduction = std::max(reduce_ratio, minDistanceReduction);
        maxDistanceReduction = std::min(reduce_ratio, maxDistanceReduction);

        if(heuristic == reduced)
            roundsWithoutReduction++;
        else
            sumDistanceReduction += reduce_ratio;
    }


    void updateRoundStats(const bool delayed, const bool eqn)
    {
        if(!delayed)
            roundsWithoutDelayedStages++;

        if(eqn)
            roundsWithEquations++;
    }

    void updateStagesCheckedForEquations(const unsigned long long stagesCheckedForEquations, const bool equationsFound)
    {
        if(equationsFound)
        {
            totalStagesCheckedForEquationsWithEquations += stagesCheckedForEquations;
            if(minStagesCheckedForEquationsWithEquations > stagesCheckedForEquations)
                minStagesCheckedForEquationsWithEquations = stagesCheckedForEquations;
            if(maxStagesCheckedForEquationsWithEquations < stagesCheckedForEquations)
                maxStagesCheckedForEquationsWithEquations = stagesCheckedForEquations;
        }
        else
        {
            totalStagesCheckedForEquationsWithoutEquations += stagesCheckedForEquations;
            if(minStagesCheckedForEquationsWithoutEquations > stagesCheckedForEquations)
                minStagesCheckedForEquationsWithoutEquations = stagesCheckedForEquations;
            if(maxStagesCheckedForEquationsWithoutEquations < stagesCheckedForEquations)
                maxStagesCheckedForEquationsWithoutEquations = stagesCheckedForEquations;
        }
    }

    void closeStatistics(const long totalRounds, const long uniqueEquations, const long duplicateEquations)
    {
        roundsTotal = totalRounds;
        eqnUniqueTotal = uniqueEquations;
        eqnDuplicates = duplicateEquations;

        avgStagesCheckedForEquationsWithEquations = 1.0 * totalStagesCheckedForEquationsWithEquations / totalRounds;
        avgStagesCheckedForEquationsWithoutEquations = 1.0 * totalStagesCheckedForEquationsWithoutEquations / totalRounds;

        avgSlightBkz = sumSlightBkz / roundsTotal;

        avgNewEnum = sumNewEnum / roundsTotal;
        avgNewEnumWE = sumNewEnumWE / roundsWithEquations;

        avgNumEqnPerRoundWithEqn = 1.0 * (eqnUniqueTotal + eqnDuplicates) / roundsWithEquations;
        avgNumUniqEqnPerRoundWithEqn =  1.0 * eqnUniqueTotal / roundsWithEquations;

        avgDistanceReduction = sumDistanceReduction / (roundsTotal - roundsWithoutReduction);
    }
};

#endif	/* STATISTICS_H */