
#include "Statistics.h"

void Statistics::newSlightBkzTime(double time)
{
    this->minSlightBkz = std::min(time, this->minSlightBkz);
    this->maxSlightBkz = std::max(time, this->maxSlightBkz);
    this->sumSlightBkz += time;
}

void Statistics::newNewEnumTime(double time, bool eqn)
{
    this->minNewEnum = std::min(time, this->minNewEnum);
    this->maxNewEnum = std::max(time, this->maxNewEnum);
    this->sumNewEnum += time;

    if(!eqn) return;

    this->minNewEnumWE = std::min(time, this->minNewEnumWE);
    this->maxNewEnumWE = std::max(time, this->maxNewEnumWE);
    this->sumNewEnumWE += time;
}

void Statistics::updateDistanceStats(const RR& theoretical, const RR&heuristic, const RR& reduced)
{
    double reduce_ratio = conv<double>(reduced / theoretical);

    this->minDistanceReduction = std::max(reduce_ratio, this->minDistanceReduction);
    this->maxDistanceReduction = std::min(reduce_ratio, this->maxDistanceReduction);

    if(heuristic == reduced)
        this->roundsWithoutReduction++;
    else
        this->sumDistanceReduction += reduce_ratio;
}

void Statistics::updateRoundStats(bool delayed, bool eqn)
{
    if(!delayed)
        this->roundsWithoutDelayedStages++;

    if(eqn)
        this->roundsWithEquations++;
}

void Statistics::closeStatistics(long totalRounds, long uniqueEquations, long duplicateEquations)
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

void Statistics::updateStagesCheckedForEquations(unsigned long long stagesCheckedForEquations, bool equationsFound)
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