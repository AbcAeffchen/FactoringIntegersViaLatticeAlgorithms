
#include "Statistics.h"


Statistics::Statistics()
{
//    this->eqnDuplicates = 0;
//    this->eqnUniqueTotal = 0;
//    this->roundsTotal = 0;
//    this->roundsWithEquations = 0;
//    this->roundsWithoutReduction = 0;
//    this->roundsWithoutDelayedStages = 0;
//    this->avgNumEqnPerRoundWithEqn = 0;
//    this->avgNumUniqEqnPerRoundWithEqn = 0;
//
//    this->minSlightBkz = INFINITY;
//    this->avgSlightBkz = 0;
//    this->maxSlightBkz = 0;
//    this->sumSlightBkz = 0;
//
//    // total
//    this->minNewEnum = INFINITY;
//    this->avgNewEnum = 0;
//    this->maxNewEnum = 0;
//    this->sumNewEnum = 0;
//
//    // With equation
//    this->minNewEnumWE = INFINITY;
//    this->avgNewEnumWE = 0;
//    this->maxNewEnumWE = 0;
//    this->sumNewEnumWE = 0;
}

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

void Statistics::updateDistanceStats(const RR& theoretical, const RR& heuristical, const RR& reduced)
{
    double reduce_ratio = conv<double>(reduced / theoretical);
    
    this->minDistanceReduction = std::max(reduce_ratio, this->minDistanceReduction);
    this->maxDistanceReduction = std::min(reduce_ratio, this->maxDistanceReduction);
    
    if(heuristical == reduced)
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

    this->avgSlightBkz = this->sumSlightBkz / this->roundsTotal;
    
    this->avgNewEnum = this->sumNewEnum / this->roundsTotal;
    this->avgNewEnumWE = this->sumNewEnumWE / this->roundsWithEquations;
    
    this->avgNumEqnPerRoundWithEqn = 1.0 * (this->eqnUniqueTotal + this->eqnDuplicates) / this->roundsWithEquations;
    this->avgNumUniqEqnPerRoundWithEqn =  1.0 * this->eqnUniqueTotal / this->roundsWithEquations;
    
    this->avgDistanceReduction = this->sumDistanceReduction / (this->roundsTotal - this->roundsWithoutReduction);
}