
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

    Statistics() {}

    void newSlightBkzTime(double time);

    void newNewEnumTime(double time, bool eqn);

    void updateDistanceStats(const RR& theoretical, const RR&heuristic, const RR& reduced);

    void updateRoundStats(bool delayed, bool eqn);

    void closeStatistics(long totalRounds, long uniqueEquations, long duplicateEquations);
};

#endif	/* STATISTICS_H */