
#ifndef FILEOUTPUT_H
#define	FILEOUTPUT_H

#include "Equation.h"
#include "Statistics.h"
#include "FactoringSettings.h"

#include <NTL/RR.h>
#include <NTL/matrix.h>
#include <NTL/version.h>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <time.h>
#include <cstdlib>

#include <unistd.h>
#include <stdlib.h>
#include <cstdio>

using namespace NTL;
using namespace std;

struct FactoringSettings;

/**
 * FileOutput This class just writes stuff into some files. Add all file output
 * you want here.
 */
class FileOutput
{

private:
    /**
     * @var formatted A LaTeX file, that will contain the settings and the found
     * equations.
     */
    fstream equations;
    /**
     * @var statistics A LaTeX file, that will contain statistics to the rounds
     */
    fstream statistics;

    string eqnName;

    string statsName;

    const string alpha_2_indicator_text[3] = {"$\\alpha_{2} \\leq 0.4$",
                                              "$0.4 < \\alpha_{2} \\leq 0.65$",
                                              "$\\alpha_{2} > 0.65$"};
    const string t_indicator[3] = {"$t < 18$",
                                   "$18 \\leq t < 40$",
                                   "$t \\geq 40$"};

    void writeEqnFormatted(const Equation& eqn, const Vec<long>& primes);

    void writeEqnStatistics(const Equation& eqn, const Vec<long>& primes);

    void createDirectory();

    string getFilePrefix();

public:

    FileOutput();

    void statisticNewRound(long round);

    void statisticSlightBKZ(double slightBkz, double newEnum);

    void statisticsWriteStagesChecked(unsigned long long stagesChecked);

    void statisticsDistances(RR theoretical, RR heuristic, RR reduced);

    void statisticsDelayedStagesOnLevel(int max_level, const vector<vector<double>> &alpha_2_min,
                                        const vector<vector<vector<unsigned long long>>> &delayedAndPerformedStages,
                                        const vector<vector<vector<unsigned long long>>> &delayedStages,
                                        unsigned long long totalDelayedAndPerformedStages);

    void statisticsNewEquations(const list<Equation>& eqns, const Vec<long>& primes);

    void statisticsWriteScaledPrimes(const vector<bool> &scaledPrimes, const Vec<long> &primes);

    void writeFormattedEquationList(std::set<Equation>& eqns, const Vec<long>& primes);

    void writeSettings(const FactoringSettings &settings, long max_prime, long long int seed);

    void statisticsStrongBkzTime(double time);

    void prepareEquationTable();

    void closeEquationFile();

    void closeStatisticsFile();

    void writeSummary(const Statistics& stats, double time, long n, std::set<Equation> &uniqueEquations);

    void texToPdf();

};

#endif	/* FILEOUTPUT_H */