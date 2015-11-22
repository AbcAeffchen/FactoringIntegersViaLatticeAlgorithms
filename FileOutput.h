
#ifndef FILEOUTPUT_H
#define	FILEOUTPUT_H

#include "Equation.h"
#include "Statistics.h"
#include "FactoringSettings.h"

#include <NTL/RR.h>
#include <NTL/matrix.h>
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
    const string t_indicator[3] = {"$t < 20$",
                                   "$20 \\leq t < 50$",
                                   "$t \\geq 50$"};

    void writeEqnFormatted(const Equation& eqn, const Vec<long>& primes);

    void writeEqnStatistics(const Equation& eqn, const Vec<long>& primes);

    void createDirectory();

    string getFilePrefix();

    void debugSeparator(string separatorText);

public:

    /**
     * @var debug A file with direct access to print debug data directly to a file.
     * equations.
     */
    fstream debug;

    FileOutput();

    void statisticNewRound(long round);

    void statisticSlightBKZ(double slightBkz, double newEnum);

    void statisticsDistances(RR theoretical, RR heuristic, RR reduced);

    void statisticsDelayedStagesOnLevel(unsigned long long max_level, const vector<vector<double>> &alpha_2_min,
                                        const vector<vector<vector<unsigned long long>>> &delayedAndPerformedStages,
                                        const vector<vector<vector<unsigned long long>>> &delayedStages,
                                        unsigned long long totalDelayedAndPerformedStages);

    void statisticsNewEquations(const list<Equation>& eqns, const Vec<long>& primes);

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