
#ifndef FILEOUTPUT_H
#define	FILEOUTPUT_H

#include "Equation.h"
#include "Statistics.h"

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

    void statisticsDistances(RR theoretical, RR heuristical, RR reduced);

    void statisticsDelayedStagesOnLevel(vector<long> delayedStagesCounter);

    void statisticsNewEquations(const list<Equation>& eqns, const Vec<long>& primes);

    void writeFormattedEquationList(std::set<Equation>& eqns, const Vec<long>& primes);

    void writeSettings(ZZ N, RR c, long accuracy_factor, int s_max, double A_start_factor, double reduce_ratio, long strong_bkz, long slight_bkz, long prime_num, long max_prime, long long int seed);

    void statisticsStrongBkzTime(double time);

    void prepareEquationTable();

    void closeEquationFile();

    void closeStatisticsFile();

    void writeSummary(const Statistics& stats, double time);

    void texToPdf();

};

#endif	/* FILEOUTPUT_H */

