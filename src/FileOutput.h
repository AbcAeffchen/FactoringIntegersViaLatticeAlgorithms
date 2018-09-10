
#ifndef FILEOUTPUT_H
#define FILEOUTPUT_H

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
#include <ctime>
#include <cstdlib>

#include <unistd.h>
#include <cstdlib>
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
     * @var summary A LaTeX file, that will contain the settings, the found equations and
     * the statistics summary.
     */
    fstream summary;
    /**
     * @var statistics A LaTeX file, that will contain statistics to the rounds
     */
    fstream statistics;

    string summaryName;

    string statsName;

    const string alpha_2_indicator_text[3] = {"$\\alpha_{2} \\leq 0.4$",
                                              "$0.4 < \\alpha_{2} \\leq 0.65$",
                                              "$\\alpha_{2} > 0.65$"};
    const string t_indicator[3] = {"$t < 18$",
                                   "$18 \\leq t < 40$",
                                   "$t \\geq 40$"};

    void writeEqnFormatted(fstream& file, const Equation& eqn, const Vec<long>& primes)
    {
        file << eqn.round << " & ";
        file << eqn.counter << " & ";

        if(eqn.fromContinuedFraction)
            file << "\\checkmark";
        else
            file << "\\cross";

        bool nEmp = false;

        file << "& $";

        for(long i = 1; i <= primes.length(); i++)
        {
            if(eqn.e(i) > 0)
            {
                if(nEmp)
                {
                    file << " \\cdot ";
                }
                file << primes(i);
                if(eqn.e(i) > 1)
                {
                    file << "^{" << eqn.e(i) << "}";
                }
                nEmp = true;
            }
        }
        if(!nEmp)
        {
            file << "1";
        }
        file << "$ & $" << eqn.v << "$ & $";
        nEmp = false;

        for(long i = 1; i <= primes.length(); i++)
        {
            if(eqn.e(i) < 0)
            {
                if(nEmp)
                {
                    file << " \\cdot ";
                }
                file << primes(i);
                if(eqn.e(i) < -1)
                {
                    file << "^{" << -eqn.e(i) << "}";
                }
                nEmp = true;
            }
        }
        if(!nEmp)
        {
            file << "1";
        }
        file << "$ & $";
        if(eqn.e(primes.length()+1) == 1)
            file << "-1";
        else
            file << "1";
        file << "$\\\\\n";
    }

    // todo unify
    void writeEqnStatistics(const Equation& eqn, const Vec<long>& primes)
    {
        bool nEmp = false;
        statistics << eqn.reduced << " & "
                   << eqn.level << " & ";

        statistics << "$";

        for(long i = 1; i <= primes.length(); i++)
        {
            if(eqn.e(i) > 0)
            {
                if(nEmp)
                    statistics << " \\cdot ";
                statistics << primes(i);
                if(eqn.e(i) > 1)
                {
                    statistics << "^{" << eqn.e(i) << "}";
                }
                nEmp = true;
            }
        }
        if(!nEmp)
        {
            statistics << "1";
        }
        statistics << "$ & $"<< eqn.v << "$ & $";
        nEmp = false;

        for(long i = 1; i <= primes.length(); i++)
        {
            if(eqn.e(i) < 0)
            {
                if(nEmp)
                {
                    statistics << " \\cdot ";
                }
                statistics << primes(i);
                if(eqn.e(i) < -1)
                {
                    statistics << "^{" << -eqn.e(i) << "}";
                }
                nEmp = true;
            }
        }
        if(!nEmp)
        {
            statistics << "1";
        }
        statistics << "$ & $" <<  eqn.time << "$s & ";

        if(eqn.fromContinuedFraction)
            statistics << "\\checkmark";
        else
            statistics << "\\cross";

        statistics << "\\\\\n";
    }

    void createDirectory()
    {
        if(system("mkdir output") == 0)
            cout << "output directory created" << endl;
    }

    /**
     * @return Current date/time, format is YYYY-MM-DD.HH:mm:ss
     */
    string getFilePrefix()
    {
        time_t now = time(nullptr);
        auto tstruct = *localtime(&now);
        char buf[80];
        strftime(buf, sizeof(buf), "%Y-%m-%d_%H-%M-%S", &tstruct);

        return buf;
    }

public:

    FileOutput()
    {
        createDirectory();

        string filePrefix = getFilePrefix();
        std::stringstream formattedName;
        formattedName << "./output/" << filePrefix << "_summary.tex";
        summaryName = filePrefix + "_summary";
        std::stringstream statisticsName;
        statisticsName << "./output/" << filePrefix << "_statistics.tex";
        statsName = filePrefix + "_statistics";

        summary.open(formattedName.str().c_str(), ios::out);
        statistics.open(statisticsName.str().c_str(), ios::out);

        // Init formatted file
        summary << "\\documentclass[a4paper,twoside,10pt]{report}\n\n"
                   "\\usepackage[a4paper, left=1cm, right=1cm, top=2cm,bottom=2cm]{geometry}\n"
                   "\\usepackage{fancyhdr,amsmath,longtable,setspace,graphicx,booktabs,amssymb,pifont}\n"
                   "\\DeclareMathOperator{\\sign}{\\text{sign}}\n"
                   "\\def\\qqquad{\\qquad\\qquad}\n"
                   "\\renewcommand{\\checkmark}{\\text{\\ding{51}}}\n"
                   "\\newcommand{\\cross}{\\text{\\ding{55}}}\n"
                   "\\pagestyle{empty}\n\n"
                   "\\begin{document}\n"
                   "\\onehalfspacing\n";

        statistics << "\\documentclass[a4paper,twoside,10pt]{report}\n\n"
                      "\\usepackage[a4paper, left=1cm, right=1cm, top=2cm,bottom=2cm]{geometry}\n"
                      "\\usepackage{fancyhdr,amsmath,longtable,setspace,booktabs,graphicx,xcolor,array,amssymb,pifont}\n"
                      "\\DeclareMathOperator{\\sign}{\\text{sign}}\n"
                      "\\def\\qqquad{\\qquad\\qquad}\n"
                      "\\pagestyle{empty}\n\n"
                      "\\newcolumntype{C}[1]{>{\\centering\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}\n"
                      "\\newcolumntype{R}[1]{>{\\raggedleft\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}\n"
                      "\\renewcommand{\\checkmark}{\\text{\\ding{51}}}\n"
                      "\\newcommand{\\cross}{\\text{\\ding{55}}}\n"
                      "\\newcommand{\\green}[1]{\\textcolor[HTML]{2B9D0E}{#1}}\n"
                      "\\newcommand{\\red}[1]{\\textcolor[HTML]{A80000}{#1}}\n"
                      "\\newcommand{\\grey}[1]{\\textcolor{black!40}{#1}}\n"
                      "\\begin{document}\n"
                      "\\onehalfspacing\n\n";
    }

    void statisticNewRound(unsigned long round)
    {
        statistics << "\\section*{Round " << round << "}\n\\begin{longtable}{p{13cm}p{4cm}}\n";
    }

    void statisticSlightBKZ(double slightBkz, double newEnum)
    {
        statistics << "& \\textbf{Distances:}\\newline\n"
                      "\\resizebox{\\linewidth}{!}{%\n"
                      "\\begin{tabular}[t]{ll}\n"
                      "Slight BKZ:&" << slightBkz << "s\\\\\n"
                      "NewEnum:&" << newEnum << "s\\\\\\\\\n";
    }

    void statisticsWriteStagesChecked(unsigned long long stagesChecked)
    {
        statistics << "Checked Stages: & " << stagesChecked << "\\\\\\\\\n";
    }

    void statisticsDistances(const RR& theoretical, const RR& heuristic, const RR& reduced)
    {
        statistics
            << "$A=\\frac{1}{4}\\sum_{i=1}^n r_{ii}^2$:&" << trunc(theoretical) << "\\\\\n"
               "$\\frac{1}{5}A$:&" << trunc(heuristic) << "\\\\\n"
               "Reduced to $A'$:&" << trunc(reduced) << "\\\\\n"
               "$A'/A$:&" << conv<double>(reduced/theoretical)
            << "\n\\end{tabular}}\\end{longtable}\n";
    }

    void statisticsDelayedStagesOnLevel(int max_level,
                                        const vector<vector<vector<unsigned long long>>>& delayedAndPerformedStages,
                                        const vector<vector<vector<unsigned long long>>>& delayedStages,
                                        unsigned long long totalDelayedAndPerformedStages)
    {
        vector<vector<long long>> sumDelayedStagesOverAlpha2(3,vector<long long>(max_level- 10,0));
        vector<vector<long long>> sumDelayedStagesOverT(3,vector<long long>(max_level- 10,0));
        vector<long long> totalSumDelayedStages = vector<long long>(max_level- 10,0);
        long long totalDelayedStages = 0;
        vector<vector<long long>> sumDelayedAndPerformedStagesOverAlpha2(3,vector<long long>(max_level- 10,0));
        vector<vector<long long>> sumDelayedAndPerformedStagesOverT(3,vector<long long>(max_level- 10,0));
        vector<long long> totalSumDelayedAndPerformedStages = vector<long long>(max_level- 10,0);

        for(long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
            for(long t_indicator = 0; t_indicator < 3; t_indicator++)
                for(long level = 0; level < max_level- 10; level++)
                {
                    sumDelayedAndPerformedStagesOverAlpha2[t_indicator][level] += delayedAndPerformedStages[alpha_2_indicator][t_indicator][level];
                    sumDelayedAndPerformedStagesOverT[alpha_2_indicator][level] += delayedAndPerformedStages[alpha_2_indicator][t_indicator][level];
                    sumDelayedStagesOverAlpha2[t_indicator][level] += delayedStages[alpha_2_indicator][t_indicator][level];
                    sumDelayedStagesOverT[alpha_2_indicator][level] += delayedStages[alpha_2_indicator][t_indicator][level];
                }

        for(long level = 0; level < max_level- 10; level++)
            for(long i = 0; i < 3; i++)
            {
                totalSumDelayedAndPerformedStages[level] += sumDelayedAndPerformedStagesOverT[i][level];
                totalSumDelayedStages[level] += sumDelayedStagesOverT[i][level];
            }

        for(long level = 0; level < max_level - 10; level++)
            totalDelayedStages += totalSumDelayedStages[level];

        statistics << "\\textbf{Delayed Stages}\\newline\\resizebox{\\linewidth}{!}{%\n"
                      "\\begin{tabular}[t]{c|C{4cm}C{4cm}C{4cm}|C{4cm}}\n"
                      "\\textbf{\\green{" << totalDelayedAndPerformedStages
                   << "}}, \\textbf{\\red{" << totalDelayedStages << "}}& " << t_indicator[0]
                   << "&" << t_indicator[1] << "&"  << t_indicator[2] << "&\\\\\\midrule\n";

        for(long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
        {
            statistics << alpha_2_indicator_text[alpha_2_indicator];
            for(long t_indicator = 0; t_indicator < 3; t_indicator++)
            {
                // delayed and performed (green)
                statistics << "& \\green{";
                statistics << delayedAndPerformedStages[alpha_2_indicator][t_indicator][0];
                for (long level = 1; level < max_level - 10; level++)
                {
                    statistics << "," <<
                               delayedAndPerformedStages[alpha_2_indicator][t_indicator][level];
                }

                statistics << "} \\newline\n\\red{";
                // delayed stages (red)
                statistics << delayedStages[alpha_2_indicator][t_indicator][0];
                for (long level = 1; level < max_level - 10; level++)
                {
                    statistics << "," <<
                               delayedStages[alpha_2_indicator][t_indicator][level];
                }

                statistics << "}";
            }

            statistics << "& \\green{";
            statistics << sumDelayedAndPerformedStagesOverT[alpha_2_indicator][0];
            for (long level = 1; level < max_level - 10; level++)
            {
                statistics << "," <<
                           sumDelayedAndPerformedStagesOverT[alpha_2_indicator][level];
            }

            statistics << "} \\newline\n\\red{" << sumDelayedStagesOverT[alpha_2_indicator][0];
            for (long level = 1; level < max_level - 10; level++)
            {
                statistics << "," << sumDelayedStagesOverT[alpha_2_indicator][level];
            }
            statistics << "}\\\\\n";
        }

        statistics << "\\midrule\n";

        for(long t_indicator = 0; t_indicator < 3; t_indicator++)
        {
            statistics << "& \\green{" << sumDelayedAndPerformedStagesOverAlpha2[t_indicator][0];
            for(long level = 1; level < max_level- 10;level++)
            {
                statistics << "," << sumDelayedAndPerformedStagesOverAlpha2[t_indicator][level];
            }
            statistics << "} \\newline\\red{" << sumDelayedStagesOverAlpha2[t_indicator][0];
            for(long level = 1; level < max_level- 10;level++)
            {
                statistics << "," << sumDelayedStagesOverAlpha2[t_indicator][level];
            }
            statistics << "}";
        }

        statistics << "& \\green{" << totalSumDelayedAndPerformedStages[0];
        for(long level = 1; level < max_level- 10; level++)
        {
            statistics << "," << totalSumDelayedAndPerformedStages[level];
        }
        statistics << "} \\newline\\red{" << totalSumDelayedStages[0];
        for(long level = 1; level < max_level- 10; level++)
        {
            statistics << "," << totalSumDelayedStages[level];
        }

        statistics << "}\\end{tabular}}\n";
    }

    void statisticsNewEquations(const list<Equation>& eqns, const Vec<long>& primes)
    {
        if(eqns.empty())
        {
            statistics << "\\subsection*{No Equations found}\n";
            return;
        }

        statistics << "\\subsection*{" << eqns.size() << " new equations}\n"
                   << "\\begin{longtable}{p{1.5cm}lp{5.0cm}R{3.5cm}p{3.0cm}rc}\n"
                   << "\\toprule\n"
                   << "dist& level & $u$ & $v$ & $|u-vN|$ & time until & CF\\\\\\midrule\n"
                   << "\\endfirsthead\n"
                   << "\\toprule\n"
                   << "dist& level & $u$ & $v$ & $|u-vN|$ & time until & CF\\\\\\midrule\n"
                   << "\\endhead\n";

        for(const auto& eqn : eqns)
        {
            writeEqnStatistics(eqn, primes);
        }

        statistics << "\\end{longtable}\n";
    }

    void statisticsWriteScaledPrimes(const vector<bool>& scaledPrimes, const Vec<long>& primes)
    {
        statistics << "\\paragraph*{Scaled prime numbers:}";
        long n = primes.length();
        for(long i = 0; i < n; i++)
        {
            if(scaledPrimes[i])
                statistics << primes[i] << " ";
        }

    }

    void writeFormattedEquationList(std::set<Equation>& eqns, const Vec<long>& primes)
    {
        list<Equation> eqnList(eqns.begin(),eqns.end());

        eqnList.sort(sort_equations);

        summary << "\\newpage\n"
                   "\\section*{Equations}\n"
                   "%\\begin{landscape}\n"
                   "\\begin{longtable}{rrcp{5cm}rp{5cm}r}\n"
                   "\\toprule\n"
                   "round & \\# & CF & $u$ & $v$ & $|u-vN|$ & $\\sign(u-vN)$\\\\\\midrule\n"
                   "\\endfirsthead\n"
                   "\\toprule\n"
                   "round & \\# & CF & $u$ & $v$ & $|u-vN|$ & $\\sign(u-vN)$\\\\\\midrule\n"
                   "\\endhead\n";

        statistics << "\\newpage\n"
                      "\\section*{Equations}\n"
                      "%\\begin{landscape}\n"
                      "\\begin{longtable}{rrcp{5cm}rp{5cm}r}\n"
                      "\\toprule\n"
                      "round & \\# & CF & $u$ & $v$ & $|u-vN|$ & $\\sign(u-vN)$\\\\\\midrule\n"
                      "\\endfirsthead\n"
                      "\\toprule\n"
                      "round & \\# & CF & $u$ & $v$ & $|u-vN|$ & $\\sign(u-vN)$\\\\\\midrule\n"
                      "\\endhead\n";

        for(const auto& it : eqnList)
        {
            writeEqnFormatted(summary, it, primes);
            writeEqnFormatted(statistics, it, primes);
        }

        summary << "\\end{longtable}\n%\\end{landscape}\n";

        statistics << "\\end{longtable}\n%\\end{landscape}\n";
    }

    void writeSettings(const FactoringSettings& settings, long max_prime, long long int seed)
    {
        // todo unify
        summary << "\\section*{Settings}\n"
                   "\\begin{tabular}{ll}\n"
                   "NTL Version: & " << NTL_VERSION << "\\\\"
                #ifdef __GNUC__
                   "GCC: &" << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << "\\\\\n"
                #endif
                   "$N$ & $" << settings.N << "\\approx 10^{" << round(log(settings.N)/log(10)) << "}$ (ca. " << round(log(settings.N)/log(2)) << " Bits)\\\\"
                   "$c$ & $" << settings.c << "$\\\\\n"
                   "$n$ & $ " << settings.n << "$\\\\\n"
                   "$p_n$ & $ " << max_prime << "$\\\\\n"
                   "Pruning Level: & $" << settings.max_level << "$\\\\\n"
                   "Basis accuracy factor: & $" << settings.accuracy_factor << "$\\\\\n"
                   "Reduce ratio: & $" << settings.reduce_ratio << "$\\\\\n"
                   "Start reduction factor: & $" << settings.A_start_factor << "$\\\\\n"
                   "Strong BKZ block size: & $" << settings.strong_bkz << "$\\\\\n"
                   "Slight BKZ block size: & $" << settings.slight_bkz << "$\\\\\n"
                   "Seed for random number generator: & $" << seed << "$\\\\\n"
                   "Scaling type: & ";

        switch(settings.scalingType)
        {
            case FS_SCALING_MIXED: summary << "Use various scaleing types:\\\\"
                                              "& (\\phantom{7}6\\% of rounds): Scale every row with probability 1/4\\\\"
                                              "& (\\phantom{7}8\\% of rounds): Scale the first n/2 rows with probability 1/4,\\\\& \\phantom{(76\\% of rounds): } the n/2 last rows with probability 1/2\\\\"
                                              "& (\\phantom{7}8\\% of rounds): Scale the first n/2 rows with probability 1/2,\\\\& \\phantom{(76\\% of rounds): } the n/2 last rows with probability 1/4\\\\"
                                              "& (78\\% of rounds): Scale every row with probability 1/2\\\\";
                break;
            case FS_SCALING_ONE_FOURTH: summary << "Scale every row with probability 1/4\\\\";
                break;
            case FS_SCALING_THREE_FORTH: summary << "Scale every row with probability 3/4\\\\";
                break;
            case FS_SCALING_ONE_FOURTH_ONE_HALF: summary << "Scale the first n/2 rows with probability 1/4,\\\\& the n/2 last rows with probability 1/2\\\\";
                break;
            case FS_SCALING_ONE_HALF_ONE_FOURTH: summary << "Scale the first n/2 rows with probability 1/2,\\\\& the n/2 last rows with probability 1/4\\\\";
                break;
            default: summary << "Scale every row with probability 1/2\\\\";   // FS_SCALING_ONE_HALF
        }

#ifndef FS_CCF
        summary << "\nUse Continued Fractions (CF): & " << (settings.useContinuedFractions ? "yes" : "no") << "\\\\\n"
#else
        summary << "\nUse Continued Fractions (CF): & " << (settings.useContinuedFractions ? "yes (centered)" : "no") << "\\\\\n"
#endif
                   "\\end{tabular}\\\\\n";

        statistics << "\\section*{Settings}\n"
                      "\\begin{tabular}{ll}\n"
                      "NTL Version: & " << NTL_VERSION << "\\\\"
                   #ifdef __GNUC__
                      "GCC: &" << __GNUC__ << "." <<__GNUC_MINOR__ << "."  << __GNUC_PATCHLEVEL__ << "\\\\\n"
                   #endif
                      "$N$ & $" << settings.N << "\\approx 10^{" << round(log(settings.N)/log(10)) << "}$ (ca. " << round(log(settings.N)/log(2)) << " Bits)\\\\"
                      "$c$ & $" << settings.c << "$\\\\\n"
                      "$n$ & $ " << settings.n << "$\\\\\n"
                      "$p_n$ & $ " << max_prime << "$\\\\\n"
                      "Pruning Level: & $" << settings.max_level << "$\\\\\n"
                      "Basis accuracy factor: & $" << settings.accuracy_factor << "$\\\\\n"
                      "Reduce ratio: & $" << settings.reduce_ratio << "$\\\\\n"
                      "Start reduction factor: & $" << settings.A_start_factor << "$\\\\\n"
                      "Strong BKZ block size: & $" << settings.strong_bkz << "$\\\\\n"
                      "Slight BKZ block size: & $" << settings.slight_bkz << "$\\\\\n"
                      "Seed for random number generator: & $" << seed << "$\\\\\n"
                      "Scaling type: & ";
        switch(settings.scalingType)
        {
            case FS_SCALING_MIXED: statistics << "Use various scaleing types:\\\\"
                                                 "& (\\phantom{7}6\\% of rounds): Scale every row with probability 1/4\\\\"
                                                 "& (\\phantom{7}8\\% of rounds): Scale the first n/2 rows with probability 1/4,\\\\& \\phantom{(76\\% of rounds): } the n/2 last rows with probability 1/2\\\\"
                                                 "& (\\phantom{7}8\\% of rounds): Scale the first n/2 rows with probability 1/2,\\\\& \\phantom{(76\\% of rounds): } the n/2 last rows with probability 1/4\\\\"
                                                 "& (78\\% of rounds): Scale every row with probability 1/2\\\\";
                break;
            case FS_SCALING_ONE_FOURTH: statistics << "Scale every row with probability 1/4\\\\";
                break;
            case FS_SCALING_THREE_FORTH: statistics << "Scale every row with probability 3/4\\\\";
                break;
            case FS_SCALING_ONE_FOURTH_ONE_HALF: statistics << "Scale the first n/2 rows with probability 1/4,\\\\& the n/2 last rows with probability 1/2\\\\";
                break;
            case FS_SCALING_ONE_HALF_ONE_FOURTH: statistics << "Scale the first n/2 rows with probability 1/2,\\\\& the n/2 last rows with probability 1/4\\\\";
                break;
            default: statistics << "Scale every row with probability 1/2\\\\";    // FS_SCALING_ONE_HALF
        }

#ifndef FS_CCF
        statistics << "\nUse Continued Fractions (CF): & " << (settings.useContinuedFractions ? "yes" : "no") << "\\\\\n"
#else
        statistics << "\nUse Continued Fractions (CF): & " << (settings.useContinuedFractions ? "yes (centered)" : "no") << "\\\\\n"
#endif
                      "\\end{tabular}\\\\\n";

        // todo add explanation for tables to statistics file
    }

    void statisticsStrongBkzTime(double time)
    {
        statistics << "\\paragraph{Strong BKZ:}" << time << "s\n";
    }

    // todo unify
    void closeEquationFile()
    {
        summary << "\\end{document}\n";
        summary.close();
    }

    void closeStatisticsFile()
    {
        statistics << "\\end{document}\n";
        statistics.close();
    }

    void writeSummary(const Statistics& stats, double time, long n, std::set<Equation>& uniqueEquations)
    {
        list<Equation> eqns(uniqueEquations.begin(),uniqueEquations.end()),
            eqnList,
            eqnList_cf;

        RR v_avg = conv<RR>(0), v_cf_avg = conv<RR>(0);
        ZZ v_min = conv<ZZ>(0), v_cf_min = conv<ZZ>(0), v_max = conv<ZZ>(0), v_cf_max = conv<ZZ>(0), v_median = conv<ZZ>(0), v_cf_median = conv<ZZ>(0);
        for(const auto& eqn : eqns)
        {
            if(eqn.fromContinuedFraction)
            {
                v_cf_avg += conv<RR>(eqn.v);
                eqnList_cf.push_front(eqn);
            }
            else
            {
                eqnList.push_front(eqn);
                v_avg += conv<RR>(eqn.v);
            }
        }

        if(!eqnList.empty())
            v_avg /= eqnList.size();
        if(!eqnList_cf.empty())
            v_cf_avg /= eqnList_cf.size();

        eqnList.sort(sort_equations_by_v);
        eqnList_cf.sort(sort_equations_by_v);

        v_min = eqnList.front().v;
        v_max = eqnList.back().v;
        if(!eqnList_cf.empty())
        {
            v_cf_min = eqnList_cf.front().v;
            v_cf_max = eqnList_cf.back().v;
        }

        if(!eqnList.empty())
        {
            if (eqnList.size() % 2 == 1)
            {
                auto it = eqnList.begin();
                std::advance(it, (eqnList.size() - 1) / 2);
                v_median = it->v;
            }
            else
            {
                auto it = eqnList.begin();
                std::advance(it, eqnList.size() / 2);
                v_median = it->v;
                v_median += (--it)->v;
                v_median /= 2;
            }
        }

        if(!eqnList_cf.empty())
        {
            if (eqnList_cf.size() % 2 == 1)
            {
                auto it = eqnList_cf.begin();
                std::advance(it, (eqnList_cf.size() - 1) / 2);
                v_cf_median = it->v;
            }
            else
            {
                auto it = eqnList_cf.begin();
                std::advance(it, eqnList_cf.size() / 2);
                v_cf_median = it->v;
                v_cf_median += (--it)->v;
                v_cf_median /= 2;
            }
        }

        statistics << "\n\n\\newpage\n"
                      "\\section*{Summary}\n"
                      "\\resizebox{\\textwidth}{!}{\n"
                      "\\begin{tabular}{llll}\n"
                      "\\toprule\n"

                      "Runtime (total): & \\textbf{" << time << "s}\\\\\n"
                      "Runtime (per unique equation): & \\textbf{" << time/stats.eqnUniqueTotal << "s}\\\\\n"
                      "Time to factor $N$: & \\textbf{" << (time/stats.eqnUniqueTotal * (n+2)) << "s}\\\\\\midrule[0.05pt]\n"

                      "Rounds (total): & $" << stats.roundsTotal << "$&\n"
                      "Unique Equations found (duplicates) & $" << stats.eqnUniqueTotal << "$ ($" << stats.eqnDuplicates  << "$)\\\\\n"
                      "Rounds without distance reduction: & " << stats.roundsWithoutReduction << " (" << round(1.0 * stats.roundsWithoutReduction/stats.roundsTotal * 10000) /100.0 << "\\%)&\n"
                      "Unique equations per round & $" << 1.0 * stats.eqnUniqueTotal/stats.roundsTotal << "$\\\\\n"
                      "Rounds without delayed stages: & " << stats.roundsWithoutDelayedStages << " (" << round(1.0 * stats.roundsWithoutDelayedStages/stats.roundsTotal * 10000) /100.0 << "\\%)&\n"
                      "Unique equations per round with 1+ eqn. & $" << stats.avgNumUniqEqnPerRoundWithEqn << "$\\\\\n"
                      "Rounds with equations: & " << stats.roundsWithEquations << " (" << round(1.0 * stats.roundsWithEquations/stats.roundsTotal * 10000) /100.0 << "\\%)&\n"
                      "Total equations per round with 1+ eqn.   & $" << stats.avgNumEqnPerRoundWithEqn << "$\\\\\\midrule[0.05pt]\n"

                      "Total stages checked for equ. (rounds w. Equ.): & " << stats.totalStagesCheckedForEquationsWithEquations << "& Total stages checked for equ. (rounds wo. Equ.): & " << stats.totalStagesCheckedForEquationsWithoutEquations << "\\\\\n"
                      "Min. stages checked for equ. (rounds w. Equ.): & " << stats.minStagesCheckedForEquationsWithEquations << "& Min. stages checked for equ. (rounds wo. Equ.): & " << stats.minStagesCheckedForEquationsWithoutEquations << "\\\\\n"
                      "Max. stages checked for equ. (rounds w. Equ.): & " << stats.maxStagesCheckedForEquationsWithEquations << "& Max. stages checked for equ. (rounds wo. Equ.): & " << stats.maxStagesCheckedForEquationsWithoutEquations << "\\\\\n"
                      "Avg. stages checked for equ. (rounds w. Equ.): & " << round(stats.avgStagesCheckedForEquationsWithEquations * 100.0)/100.0 << "& Avg. stages checked for equ. (rounds wo. Equ.): & " << round(stats.avgStagesCheckedForEquationsWithoutEquations * 100.0)/100.0 << "\\\\\\midrule[0.05pt]\n"

                      "Min. time NewEnum (total): & " << stats.minNewEnum << "s&\n"
                      "Min. time NewEnum (w.eqn.): & " << stats.minNewEnumWE << "s\\\\\n"
                      "Max. time NewEnum: & " << stats.maxNewEnum << "s&\n"
                      "Max. time NewEnum: & " << stats.maxNewEnumWE << "s\\\\\n"
                      "Avg. time NewEnum: & " << stats.avgNewEnum << "s&\n"
                      "Avg. time NewEnum: & " << stats.avgNewEnumWE << "s\\\\\\midrule[0.05pt]\n"

                      "Min. time slight BKZ: & " << stats.minSlightBkz << "s&\n"
                      "Max. distance reduction: & " << stats.maxDistanceReduction << "\\\\\n"
                      "Max. time slight BKZ: & " << stats.maxSlightBkz << "s&\n"
                      "Min. distance reduction: & " << stats.minDistanceReduction << "\\\\\n"
                      "Avg. time slight BKZ: & " << stats.avgSlightBkz << "s&\n"
                      "Avg. distance reduction: & " << stats.avgDistanceReduction << "\\\\\\midrule[0.05pt]\n"

                      "Equations from NewEnum only: &" << eqnList.size() << " (" << round(eqnList.size() * 1000.0 / stats.eqnUniqueTotal) / 10.0
                   << "\\%) & Equations from CF: & " << eqnList_cf.size() << " (" << round(eqnList_cf.size() * 1000.0 / stats.eqnUniqueTotal)/10.0 << "\\%)\\\\\n"
                      "Min. $v$ value (without CF): & " << v_min << " & Min. $v$ value (CF only): & " << v_cf_min << "\\\\\n"
                      "Max. $v$ value (without CF): & " << v_max << " & Max. $v$ value (CF only): & " << v_cf_max << "\\\\\n"
                      "Avg. $v$ value (without CF): & " << (round(v_avg * 100.0)/100.0) << " & Avg. $v$ value (CF only): & " << (round(v_cf_avg * 100.0)/100.0) << "\\\\\n"
                      "Median $v$ value (without CF): & " << v_median << " & Median $v$ value (CF only): & " << v_cf_median << "\\\\\\bottomrule\n"

                      "\\end{tabular}}\n";
        // todo unify
        summary << "\n\n\\newpage\n"
                   "\\section*{Summary}\n"
                   "\\resizebox{\\textwidth}{!}{\n"
                   "\\begin{tabular}{llll}\n"
                   "\\toprule\n"

                   "Runtime (total): & \\textbf{" << time << "s}\\\\\n"
                   "Runtime (per unique equation): & \\textbf{" << time/stats.eqnUniqueTotal << "s}\\\\\n"
                   "Time to factor $N$: & \\textbf{" << (time/stats.eqnUniqueTotal * (n+2)) << "s}\\\\\\midrule[0.05pt]\n"

                   "Rounds (total): & $" << stats.roundsTotal << "$&\n"
                   "Unique Equations found (duplicates) & $" << stats.eqnUniqueTotal << "$ ($" << stats.eqnDuplicates  << "$)\\\\\n"
                   "Rounds without distance reduction: & " << stats.roundsWithoutReduction << " (" << round(1.0 * stats.roundsWithoutReduction/stats.roundsTotal * 10000) /100.0 << "\\%)&\n"
                   "Unique equations per round & $" << 1.0 * stats.eqnUniqueTotal/stats.roundsTotal << "$\\\\\n"
                   "Rounds without delayed stages: & " << stats.roundsWithoutDelayedStages << " (" << round(1.0 * stats.roundsWithoutDelayedStages/stats.roundsTotal * 10000) /100.0 << "\\%)&\n"
                   "Unique equations per round with 1+ eqn. & $" << stats.avgNumUniqEqnPerRoundWithEqn << "$\\\\\n"
                   "Rounds with equations: & " << stats.roundsWithEquations << " (" << round(1.0 * stats.roundsWithEquations/stats.roundsTotal * 10000) /100.0 << "\\%)&\n"
                   "Total equations per round with 1+ eqn.   & $" << stats.avgNumEqnPerRoundWithEqn << "$\\\\\\midrule[0.05pt]\n"

                   "Total stages checked for equ. (rounds w. Equ.): & " << stats.totalStagesCheckedForEquationsWithEquations << "& Total stages checked for equ. (rounds wo. Equ.): & " << stats.totalStagesCheckedForEquationsWithoutEquations << "\\\\\n"
                   "Min. stages checked for equ. (rounds w. Equ.): & " << stats.minStagesCheckedForEquationsWithEquations << "& Min. stages checked for equ. (rounds wo. Equ.): & " << stats.minStagesCheckedForEquationsWithoutEquations << "\\\\\n"
                   "Max. stages checked for equ. (rounds w. Equ.): & " << stats.maxStagesCheckedForEquationsWithEquations << "& Max. stages checked for equ. (rounds wo. Equ.): & " << stats.maxStagesCheckedForEquationsWithoutEquations << "\\\\\n"
                   "Avg. stages checked for equ. (rounds w. Equ.): & " << round(stats.avgStagesCheckedForEquationsWithEquations * 100.0)/100.0 << "& Avg. stages checked for equ. (rounds wo. Equ.): & " << round(stats.avgStagesCheckedForEquationsWithoutEquations * 100.0)/100.0 << "\\\\\\midrule[0.05pt]\n"

                   "Min. time NewEnum (total): & " << stats.minNewEnum << "s&\n"
                   "Min. time NewEnum (w.eqn.): & " << stats.minNewEnumWE << "s\\\\\n"
                   "Max. time NewEnum: & " << stats.maxNewEnum << "s&\n"
                   "Max. time NewEnum: & " << stats.maxNewEnumWE << "s\\\\\n"
                   "Avg. time NewEnum: & " << stats.avgNewEnum << "s&\n"
                   "Avg. time NewEnum: & " << stats.avgNewEnumWE << "s\\\\\\midrule[0.05pt]\n"

                   "Min. time slight BKZ: & " << stats.minSlightBkz << "s&\n"
                   "Max. distance reduction: & " << stats.maxDistanceReduction << "\\\\\n"
                   "Max. time slight BKZ: & " << stats.maxSlightBkz << "s&\n"
                   "Min. distance reduction: & " << stats.minDistanceReduction << "\\\\\n"
                   "Avg. time slight BKZ: & " << stats.avgSlightBkz << "s&\n"
                   "Avg. distance reduction: & " << stats.avgDistanceReduction << "\\\\\\midrule[0.05pt]\n"

                   "Equations from NewEnum only: &" << eqnList.size() << " (" << round(eqnList.size() * 1000.0 / stats.eqnUniqueTotal) / 10.0
                << "\\%) & Equations from CF: & " << eqnList_cf.size() << " (" << round(eqnList_cf.size() * 1000.0 / stats.eqnUniqueTotal)/10.0 << "\\%)\\\\\n"
                   "Min. $v$ value (without CF): & " << v_min << " & Min. $v$ value (CF only): & " << v_cf_min << "\\\\\n"
                   "Max. $v$ value (without CF): & " << v_max << " & Max. $v$ value (CF only): & " << v_cf_max << "\\\\\n"
                   "Avg. $v$ value (without CF): & " << (round(v_avg * 100.0)/100.0) << " & Avg. $v$ value (CF only): & " << (round(v_cf_avg * 100.0)/100.0) << "\\\\\n"
                   "Median $v$ value (without CF): & " << v_median << " & Median $v$ value (CF only): & " << v_cf_median << "\\\\\\bottomrule\n"

                   "\\end{tabular}}\n";
    }

    void texToPdf()
    {
        if(chdir("output") != 0)
        {
            cerr << "Couldn't change directory to output." << endl;
            return;
        }

        // run pdf latex twice on summary and equation list.
        string pdftex = "pdflatex " + summaryName + ".tex && pdflatex " + statsName + ".tex";
        if(system(pdftex.c_str()) != 0 || system(pdftex.c_str()) != 0)
        {
            cerr << "Couldn't run pdflatex." << endl;
        }

        string trash;
        trash = summaryName + ".log";
        remove(trash.c_str());
        trash = summaryName + ".aux";
        remove(trash.c_str());
        trash = statsName + ".log";
        remove(trash.c_str());
        trash = statsName + ".aux";
        remove(trash.c_str());

        if(chdir("..") != 0)
            cerr << "Couldn't go out of output directory." << endl;
    }

};

#endif    /* FILEOUTPUT_H */