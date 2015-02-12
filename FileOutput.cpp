
#include "FileOutput.h"

using namespace NTL;
using namespace std;

void FileOutput::writeEqnFormatted(const Equation& eqn, const Vec<long>& primes)
{
    this->equations << eqn.round << " & ";
    this->equations << eqn.counter << " & ";

    bool nEmp = false;

    this->equations << "$";

    for(long i = 1; i <= primes.length(); i++)
    {
        if(eqn.e(i) > 0)
        {
            if(nEmp)
            {
                this->equations << " \\cdot ";
            }
            this->equations << primes(i);
            if(eqn.e(i) > 1)
            {
                this->equations << "^{" << eqn.e(i) << "}";
            }
            nEmp = true;
        }
    }
    if(!nEmp)
    {
        this->equations << "1";
    }
    this->equations << "$ & $"<< eqn.v << "$ & $";
    nEmp = false;

    for(long i = 1; i <= primes.length(); i++)
    {
        if(eqn.e(i) < 0)
        {
            if(nEmp)
            {
                this->equations << " \\cdot ";
            }
            this->equations << primes(i);
            if(eqn.e(i) < -1)
            {
                this->equations << "^{" << -eqn.e(i) << "}";
            }
            nEmp = true;
        }
    }
    if(!nEmp)
    {
        this->equations << "1";
    }
    this->equations << "$ & $";
    if(eqn.e(primes.length()+1) == 1)
        this->equations << "-1";
    else
        this->equations << "1";
    this->equations << "$\\\\" << endl;
}

void FileOutput::writeEqnStatistics(const Equation& eqn, const Vec<long>& primes)
{
    bool nEmp = false;
    this->statistics << eqn.reduced << " & "
                     << eqn.level << " & ";

    this->statistics << "$";

    for(long i = 1; i <= primes.length(); i++)
    {
        if(eqn.e(i) > 0)
        {
            if(nEmp)
                this->statistics << " \\cdot ";
            this->statistics << primes(i);
            if(eqn.e(i) > 1)
            {
                this->statistics << "^{" << eqn.e(i) << "}";
            }
            nEmp = true;
        }
    }
    if(!nEmp)
    {
        this->statistics << "1";
    }
    this->statistics << "$ & $"<< eqn.v << "$ & $";
    nEmp = false;

    for(long i = 1; i <= primes.length(); i++)
    {
        if(eqn.e(i) < 0)
        {
            if(nEmp)
            {
                this->statistics << " \\cdot ";
            }
            this->statistics << primes(i);
            if(eqn.e(i) < -1)
            {
                this->statistics << "^{" << -eqn.e(i) << "}";
            }
            nEmp = true;
        }
    }
    if(!nEmp)
    {
        this->statistics << "1";
    }
    this->statistics << "$ & $";
    this->statistics << "$" << eqn.time << "s \\\\" << endl;
}

void FileOutput::createDirectory()
{
    system("mkdir output"); 
}

/**
 * @return Current date/time, format is YYYY-MM-DD.HH:mm:ss
 */
string FileOutput::getFilePrefix()
{
    time_t now = time(0);
    struct tm  tstruct;
    tstruct = *localtime(&now);
    char buf[80];
    strftime(buf, sizeof(buf), "%Y-%m-%d_%H-%M-%S", &tstruct);
    
    return buf;
}

void FileOutput::debugSeparator(string separatorText)
{
    this->debug << endl << endl
                << "+---------------------------------------------------------------------------------------+" << endl
                << "|                                 " << separatorText << "                                   |" << endl
                << "+---------------------------------------------------------------------------------------+" << endl
                << endl << endl;
}

FileOutput::FileOutput()
{
    this->createDirectory();
    
    string filePrefix = this->getFilePrefix();
    std::stringstream formattedName;
    formattedName << "./output/" << filePrefix << "_equations.tex";
    eqnName = filePrefix + "_equations";
    std::stringstream statisticsName;
    statisticsName << "./output/" << filePrefix << "_statistics.tex";
    statsName = filePrefix + "_statistics";
    
    this->debug.open("./output/debug.txt", ios::out | ios::app);
    this->debugSeparator(filePrefix);
    
    this->equations.open(formattedName.str().c_str(), ios::out);
    this->statistics.open(statisticsName.str().c_str(), ios::out);

    // Init formatted file
    this->equations  << "\\documentclass[a4paper,twoside,10pt]{report}" << endl << endl
                     << "\\usepackage[a4paper, left=2cm, right=2cm, top=2cm,bottom=2cm]{geometry}" << endl
                     << "\\usepackage{fancyhdr}" << endl
                     << "\\usepackage{amsmath}" << endl
                     << "\\usepackage{longtable}" << endl
                     << "\\usepackage{setspace}" << endl
                     << "\\usepackage{booktabs}" << endl
                     << "\\pagestyle{empty}" << endl
                     << "\\DeclareMathOperator{\\sign}{\\text{sign}}" << endl
                     << "\\def\\qqquad{\\qquad\\qquad}" << endl
                     << "\\pagestyle{empty}" << endl << endl
                     << "\\begin{document}" << endl
                     << "\\onehalfspacing" << endl << endl;
    this->statistics << "\\documentclass[a4paper,twoside,10pt]{report}" << endl  << endl
                     << "\\usepackage[a4paper, left=2cm, right=2cm, top=2cm,bottom=2cm]{geometry}" << endl
                     << "\\usepackage{fancyhdr}" << endl
                     << "\\usepackage{amsmath}" << endl
                     << "\\usepackage{longtable}" << endl
                     << "\\usepackage{setspace}" << endl
                     << "\\usepackage{booktabs}" << endl
                     << "\\pagestyle{empty}" << endl
                     << "\\DeclareMathOperator{\\sign}{\\text{sign}}" << endl
                     << "\\def\\qqquad{\\qquad\\qquad}" << endl
                     << "\\pagestyle{empty}" << endl << endl
                     << "\\begin{document}" << endl
                     << "\\onehalfspacing" << endl << endl;
}

void FileOutput::statisticNewRound(long round)
{
    this->statistics << "\\section*{Round " << round << "}" << endl
                     << "\\begin{longtable}{p{7.3cm}p{5.3cm}p{4.3cm}}" << endl;        
}

void FileOutput::statisticSlightBKZ(double slightBkz, double newEnum)
{
    this->statistics << "& Slight BKZ: " << slightBkz << "s\\newline " << endl
        << "NewEnum: " << newEnum << "s ";

    this->statistics << "\\end{longtable}" << endl;
}

void FileOutput::statisticsDistances(RR theoretical, RR heuristical, RR reduced)
{
    this->statistics << "\\textbf{Distances:} \\newline " << endl;
    this->statistics << "Theoretical: " << trunc(theoretical) << "\\newline" << endl;        
    this->statistics << "Heuristical: " << trunc(heuristical) << "\\newline" << endl;        
    this->statistics << "Reduced: " << trunc(reduced) << "\\newline" << endl;        
    this->statistics << "Ratio: " << conv<double>(reduced/theoretical) << endl;        
}

void FileOutput::statisticsDelayedStagesOnLevel(vector<long> delayedStagesCounter)
{
    this->statistics << "& \\textbf{Delayed Stages:} " << endl;
    for(unsigned long i = 1; i < delayedStagesCounter.size(); i++)
        this->statistics << "\\newline " << "Lvl " << (i + 10) << ": " << delayedStagesCounter[i] << endl;

    this->statistics << "\\newline " << "Total: " << delayedStagesCounter[0] << endl;
}

void FileOutput::statisticsNewEquations(const list<Equation>& eqns, const Vec<long>& primes)
{
    if(eqns.size() == 0)
    {
        this->statistics << "\\subsection*{No Equations found}" << endl;
        return;
    }        

    this->statistics << "\\subsection*{Equations}" << endl;
    this->statistics << "\\begin{longtable}{llp{6cm}rp{4cm}rcr}" << endl
                    << "\\toprule" << endl
                    << "dist& level & $u$ & $v$ & $|u-vN|$ & time until\\\\\\midrule" << endl
                    << "\\endfirsthead" << endl
                    << "\\toprule" << endl
                    << "dist& level & $u$ & $v$ & $|u-vN|$ & time until\\\\\\midrule" << endl
                    << "\\endhead" << endl;
    for(std::list<Equation>::const_iterator it = eqns.begin(); it != eqns.end(); ++it)
    {
        this->writeEqnStatistics(*it, primes);
    }
    this->statistics << "\\end{longtable}" << endl;
}    

void FileOutput::writeFormattedEquationList(std::set<Equation>& eqns, const Vec<long>& primes)
{
    list<Equation> eqnList(eqns.begin(),eqns.end());

    eqnList.sort(sort_equations);
    for(std::list<Equation>::iterator it = eqnList.begin(); it != eqnList.end(); ++it)
    {
        this->writeEqnFormatted(*it, primes);
    }

}

void FileOutput::writeSettings(ZZ N, RR c, long accuracy_factor, int s_max, double A_start_factor, double reduce_ratio, long strong_bkz, long slight_bkz, long prime_num, long max_prime, long long int seed)
{
    this->equations << "\\begin{tabular}{ll}" << endl
                    << "$N$ & $" << N << "\\approx 10^{" << round(log(N)/log(10)) << "}$ (ca. " << round(log(N)/log(2)) << " Bits)\\\\"
                    << "$c$ & $" << c << "$\\\\" << endl
                    << "$n$ & $ " << prime_num << "$\\\\" << endl
                    << "$p_n$ & $ " << max_prime << "$\\\\" << endl
                    << "Pruning Level: & $" << s_max << "$\\\\" << endl
                    << "Basis accuracy factor: & $" << accuracy_factor << "$\\\\" << endl
                    << "Reduce ratio: & $" << reduce_ratio << "$\\\\" << endl
                    << "Start reduction factor: & $" << A_start_factor << "$\\\\" << endl
                    << "Strong BKZ block size: & $" << strong_bkz << "$\\\\" << endl
                    << "Slight BKZ block size: & $" << slight_bkz << "$\\\\" << endl
                    << "Seed for random number generator: & $" << seed << "$" << endl
                    << "\\end{tabular}\\\\" << endl;
    
    this->statistics << "\\begin{tabular}{ll}" << endl
                     << "$N$ & $" << N << "\\approx 10^{" << round(log(N)/log(10)) << "}$ (ca. " << round(log(N)/log(2)) << " Bits)\\\\"
                     << "$c$ & $" << c << "$\\\\" << endl
                     << "$n$ & $ " << prime_num << "$\\\\" << endl
                     << "$p_n$ & $ " << max_prime << "$\\\\" << endl
                     << "Pruning Level: & $" << s_max << "$\\\\" << endl
                     << "Basis accuracy factor: & $" << accuracy_factor << "$\\\\" << endl
                     << "Reduce ratio: & $" << reduce_ratio << "$\\\\" << endl
                     << "Start reduction factor: & $" << A_start_factor << "$\\\\" << endl
                     << "Strong BKZ block size: & $" << strong_bkz << "$\\\\" << endl
                     << "Slight BKZ block size: & $" << slight_bkz << "$\\\\" << endl
                     << "Seed for random number generator: & $" << seed << "$" << endl
                     << "\\end{tabular}\\\\" << endl;
}

void FileOutput::statisticsStrongBkzTime(double time)
{
    this->statistics << "\\paragraph{Strong BKZ:}" << time << "s" << endl;
}

void FileOutput::prepareEquationTable()
{
    // start equationtable
    this->equations << "%\\begin{landscape}" << endl
                    << "\\begin{longtable}{rrp{6cm}rp{5cm}r}" << endl
                    << "\\toprule" << endl
                    << "round & \\# & $u$ & $v$ & $|u-vN|$ & $\\sign(u-vN)$\\\\\\midrule" << endl
                    << "\\endfirsthead" << endl
                    << "\\toprule" << endl
                    << "round & \\# & $u$ & $v$ & $|u-vN|$ & $\\sign(u-vN)$\\\\\\midrule" << endl
                    << "\\endhead" << endl;
}

void FileOutput::closeEquationFile()
{
    this->equations << "\\end{longtable}" << endl
                    << "%\\end{landscape}" << endl
                    << "\\end{document}" << endl;
    this->equations.close();
}

void FileOutput::closeStatisticsFile()
{
    this->statistics << "\\end{document}" << endl;
    this->statistics.close();
}

void FileOutput::writeSummary(const Statistics& stats, double time)
{
    this->statistics << endl << endl 
                    << "\\section*{Summary}" << endl
                    << "\\begin{tabular}{ll}" << endl
                    << "Rounds (total): & $" << stats.roundsTotal << "$\\\\" << endl
                    << "Rounds without distance reduction: & " << stats.roundsWithoutReduction << " (" << round(1.0 * stats.roundsWithoutReduction/stats.roundsTotal * 10000) /100.0 << "\\%)\\\\" << endl
                    << "Rounds without delayed stages: & " << stats.roundsWithoutDelayedStages << " (" << round(1.0 * stats.roundsWithoutDelayedStages/stats.roundsTotal * 10000) /100.0 << "\\%)\\\\" << endl
                    << "Rounds with equations: & " << stats.roundsWithEquations << " (" << round(1.0 * stats.roundsWithEquations/stats.roundsTotal * 10000) /100.0 << "\\%)\\\\\\midrule" << endl
                    << "Equations found (unique) & $" << stats.eqnUniqueTotal << "$\\\\" << endl
                    << "Duplicate Equations & $" << stats.eqnDuplicates << "$\\\\" << endl
                    << "Unique equations per round & $" << 1.0 * stats.eqnUniqueTotal/stats.roundsTotal << "$\\\\" << endl
                    << "Unique equations per round with 1+ eqn. & $" << stats.avgNumUniqEqnPerRoundWithEqn << "$\\\\" << endl
                    << "Total equations per round with 1+ eqn.   & $" << stats.avgNumEqnPerRoundWithEqn << "$\\\\\\midrule" << endl 
                    << "Runtime (total): & $" << time << "$s\\\\" << endl
                    << "Runtime (per unique equation): & $" << time/stats.eqnUniqueTotal << "$s\\\\\\midrule[0.05pt]" << endl
                    << "Min. time slight BKZ: & " << stats.minSlightBkz << "s\\\\" << endl
                    << "Max. time slight BKZ: & " << stats.maxSlightBkz << "s\\\\" << endl
                    << "Avg. time slight BKZ: & " << stats.avgSlightBkz << "s\\\\\\midrule[0.05pt]" << endl
                    << "Min. time NewEnum (total): & " << stats.minNewEnum << "s\\\\" << endl
                    << "Max. time NewEnum: & " << stats.maxNewEnum << "s\\\\" << endl
                    << "Avg. time NewEnum: & " << stats.avgNewEnum << "s\\\\\\midrule[0.05pt]" << endl
                    << "Min. time NewEnum (w.eqn.): & " << stats.minNewEnumWE << "s\\\\" << endl
                    << "Max. time NewEnum: & " << stats.maxNewEnumWE << "s\\\\" << endl
                    << "Avg. time NewEnum: & " << stats.avgNewEnumWE << "s\\\\\\midrule" << endl
                    << "Max. distance reduction: & " << stats.maxDistanceReduction << "\\\\" << endl
                    << "Min. distance reduction: & " << stats.minDistanceReduction << "\\\\" << endl
                    << "Avg. distance reduction: & " << stats.avgDistanceReduction << "\\\\" << endl
                    <<  "\\end{tabular}" << endl;
}

void FileOutput::texToPdf()
{
    string pdftex;
    chdir("output");
    pdftex = "pdflatex " + this->eqnName + ".tex && pdflatex " + this->statsName + ".tex";
    system(pdftex.c_str());
    system(pdftex.c_str());
    string trash;
    trash = this->eqnName + ".log";
    remove(trash.c_str());
    trash = this->eqnName + ".aux";
    remove(trash.c_str());
    trash = this->statsName + ".log";
    remove(trash.c_str());
    trash = this->statsName + ".aux";
    remove(trash.c_str());
    chdir("..");
}