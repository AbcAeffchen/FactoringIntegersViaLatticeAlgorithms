
#include "FileOutput.h"

using namespace NTL;
using namespace std;

void FileOutput::writeEqnFormatted(const Equation& eqn, const Vec<long>& primes)
{
    this->equations << eqn.round << " & ";
    this->equations << eqn.counter << " & ";

    if(eqn.fromContinuedFraction)
        this->equations << "\\checkmark";
    else
        this->equations << "\\cross";

    bool nEmp = false;

    this->equations << "& $";

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
    this->statistics << "$ & $" <<  eqn.time << "$s & ";

    if(eqn.fromContinuedFraction)
        this->statistics << "\\checkmark";
    else
        this->statistics << "\\cross";

    this->statistics << "\\\\" << endl;
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
                     << "\\usepackage[a4paper, left=1cm, right=1cm, top=2cm,bottom=2cm]{geometry}" << endl
                     << "\\usepackage{fancyhdr}" << endl
                     << "\\usepackage{amsmath}" << endl
                     << "\\usepackage{longtable}" << endl
                     << "\\usepackage{setspace}" << endl
                     << "\\usepackage{booktabs}" << endl
                     << "\\usepackage{amssymb}" << endl
                     << "\\usepackage{pifont}" << endl
                     << "\\DeclareMathOperator{\\sign}{\\text{sign}}" << endl
                     << "\\def\\qqquad{\\qquad\\qquad}" << endl
                     << "\\renewcommand{\\checkmark}{\\text{\\ding{51}}}" << endl
                     << "\\newcommand{\\cross}{\\text{\\ding{55}}}" << endl
                     << "\\pagestyle{empty}" << endl << endl
                     << "\\begin{document}" << endl
                     << "\\onehalfspacing" << endl << endl;

    this->statistics << "\\documentclass[a4paper,twoside,10pt]{report}" << endl  << endl
                     << "\\usepackage[a4paper, left=1cm, right=1cm, top=2cm,bottom=2cm]{geometry}" << endl
                     << "\\usepackage{fancyhdr}" << endl
                     << "\\usepackage{amsmath}" << endl
                     << "\\usepackage{longtable}" << endl
                     << "\\usepackage{setspace}" << endl
                     << "\\usepackage{booktabs}" << endl
                     << "\\usepackage{graphicx}" << endl
                     << "\\usepackage{xcolor}" << endl
                     << "\\usepackage{array}" << endl
                     << "\\usepackage{amssymb}" << endl
                     << "\\usepackage{pifont}" << endl
                     << "\\DeclareMathOperator{\\sign}{\\text{sign}}" << endl
                     << "\\def\\qqquad{\\qquad\\qquad}" << endl
                     << "\\pagestyle{empty}" << endl << endl
                     << "\\newcolumntype{C}[1]{>{\\centering\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}" << endl
                     << "\\newcolumntype{R}[1]{>{\\raggedleft\\let\\newline\\\\\\arraybackslash\\hspace{0pt}}m{#1}}" << endl
                     << "\\renewcommand{\\checkmark}{\\text{\\ding{51}}}" << endl
                     << "\\newcommand{\\cross}{\\text{\\ding{55}}}" << endl
                     << "\\newcommand{\\grey}[1]{\\textcolor{black!40}{#1}}" << endl
                     << "\\begin{document}" << endl
                     << "\\onehalfspacing" << endl << endl;
}

void FileOutput::statisticNewRound(long round)
{
    this->statistics << "\\section*{Round " << round << "}" << endl
                     << "\\begin{longtable}{p{13cm}p{4cm}}" << endl;
}

void FileOutput::statisticSlightBKZ(double slightBkz, double newEnum)
{
    this->statistics << "& \\textbf{Distances:}\\newline" << endl
                     << "\\resizebox{\\linewidth}{!}{%" << endl
                     << "\\begin{tabular}[t]{ll}" << endl
                     << "Slight BKZ:&" << slightBkz << "s\\\\" << endl
                     << "NewEnum:&" << newEnum << "s\\\\\\\\" << endl;

}

void FileOutput::statisticsWriteStagesChecked(unsigned long long stagesChecked)
{
    this->statistics << "Checked Stages: & " << stagesChecked << "\\\\\\\\" << endl;
}

void FileOutput::statisticsDistances(RR theoretical, RR heuristic, RR reduced)
{
    this->statistics
                     << "$A=\\frac{1}{4}\\sum_{i=1}^n r_{ii}^2$:&" << trunc(theoretical) << "\\\\" << endl
                     << "$\\frac{1}{5}A$:&" << trunc(heuristic) << "\\\\" << endl
                     << "Reduced to $A'$:&" << trunc(reduced) << "\\\\" << endl
                     << "$A'/A$:&" << conv<double>(reduced/theoretical) << endl
                     << "\\end{tabular}}\\end{longtable}" << endl;
}

void FileOutput::statisticsDelayedStagesOnLevel(int max_level, const vector<vector<double>> &alpha_2_min,
                                                const vector<vector<vector<unsigned long long>>> &delayedAndPerformedStages,
                                                const vector<vector<vector<unsigned long long>>> &delayedStages,
                                                unsigned long long totalDelayedAndPerformedStages)
{
    vector<vector<long long>> sumDelayedStagesOverAlpha2 = vector<vector<long long>>(3,vector<long long>(max_level- 10,0));
    vector<vector<long long>> sumDelayedStagesOverT = vector<vector<long long>>(3,vector<long long>(max_level- 10,0));
    vector<long long> totalSumDelayedStages = vector<long long>(max_level- 10,0);
    long long totalDelayedStages = 0;
    vector<vector<long long>> sumDelayedAndPerformedStagesOverAlpha2 = vector<vector<long long>>(3,vector<long long>(max_level- 10,0));
    vector<vector<long long>> sumDelayedAndPerformedStagesOverT = vector<vector<long long>>(3,vector<long long>(max_level- 10,0));
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

    this->statistics << "\\textbf{Delayed Stages}\\newline\\resizebox{\\linewidth}{!}{%" << endl
                     << "\\begin{tabular}[t]{c|C{4cm}C{4cm}C{4cm}|C{4cm}}" << endl
                     << "\\textbf{" << totalDelayedAndPerformedStages << "}, \\textbf{\\grey{" << totalDelayedStages <<"}}& " << this->t_indicator[0] << "&" << this->t_indicator[1] << "&"  << this->t_indicator[2] << "&\\\\\\midrule" << endl;

    for(long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
    {
        this->statistics << this->alpha_2_indicator_text[alpha_2_indicator];
        for(long t_indicator = 0; t_indicator < 3; t_indicator++)
        {
            this->statistics << "&";
            this->statistics << delayedAndPerformedStages[alpha_2_indicator][t_indicator][0];
            for(long level = 1; level < max_level- 10; level++)
            {
                this->statistics << "," << delayedAndPerformedStages[alpha_2_indicator][t_indicator][level];
            }

            if(alpha_2_min[alpha_2_indicator][t_indicator] > 1)
                this->statistics << " (-)";
            else
                this->statistics << " (" << alpha_2_min[alpha_2_indicator][t_indicator] << ")";

            this->statistics << "\\newline\\grey{";

            this->statistics << delayedStages[alpha_2_indicator][t_indicator][0];
            for(long level = 1; level < max_level- 10; level++)
            {
                this->statistics << "," << delayedStages[alpha_2_indicator][t_indicator][level];
            }
            this->statistics << "}";
        }

        this->statistics << "&";
        this->statistics << sumDelayedAndPerformedStagesOverT[alpha_2_indicator][0];
        for(long level = 1; level < max_level- 10; level++)
        {
            this->statistics << "," << sumDelayedAndPerformedStagesOverT[alpha_2_indicator][level];
        }

        this->statistics << "\\newline\\grey{" << sumDelayedStagesOverT[alpha_2_indicator][0];
        for(long level = 1; level < max_level- 10; level++)
        {
            this->statistics << "," << sumDelayedStagesOverT[alpha_2_indicator][level];
        }
        this->statistics << "}\\\\" << endl;
    }

    this->statistics << "\\midrule" << endl;

    for(long t_indicator = 0; t_indicator < 3; t_indicator++)
    {
        this->statistics << "&" << sumDelayedAndPerformedStagesOverAlpha2[t_indicator][0];
        for(long level = 1; level < max_level- 10;level++)
        {
            this->statistics << "," << sumDelayedAndPerformedStagesOverAlpha2[t_indicator][level];
        }
        this->statistics << "\\newline\\grey{" << sumDelayedStagesOverAlpha2[t_indicator][0];
        for(long level = 1; level < max_level- 10;level++)
        {
            this->statistics << "," << sumDelayedStagesOverAlpha2[t_indicator][level];
        }
        this->statistics << "}";
    }

    this->statistics << "&" << totalSumDelayedAndPerformedStages[0];
    for(long level = 1; level < max_level- 10; level++)
    {
        this->statistics << "," << totalSumDelayedAndPerformedStages[level];
    }
    this->statistics << "\\newline\\grey{" << totalSumDelayedStages[0];
    for(long level = 1; level < max_level- 10; level++)
    {
        this->statistics << "," << totalSumDelayedStages[level];
    }

    this->statistics << "}\\end{tabular}}" << endl;
}

void FileOutput::statisticsNewEquations(const list<Equation>& eqns, const Vec<long>& primes)
{
    if(eqns.size() == 0)
    {
        this->statistics << "\\subsection*{No Equations found}" << endl;
        return;
    }

    this->statistics << "\\subsection*{" << eqns.size() << " new equations}" << endl;
    this->statistics << "\\begin{longtable}{p{1.5cm}lp{5.0cm}R{3.5cm}p{3.0cm}rc}" << endl
                     << "\\toprule" << endl
                     << "dist& level & $u$ & $v$ & $|u-vN|$ & time until & CF\\\\\\midrule" << endl
                     << "\\endfirsthead" << endl
                     << "\\toprule" << endl
                     << "dist& level & $u$ & $v$ & $|u-vN|$ & time until & CF\\\\\\midrule" << endl
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

void FileOutput::writeSettings(const FactoringSettings &settings, long max_prime, long long int seed)
{
    this->equations << "\\section*{Settings}" << endl
                    << "\\begin{tabular}{ll}" << endl
                    << "NTL Version: & " << NTL_VERSION << "\\\\"
                #ifdef __GNUC__
                    << "GCC: &" << __GNUC__ << "." <<__GNUC_MINOR__ << "."  << __GNUC_PATCHLEVEL__ << "\\\\" << endl
                #endif
                    << "$N$ & $" << settings.N << "\\approx 10^{" << round(log(settings.N)/log(10)) << "}$ (ca. " << round(log(settings.N)/log(2)) << " Bits)\\\\"
                    << "$c$ & $" << settings.c << "$\\\\" << endl
                    << "$n$ & $ " << settings.n << "$\\\\" << endl
                    << "$p_n$ & $ " << max_prime << "$\\\\" << endl
                    << "Pruning Level: & $" << settings.max_level << "$\\\\" << endl
                    << "Basis accuracy factor: & $" << settings.accuracy_factor << "$\\\\" << endl
                    << "Reduce ratio: & $" << settings.reduce_ratio << "$\\\\" << endl
                    << "Start reduction factor: & $" << settings.A_start_factor << "$\\\\" << endl
                    << "Strong BKZ block size: & $" << settings.strong_bkz << "$\\\\" << endl
                    << "Slight BKZ block size: & $" << settings.slight_bkz << "$\\\\" << endl
                    << "Seed for random number generator: & $" << seed << "$\\\\" << endl
                    << "Scaling type: & ";
    switch(settings.scalingType)
    {
        case 0: this->equations << "Use various scaleing types:\\\\"
            << "& (\\phantom{7}4\\% of rounds): Scale every row with propability 1/4\\\\"
            << "& (\\phantom{7}4\\% of rounds): Scale every row with propability 3/4\\\\"
            << "& (\\phantom{7}8\\% of rounds): Scale the first n/2 rows with propability 1/4,\\\\& \\phantom{(76\\% of rounds): } the n/2 last rows with prpability 1/2\\\\"
            << "& (\\phantom{7}8\\% of rounds): Scale the first n/2 rows with propability 1/2,\\\\& \\phantom{(76\\% of rounds): } the n/2 last rows with prpability 1/4\\\\"
            << "& (76\\% of rounds): Scale every row with propability 1/2\\\\";
            break;
        case 2: this->equations << "Scale every row with propability 1/4\\\\";
            break;
        case 3: this->equations << "Scale every row with propability 3/4\\\\";
            break;
        case 4: this->equations << "Scale the first n/2 rows with propability 1/4,\\\\& the n/2 last rows with prpability 1/2\\\\";
            break;
        case 5: this->equations << "Scale the first n/2 rows with propability 1/2,\\\\& the n/2 last rows with prpability 1/4\\\\";
            break;
        default: this->equations << "Scale every row with propability 1/2\\\\";
    }

    this->equations << endl << "Use Continued Fractions (CF): & " << (settings.useContinuedFractions ? "yes" : "no") << "\\\\" << endl;

    this->equations << "\\end{tabular}\\\\" << endl;

    this->statistics << "\\section*{Settings}" << endl
                     << "\\begin{tabular}{ll}" << endl
                     << "NTL Version: & " << NTL_VERSION << "\\\\"
                #ifdef __GNUC__
                     << "GCC: &" << __GNUC__ << "." <<__GNUC_MINOR__ << "."  << __GNUC_PATCHLEVEL__ << "\\\\" << endl
                #endif
                     << "$N$ & $" << settings.N << "\\approx 10^{" << round(log(settings.N)/log(10)) << "}$ (ca. " << round(log(settings.N)/log(2)) << " Bits)\\\\"
                     << "$c$ & $" << settings.c << "$\\\\" << endl
                     << "$n$ & $ " << settings.n << "$\\\\" << endl
                     << "$p_n$ & $ " << max_prime << "$\\\\" << endl
                     << "Pruning Level: & $" << settings.max_level << "$\\\\" << endl
                     << "Basis accuracy factor: & $" << settings.accuracy_factor << "$\\\\" << endl
                     << "Reduce ratio: & $" << settings.reduce_ratio << "$\\\\" << endl
                     << "Start reduction factor: & $" << settings.A_start_factor << "$\\\\" << endl
                     << "Strong BKZ block size: & $" << settings.strong_bkz << "$\\\\" << endl
                     << "Slight BKZ block size: & $" << settings.slight_bkz << "$\\\\" << endl
                     << "Seed for random number generator: & $" << seed << "$\\\\" << endl
                     << "Scaling type: & ";
    switch(settings.scalingType)
    {
        case 0: this->statistics << "Use various scaleing types:\\\\"
                << "& (\\phantom{7}4\\% of rounds): Scale every row with propability 1/4\\\\"
                << "& (\\phantom{7}4\\% of rounds): Scale every row with propability 3/4\\\\"
                << "& (\\phantom{7}8\\% of rounds): Scale the first n/2 rows with propability 1/4,\\\\& \\phantom{(76\\% of rounds): } the n/2 last rows with prpability 1/2\\\\"
                << "& (\\phantom{7}8\\% of rounds): Scale the first n/2 rows with propability 1/2,\\\\& \\phantom{(76\\% of rounds): } the n/2 last rows with prpability 1/4\\\\"
                << "& (76\\% of rounds): Scale every row with propability 1/2\\\\";
            break;
        case 2: this->statistics << "Scale every row with propability 1/4\\\\";
            break;
        case 3: this->statistics << "Scale every row with propability 3/4\\\\";
            break;
        case 4: this->statistics << "Scale the first n/2 rows with propability 1/4,\\\\& the n/2 last rows with prpability 1/2\\\\";
            break;
        case 5: this->statistics << "Scale the first n/2 rows with propability 1/2,\\\\& the n/2 last rows with prpability 1/4\\\\";
            break;
        default: this->statistics << "Scale every row with propability 1/2\\\\";
    }

    this->statistics << endl << "Use Continued Fractions (CF): & " << (settings.useContinuedFractions ? "yes" : "no") << "\\\\" << endl;

    this->statistics << "\\end{tabular}\\\\" << endl;

    // todo add explanation for tables to statistics file
}

void FileOutput::statisticsStrongBkzTime(double time)
{
    this->statistics << "\\paragraph{Strong BKZ:}" << time << "s" << endl;
}

void FileOutput::prepareEquationTable()
{
    // start equation table
    this->equations << "%\\begin{landscape}" << endl
                    << "\\begin{longtable}{rrcp{5cm}rp{5cm}r}" << endl
                    << "\\toprule" << endl
                    << "round & \\# & CF & $u$ & $v$ & $|u-vN|$ & $\\sign(u-vN)$\\\\\\midrule" << endl
                    << "\\endfirsthead" << endl
                    << "\\toprule" << endl
                    << "round & \\# & CF & $u$ & $v$ & $|u-vN|$ & $\\sign(u-vN)$\\\\\\midrule" << endl
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

void FileOutput::writeSummary(const Statistics& stats, double time, long n, std::set<Equation> &uniqueEquations)
{
    list<Equation> eqns(uniqueEquations.begin(),uniqueEquations.end()),
                   eqnList,
                   eqnList_cf;

    RR v_avg = conv<RR>(0), v_cf_avg = conv<RR>(0);
    ZZ v_min = conv<ZZ>(0), v_cf_min = conv<ZZ>(0), v_max = conv<ZZ>(0), v_cf_max = conv<ZZ>(0), v_median = conv<ZZ>(0), v_cf_median = conv<ZZ>(0);
    for(std::list<Equation>::iterator it = eqns.begin(); it != eqns.end(); ++it)
    {
        if(it->fromContinuedFraction)
        {
            v_cf_avg += conv<RR>(it->v);
            eqnList_cf.push_front(*it);
        }
        else
        {
            eqnList.push_front(*it);
            v_avg += conv<RR>(it->v);
        }
    }

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

    if(eqnList.size() % 2 == 1)
    {
        std::list<Equation>::iterator it = eqnList.begin();
        std::advance(it,(eqnList.size()-1)/2);
        v_median = it->v;
    }
    else
    {
        std::list<Equation>::iterator it = eqnList.begin();
        std::advance(it,eqnList.size()/2);
        v_median = it->v;
        v_median += (--it)->v;
        v_median /= 2;
    }

    if(!eqnList_cf.empty())
    {
        if (eqnList_cf.size() % 2 == 1)
        {
            std::list<Equation>::iterator it = eqnList_cf.begin();
            std::advance(it, (eqnList_cf.size() - 1) / 2);
            v_cf_median = it->v;
        }
        else
        {
            std::list<Equation>::iterator it = eqnList_cf.begin();
            std::advance(it, eqnList_cf.size() / 2);
            v_cf_median = it->v;
            v_cf_median += (--it)->v;
            v_cf_median /= 2;
        }
    }

    this->statistics << endl << endl
                    << "\\section*{Summary}" << endl
                    << "\\resizebox{\\textwidth}{!}{" << endl
                    << "\\begin{tabular}{llll}" << endl
                    << "\\toprule" << endl

                    << "Runtime (total): & \\textbf{" << time << "s}\\\\" << endl
                    << "Runtime (per unique equation): & \\textbf{" << time/stats.eqnUniqueTotal << "s}\\\\" << endl
                    << "Time to factor $N$: & \\textbf{" << (time/stats.eqnUniqueTotal * (n+1)) << "s}\\\\\\midrule[0.05pt]" << endl

                    << "Rounds (total): & $" << stats.roundsTotal << "$&" << endl
                    << "Unique Equations found (duplicates) & $" << stats.eqnUniqueTotal << "$ ($" << stats.eqnDuplicates  << "$)\\\\" << endl
                    << "Rounds without distance reduction: & " << stats.roundsWithoutReduction << " (" << round(1.0 * stats.roundsWithoutReduction/stats.roundsTotal * 10000) /100.0 << "\\%)&" << endl
                    << "Unique equations per round & $" << 1.0 * stats.eqnUniqueTotal/stats.roundsTotal << "$\\\\" << endl
                    << "Rounds without delayed stages: & " << stats.roundsWithoutDelayedStages << " (" << round(1.0 * stats.roundsWithoutDelayedStages/stats.roundsTotal * 10000) /100.0 << "\\%)&" << endl
                    << "Unique equations per round with 1+ eqn. & $" << stats.avgNumUniqEqnPerRoundWithEqn << "$\\\\" << endl
                    << "Rounds with equations: & " << stats.roundsWithEquations << " (" << round(1.0 * stats.roundsWithEquations/stats.roundsTotal * 10000) /100.0 << "\\%)&" << endl
                    << "Total equations per round with 1+ eqn.   & $" << stats.avgNumEqnPerRoundWithEqn << "$\\\\\\midrule[0.05pt]" << endl

                    << "Total stages checked for equ. (rounds w. Equ.): & " << stats.totalStagesCheckedForEquationsWithEquations << "& Total stages checked for equ. (rounds wo. Equ.): & " << stats.totalStagesCheckedForEquationsWithoutEquations << "\\\\" << endl
                    << "Min. stages checked for equ. (rounds w. Equ.): & " << stats.minStagesCheckedForEquationsWithEquations << "& Min. stages checked for equ. (rounds wo. Equ.): & " << stats.minStagesCheckedForEquationsWithoutEquations << "\\\\" << endl
                    << "Max. stages checked for equ. (rounds w. Equ.): & " << stats.maxStagesCheckedForEquationsWithEquations << "& Max. stages checked for equ. (rounds wo. Equ.): & " << stats.maxStagesCheckedForEquationsWithoutEquations << "\\\\" << endl
                    << "Avg. stages checked for equ. (rounds w. Equ.): & " << round(stats.avgStagesCheckedForEquationsWithEquations * 100.0)/100.0 << "& Avg. stages checked for equ. (rounds wo. Equ.): & " << round(stats.avgStagesCheckedForEquationsWithoutEquations * 100.0)/100.0 << "\\\\\\midrule[0.05pt]" << endl

                    << "Min. time NewEnum (total): & " << stats.minNewEnum << "s&" << endl
                    << "Min. time NewEnum (w.eqn.): & " << stats.minNewEnumWE << "s\\\\" << endl
                    << "Max. time NewEnum: & " << stats.maxNewEnum << "s&" << endl
                    << "Max. time NewEnum: & " << stats.maxNewEnumWE << "s\\\\" << endl
                    << "Avg. time NewEnum: & " << stats.avgNewEnum << "s&" << endl
                    << "Avg. time NewEnum: & " << stats.avgNewEnumWE << "s\\\\\\midrule[0.05pt]" << endl

                    << "Min. time slight BKZ: & " << stats.minSlightBkz << "s&" << endl
                    << "Max. distance reduction: & " << stats.maxDistanceReduction << "\\\\" << endl
                    << "Max. time slight BKZ: & " << stats.maxSlightBkz << "s&" << endl
                    << "Min. distance reduction: & " << stats.minDistanceReduction << "\\\\" << endl
                    << "Avg. time slight BKZ: & " << stats.avgSlightBkz << "s&" << endl
                    << "Avg. distance reduction: & " << stats.avgDistanceReduction << "\\\\\\midrule[0.05pt]" << endl

                    << "Equations from NewEnum only: &" << eqnList.size() << " (" << round(eqnList.size() * 1000.0 / stats.eqnUniqueTotal) / 10.0
                    << "\\%) & Equations from CF: & " << eqnList_cf.size() << " (" << round(eqnList_cf.size() * 1000.0 / stats.eqnUniqueTotal)/10.0 << "\\%)\\\\" << endl
                    << "Min. $v$ value (without CF): & " << v_min << " & Min. $v$ value (CF only): & " << v_cf_min << "\\\\" << endl
                    << "Max. $v$ value (without CF): & " << v_max << " & Max. $v$ value (CF only): & " << v_cf_max << "\\\\" << endl
                    << "Avg. $v$ value (without CF): & " << (round(v_avg * 100.0)/100.0) << " & Avg. $v$ value (CF only): & " << (round(v_cf_avg * 100.0)/100.0) << "\\\\" << endl
                    << "Median $v$ value (without CF): & " << v_median << " & Median $v$ value (CF only): & " << v_cf_median << "\\\\\\bottomrule" << endl

                    <<"\\end{tabular}}" << endl;
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

void FileOutput::statisticsWriteScaledPrimes(const vector<bool> &scaledPrimes,
                                             const Vec<long> &primes)
{
    this->statistics << "\\paragraph*{Scaled prime numbers:}";
    long n = primes.length();
    for(long i = 0; i < n; i++)
    {
        if(scaledPrimes[i])
            this->statistics << primes[i] << " ";
    }

}
