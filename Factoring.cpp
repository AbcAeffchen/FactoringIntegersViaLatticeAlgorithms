
#include "Factoring.h"

using namespace std;
using namespace NTL;

// first 300 prime numbers
const vector<long> _primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
                              61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
                              131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
                              197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269,
                              271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
                              353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431,
                              433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
                              509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599,
                              601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673,
                              677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761,
                              769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857,
                              859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
                              953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031,
                              1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097,
                              1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187,
                              1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277,
                              1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327,
                              1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439,
                              1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499,
                              1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583,
                              1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663,
                              1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747,
                              1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847,
                              1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931,
                              1933, 1949, 1951, 1973, 1979, 1987};

void Factoring::randomScale1(Mat<ZZ> &basis)
{
    uniform_int_distribution<int> dist(0,1);   // [0,1] random numbers uniform distribution

    for(long i = 1; i <= basis.NumCols(); i++)
    {
        if(dist(this->rgen) == 1)       // P(scale row) = 1/2
        {
            basis(i) *= 2;     // multiply the whole row by 2
        }
    }
}

void Factoring::randomScale2(Mat<ZZ> &basis)
{
    uniform_int_distribution<int> dist(0,3);   // [0,3] random numbers uniform distribution

    for(long i = 1; i <= basis.NumCols(); i++)
    {
        if(dist(this->rgen) == 1)       // P(scale row) = 1/4
        {
            basis(i) *= 2;     // multiply the whole row by 2
        }
    }
}

void Factoring::setScaledAndReducedBasis()
{
    // scale
    this->B_scaled_transposed = this->B;    // not yet transposed

    switch(this->settings.scalingType)
    {
        case 2:
            this->randomScale2(this->B_scaled_transposed);
            break;
        default:
            this->randomScale1(this->B_scaled_transposed);
            break;
    }

    // reduce
    transpose(this->B_scaled_transposed, this->B_scaled_transposed);    // now it is transposed
    BKZ_FP(this->B_scaled_transposed, this->U_scaled, 0.99, this->settings.slight_bkz);

    transpose(this->U_scaled,this->U_scaled);
    inv(this->U_scaled_inv, this->U_scaled);

    // apply the reduction to the coordinate vector
    mul(this->target_scaled_coordinates, conv<Mat<RR>>(this->U_scaled_inv), this->target_coordinates);

    return;
}

void Factoring::setBasis(long accuracy_factor)
{
    cout << "Setting Prime Lattice Basis: ";
    long n = this->primes.length();

    this->B.SetDims(n, n + 1);    // transposed
    // Setting the basis
    for(long i = 1; i <= n; i++)
    {
        this->B(i,i) = conv<ZZ>(accuracy_factor * SqrRoot(log(conv<RR>(this->primes(i)))));
        if(this->c < 1)
            this->B(i, n + 1) = conv<ZZ>(pow(conv<RR>(this->N),this->c) * accuracy_factor * log(conv<RR>(this->primes(i))));        // todo N^c can be precomputed
        else
            this->B(i, n + 1) = conv<ZZ>(conv<RR>(this->N * accuracy_factor) * log(conv<RR>(this->primes(i))));
    }

    cout << "finished" << endl;
}

Vec<RR> Factoring::orthogonalProjection_pl()
{
    Vec<RR> orth_target;
    long n = this->primes.length();
    orth_target.SetLength(n);

    RR N_RR = conv<RR>(this->N);
    RR N_RR_2c = pow(N_RR, 2 * this->c);

    double log_prim_prod = 0;
    for(long i = 1; i <= n; i++)
        log_prim_prod += log(this->primes(i));

    RR a = ( N_RR_2c * log(N_RR) ) / (1 + N_RR_2c * log_prim_prod);

    for(long i = 1; i <= n; i++)
        orth_target(i) = a;

    return orth_target;
}

void Factoring::setTargetCoordinates()
{
    // getting the coordinates respecting the prime number lattice
    Vec<RR> target_coordinates = this->orthogonalProjection_pl();

    // taking care of the strong reduction
    mul(this->target_coordinates, conv<Mat<RR> >(this->U_inv), target_coordinates);

    // make the shift
    long n = this->primes.length();
    this->shift.SetLength(n);

    for(long i = 1; i <= n; i++)
    {
        // the shift could also be done by using floor or ceil.
        this->shift(i) = round(this->target_coordinates(i));
        this->target_coordinates(i) -= this->shift(i);
    }
}

void Factoring::reduceBasis(long strong_bkz)
{
    this->timer.startTimer();

    // this->B is still transposed

    cout << "Strong BKZ reduction: ";
    BKZ_FP(this->B, this->U, 0.99, strong_bkz);                 // strong reducing

    transpose(this->B,this->B);
    transpose(this->U,this->U);

    inv(this->U_inv, this->U);

    cout << "finished" << endl;

    this->file.statisticsStrongBkzTime(this->timer.stopTimer());
}

void Factoring::search()
{
    list<Equation> newEquations;
    double newEnumTime;
    double slightBkzTime;
    long newDuplicates;
    unsigned long round = 0;

    NewEnum newEnum = NewEnum(this->settings, this->timer, this->file, this->stats,
                              this->primes, this->U, this->shift);

    RR theoretical, heuristic, reduced;

    while(this->uniqueEquations.size() < this->settings.min_eqns)
    {
        round++;
        cout << "Round " << round << endl;
        this->file.statisticNewRound(round);

        this->timer.startTimer();
        this->setScaledAndReducedBasis();    // scale and slight bkz
        slightBkzTime = this->timer.stopTimer();

        this->timer.startTimer();

        newEnum.run(round, this->B_scaled_transposed, this->U_scaled, this->target_scaled_coordinates);
        newEquations = newEnum.getEquations();
        newDuplicates = this->addEquations(newEquations);

        newEnumTime = this->timer.stopTimer();

        newEnum.getDistances(theoretical,heuristic,reduced);

        // statistics
        this->stats.updateRoundStats(newEnum.L.totalDelayedAndPerformedStages > 0, newEquations.size() > 0);
        this->stats.updateDistanceStats(theoretical,heuristic,reduced);
        this->file.statisticsDelayedStagesOnLevel(this->settings.max_level,newEnum.L.alpha_2_min,
                                                  newEnum.L.maxDelayedAndPerformedStages,
                                                  newEnum.L.delayedStages,
                                                  newEnum.L.totalDelayedAndPerformedStages);
        this->file.statisticsDistances(theoretical,heuristic,reduced);
        this->file.statisticSlightBKZ(slightBkzTime, newEnumTime);
        this->file.statisticsNewEquations(newEquations,this->primes);
        this->stats.newSlightBkzTime(slightBkzTime);
        this->stats.newNewEnumTime(newEnumTime, newEquations.size() > 0);

        // display round results
        cout << " -> total equations (new | total):       " << newEquations.size() << "  |  " << this->uniqueEquations.size() << " (" << this->settings.min_eqns << ")" << endl;
        cout << " -> duplicate equations (new | total):   " << newDuplicates << "  |  " << this->eqnDuplicates << endl;
    }

    this->stats.closeStatistics(round,this->uniqueEquations.size(),this->eqnDuplicates);

    this->file.writeFormattedEquationList(this->uniqueEquations, this->primes);

    this->file.writeSummary(this->stats, this->timer.getTotalTime(), this->settings.n, this->uniqueEquations);
    this->file.closeEquationFile();
    this->file.closeStatisticsFile();

    this->file.texToPdf();
}

long Factoring::addEquations(list<Equation>& newEquations)
{
    long duplicates = 0;
    for (std::list<Equation>::iterator it=newEquations.begin(); it != newEquations.end(); ++it)
    {
        std::pair<std::set<Equation>::iterator,bool> res = this->uniqueEquations.insert(*it);
        if(!res.second)  // equation was not inserted
        {
            duplicates++;
            // objects in sets cannot be changed. they have to be erased and
            // reinserted after the change, due to sorting reasons.
            Equation upCounted = *(res.first);
            upCounted.counter++;
            this->uniqueEquations.erase(res.first);
            this->uniqueEquations.insert(upCounted);
        }
    }
    this->eqnDuplicates += duplicates;

    return duplicates;
}

void Factoring::setPrimes(long n)
{
    this->primes.SetLength(n);

    for(long i = 0; i < min(n,300); ++i)
        this->primes[i] = _primes[i];

    for(long i = 301; i <= n; i++)
        NextPrime(this->primes(i),this->primes(i-1)+2);
}

Factoring::Factoring(const FactoringSettings &settings)
        : settings(settings), N(settings.N), c(settings.c)
{
    this->timer = Timer();

    this->setPrimes(settings.n);

    // setup random number generator
    long long int seed;
    if(settings.seed_type <= -2)
        seed = time(NULL);
    else if(settings.seed_type == -1)
    {
        random_device rd;
        seed = rd();
    }
    else
        seed = settings.seed_type;

    this->rgen = mt19937(seed);

    // write the current settings
    this->file.writeSettings(this->settings, this->primes(settings.n), seed);

    this->file.prepareEquationTable();
    this->setBasis(settings.accuracy_factor);
    this->reduceBasis(settings.strong_bkz);
    this->setTargetCoordinates();
    this->search();
}

ZZ getN(long e)
{
    ZZ p1 = NextPrime(SqrRoot(power(conv<ZZ>(10), e)));
    ZZ p2 = NextPrime(p1 + 2);
    return p1 * p2;
}