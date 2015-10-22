
#include "Factoring.h"

using namespace std;
using namespace NTL;

/**
 * scales the lattice basis, BKZ-reduces it slightly, and also converting the
 * target coordinates
 */
void Factoring::randomScale()
{
    // scale
    this->B_scaled = this->B;

    uniform_int_distribution<int> dist(0,1);   // [0,1] random numbers uniform distribution

    for(long i = 1; i <= this->B_scaled.NumCols(); i++)
    {
        if(dist(this->rgen) == 1)       // P(scale row) = 1/2
        {
            this->B_scaled(i) *= 2;     // multiply the whole row by 2
        }
    }
    // reduce
    Mat<ZZ> U_scaled;  // Transition Matrix for the following BKZ-reduction
    this->B_scaled = transpose(this->B_scaled);
    BKZ_FP(this->B_scaled, U_scaled, 0.99, this->slight_bkz);
    this->B_scaled = transpose(this->B_scaled);

    this->U_scaled = transpose(U_scaled);
    inv(this->U_scaled_inv, this->U_scaled);

    // apply the reduction to the coordinate vector
    mul(this->target_scaled_coordinates, conv<Mat<RR> >(this->U_scaled_inv), this->target_coordinates);

    return;
}

/**
 * Generates the prime lattice basis and make a strong BKZ reduction with
 * block size strong_bkz. The entrys of the basis have to be (big) integers
 * so the basis is multiplied by accurancy_factor bevor rounding.
 * @param accuracy_factor
 */
void Factoring::setBasis(long accuracy_factor)
{
    cout << "Setting Prime Lattice Basis" << endl;
    RR c = this->c;
    long n = this->primes.length();

    Mat<ZZ> basis;

    basis.SetDims(n, n + 1);    // transposed
    // Setting the basis
    for(long i = 1; i <= n; i++)
    {
        basis(i,i) = conv<ZZ>(accuracy_factor * SqrRoot(log(conv<RR>(this->primes(i)))));
        if(c < 1)
            basis(i, n + 1) = conv<ZZ>(pow(conv<RR>(this->N),c) * accuracy_factor * log(conv<RR>(this->primes(i))));        // todo N^c can be precomputed
        else
            basis(i, n + 1) = conv<ZZ>(conv<RR>(this->N * accuracy_factor) * log(conv<RR>(this->primes(i))));
    }

    this->B = transpose(basis);
}

/**
 * Computes the coordinates of the target vector projected orthogonally into
 * the lattice plane. This method runs fast since the coordinates of the
 * projection are all the same an can be computed by a formula. NOTICE:
 * This method works only for the standard target vector.
 * @return Coordinates of the projected target vector.
 */
Vec<RR> Factoring::orthogonalProjection_pl()
{
    Vec<RR> orth_target;
    long n = this->primes.length();
    orth_target.SetLength(n);

    RR N_RR = conv<RR>(this->N);
    RR N_RR_2c = pow(N_RR, 2 * this->c);

    RR prime_prod;
    NTL::set(prime_prod); // prim_prod = 1;

    for(long i = 1; i <= n; i++)
        prime_prod *= conv<RR>(this->primes(i));

    RR a = ( N_RR_2c * log(N_RR) ) / (1 + N_RR_2c * log(prime_prod) );

    for(long i = 1; i <= n; i++)
        orth_target(i) = a;

    return orth_target;
}

/**
 * Computes the coordinates of the target vector respecting the strong reduced
 * basis and shifted to the ground mesh. Requires the the transition matrix
 * U_inv to be set.
 */
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

/**
 * Performs a strong BKZ reduction and sets this->U
 * @param strong_bkz Block size of the BKZ reduction
 */
void Factoring::reduceBasis(long strong_bkz)
{
    this->timer.startTimer();
    Mat<ZZ> U,      // Transition matrix
            basis = transpose(this->B);
    cout << "starting strong BKZ" << endl;
    BKZ_FP(basis, U, 0.99, strong_bkz);                 // strong reducing

    this->B = transpose(basis);

    this->U = transpose(U);
    inv(this->U_inv, this->U);

    cout << "finished strong BKZ" << endl;

    this->file.statisticsStrongBkzTime(this->timer.stopTimer());
}
/**
 * Runs the search for relations. This method runs a loop until 3*n relations
 * were found
 */
void Factoring::search()
{
//        NewEnum newEnum;
    list<Equation> newEquations;
    double newEnumTime;
    double slightBkzTime;
    long newDuplicates;
    long round = 0;

    while(this->uniqueEquations.size() < this->min_eqns)
    {
        round++;
        cout << "Round " << round << endl;
        this->file.statisticNewRound(round);

        this->timer.startTimer();
        this->randomScale();    // scale and slight bkz
        slightBkzTime = this->timer.stopTimer();


        this->timer.startTimer();
        NewEnum newEnum = NewEnum(this->timer, this->file, this->stats, round, this->N, this->primes, this->B_scaled,
                          this->U, this->U_scaled, this->target_scaled_coordinates,
                          this->shift, this->s_max, this->restart_ratio,
                          this->reduce_ratio, this->A_start_factor);

        newEquations = newEnum.getEquations();
        newEnumTime = this->timer.stopTimer();

        // statistics
        this->file.statisticSlightBKZ(slightBkzTime, newEnumTime);
        this->file.statisticsNewEquations(newEquations,this->primes);
        this->stats.newSlightBkzTime(slightBkzTime);
        this->stats.newNewEnumTime(newEnumTime, newEquations.size() > 0);

        newDuplicates = this->addEquations(newEquations);
        // display round results
        cout << " -> new equations:         " << newEquations.size() << "  |  " << this->uniqueEquations.size() << endl;
        cout << " -> duplicate equations:   " << newDuplicates << "  |  " << this->eqnDuplicates << endl;
    }
    this->stats.closeStatistics(round,this->uniqueEquations.size(),this->eqnDuplicates);

    this->file.writeFormattedEquationList(this->uniqueEquations, this->primes);

    this->file.writeSummary(this->stats, this->timer.getTotalTime());
    this->file.closeEquationFile();
    this->file.closeStatisticsFile();

    this->file.texToPdf();
}

/**
 * Adds new equations to the equation set and count duplicate equations.
 * @param newEquations
 * @return the number of duplicate equations from the input list.
 */
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
    long prime_num = n;

    this->primes.SetLength(prime_num);

    long prime = 1;
    for(long i = 1; i <= prime_num; i++)
    {
        prime = NextPrime(prime+1);
        this->primes(i) = prime;
    }
}

/**
 * Starts the program
 *
 * @param N This is going to be factorized
 * @param n The dimension
 * @param c The c of the prime number lattice
 * @param s_max The maximum pruning level
 * @param A_start_factor The upper bound of the distance is reduced by this factor bevor starting NewEnum
 * @param restart_ratio If the ratio of distances of a old close vector and a new close vector is lower
 * than this value, NewEnum will e restarted with a new maximal distance.
 * @param accuracy_factor Factor for the accurancy of the Basis
 * @param strong_bkz Block size for the string BKZ reduction
 * @param slight_bkz Block size for the slight BKZ reduction
 * @param seed_type -1 = random, -2 = timestamp (default), >= 0 cosntant
 */
Factoring::Factoring(ZZ N, long n, RR c, int s_max, double A_start_factor, double restart_ratio, double reduce_ratio, long accuracy_factor, long strong_bkz, long slight_bkz, unsigned long min_eqns, long long int seed_type)
{
    this->timer = Timer();
    // save the programm parameters
    this->N = N;
    this->c = c;
    this->A_start_factor = A_start_factor;
    this->restart_ratio = restart_ratio;
    this->reduce_ratio = reduce_ratio;
    this->slight_bkz = slight_bkz;
    this->s_max = s_max;
    // write the current settings

    this->setPrimes(n);

    // setup random number generator
    long long int seed;
    if(seed_type <= -2)
        seed = time(NULL);
    else if(seed_type == -1)
    {
        random_device rd;
        seed = rd();
    }
    else
        seed = seed_type;

    this->rgen = mt19937(seed);

    this->file.writeSettings(N, c, accuracy_factor, s_max, A_start_factor, reduce_ratio, strong_bkz, slight_bkz, n, this->primes(n), seed);

    if(min_eqns <= 0)
        this->min_eqns = n + 1;
    else
        this->min_eqns = min_eqns;

    this->file.prepareEquationTable();
    this->setBasis(accuracy_factor);
    this->reduceBasis(strong_bkz);
    this->setTargetCoordinates();
    this->search();
}

ZZ getN(long e)
{
    ZZ p1 = NextPrime(SqrRoot(power(conv<ZZ>(10), e)));
    ZZ p2 = NextPrime(p1 + 2);
    return p1 * p2;
}

