
#ifndef FACTORINGINTEGERS_H
#define	FACTORINGINTEGERS_H

#include "NewEnum.h"
#include "Equation.h"
#include "Timer.h"
#include "FileOutput.h"
#include "Statistics.h"
#include "FactoringSettings.h"

#include <NTL/LLL.h>

#include <random>

using namespace std;
using namespace NTL;

class Factoring
{
private:
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

protected:

    const FactoringSettings settings;
    const ZZ N;                         /**< the number which is going to be factorized */
    const RR c;                         /**< */

    Vec<long> primes;                   /**< the primes used for factoring */

    mt19937 rgen;

    // Data often used
    Mat<ZZ> B,                          /**< The used strong reduced prim lattice basis */
            U,                          /**< The transition matrix with B_old*U = BKZ */
            U_inv;                      /**< The inverse of U */
    Vec<RR> target_coordinates,         /**< The coordinates of th target vector (reduced and shifted) */
            shift;                      /**< The shift */

    // Buffer
    Mat<ZZ> B_scaled_transposed,        /**< The scaled BKZ-basis */
            U_scaled,                   /**< The transition matrix to the scaled BKZ basis */
            U_scaled_inv;               /**< The inverse of U_scale */
    Vec<RR> target_scaled_coordinates;  /**< The scaled shifted target vector coordinates */
    vector<bool> scaled_primes;

    // Output
    Timer timer;
    FileOutput file;
    std::set<Equation> uniqueEquations;         /**< contains the equations that will be found */


    Statistics stats;                           /**< all the statistics */
    long eqnDuplicates = 0;

    /**
     * scales the lattice basis, BKZ-reduces it slightly, and also converting the
     * target coordinates
     */
    void setScaledAndReducedBasis()
    {
        // scale
        this->B_scaled_transposed = this->B;    // not yet transposed

        long threshold1,threshold2;
        int scalingType;
        if(this->settings.scalingType == FS_SCALING_MIXED)
        {
            uniform_int_distribution<int> scalingType_dist(0,49);
            int choice = scalingType_dist(this->rgen);
            if(choice <= 2)         // 6%
                scalingType = FS_SCALING_ONE_FOURTH;
            else if(choice <= 6)    // 8%
                scalingType = FS_SCALING_ONE_FOURTH_ONE_HALF;
            else if(choice <= 10)   // 8%
                scalingType = FS_SCALING_ONE_HALF_ONE_FOURTH;
            else                    // 78%
                scalingType = FS_SCALING_ONE_HALF;
        }
        else
        {
            scalingType = this->settings.scalingType;
        }

        switch(scalingType)
        {
            case FS_SCALING_ONE_FOURTH:
                threshold1 = 0;
                threshold2 = 0;
                break;
            case FS_SCALING_THREE_FORTH:
                threshold1 = 2;
                threshold2 = 2;
                break;
            case FS_SCALING_ONE_FOURTH_ONE_HALF:
                threshold1 = 0;
                threshold2 = 1;
                break;
            case FS_SCALING_ONE_HALF_ONE_FOURTH:
                threshold1 = 1;
                threshold2 = 0;
                break;
            default:    // FS_SCALING_ONE_HALF
                threshold1 = 1;
                threshold2 = 1;
                break;
        }

        uniform_int_distribution<int> dist(0,3);   // [0,3] random numbers uniform distribution

        for(long i = 1; i <= this->settings.n/2; i++)
        {
            if(dist(this->rgen) <= threshold1)       // P(scale row) = 1/2
            {
                this->B_scaled_transposed(i) *= 2;     // multiply the whole row by 2
                this->scaled_primes[i-1] = true;
            }
            else
            {
                this->scaled_primes[i-1] = false;
            }
        }

        for(long i = this->settings.n/2+1; i <= this->settings.n; i++)
        {
            if(dist(this->rgen) <= threshold2)       // P(scale row) = 1/4
            {
                this->B_scaled_transposed(i) *= 2;     // multiply the whole row by 2
                this->scaled_primes[i-1] = true;
            }
            else
            {
                this->scaled_primes[i-1] = false;
            }
        }

        // reduce
        transpose(this->B_scaled_transposed, this->B_scaled_transposed);    // now it is transposed
        BKZ_FP(this->B_scaled_transposed, this->U_scaled, 0.99, this->settings.slight_bkz);

        transpose(this->U_scaled,this->U_scaled);
        inv(this->U_scaled_inv, this->U_scaled);

        // apply the reduction to the coordinate vector
        mul(this->target_scaled_coordinates, conv<Mat<RR>>(this->U_scaled_inv), this->target_coordinates);
    }


    /**
     * Generates the prime lattice basis and make a strong BKZ reduction with
     * block size strong_bkz. The entries of the basis have to be (big) integers
     * so the basis is multiplied by accuracy_factor before rounding.
     * @param accuracy_factor
     */
    void setBasis(const long accuracy_factor)
    {
        cout << "Setting Prime Lattice Basis: ";
        long n = this->primes.length();

        RR NpC;
        conv(NpC, this->N);
        pow(NpC, NpC, this->c);
        this->B.SetDims(n, n + 1);    // transposed
        // Setting the basis
        for(long i = 1; i <= n; i++)
        {
            conv(this->B(i, i), accuracy_factor * sqrt(log(this->primes(i))));
            conv(this->B(i, n + 1), NpC * accuracy_factor * log(this->primes(i)));
        }

        cout << "finished" << endl;
    }

    /**
     * Computes the coordinates of the target vector projected orthogonally into
     * the lattice plane. This method runs fast since the coordinates of the
     * projection are all the same an can be computed by a formula. NOTICE:
     * This method works only for the standard target vector.
     * @return Coordinates of the projected target vector.
     */
    Vec<RR> orthogonalProjection_pl() const
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

    /**
     * Computes the coordinates of the target vector respecting the strong reduced
     * basis and shifted to the ground mesh. Requires the the transition matrix
     * U_inv to be set.
     */
    void setTargetCoordinates()
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
    void reduceBasis(const long strong_bkz)
    {
        this->timer.startTimer();

        // this->B is still transposed

        cout << "Strong BKZ reduction: ";
        BKZ_FP(this->B, this->U, 0.99, strong_bkz);                 // strong reducing
//    BKZ_QP(this->B, this->U, 0.99, strong_bkz);                 // strong reducing

        transpose(this->B,this->B);
        transpose(this->U,this->U);

        inv(this->U_inv, this->U);

        cout << "finished" << endl;

        this->file.statisticsStrongBkzTime(this->timer.stopTimer());
    }

    /**
     * Runs the search for relations. This method runs a loop until 3*n relations
     * were found
     */
    void search()
    {
        list<Equation> newEquations;
        double newEnumTime;
        double slightBkzTime;
        long newDuplicates;
        unsigned long round = 0;

        NewEnum newEnum = NewEnum(this->settings, this->timer, this->primes, this->U, this->shift);

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
            this->stats.updateRoundStats(newEnum.L.totalDelayedAndPerformedStages > 0, !newEquations.empty());
            this->stats.updateDistanceStats(theoretical,heuristic,reduced);
            this->stats.updateStagesCheckedForEquations(newEnum.stagesCheckedForEquations, !newEquations.empty());
            this->stats.newSlightBkzTime(slightBkzTime);
            this->stats.newNewEnumTime(newEnumTime, !newEquations.empty());
            this->file.statisticsDelayedStagesOnLevel(this->settings.max_level,
                                                      newEnum.L.maxDelayedAndPerformedStages,
                                                      newEnum.L.delayedStages,
                                                      newEnum.L.totalDelayedAndPerformedStages);
            this->file.statisticSlightBKZ(slightBkzTime, newEnumTime);
            this->file.statisticsWriteStagesChecked(newEnum.stagesCheckedForEquations);
            this->file.statisticsDistances(theoretical,heuristic,reduced);
            this->file.statisticsWriteScaledPrimes(this->scaled_primes,this->primes);
            this->file.statisticsNewEquations(newEquations,this->primes);

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


    /**
     * Adds new equations to the equation set and count duplicate equations.
     * @param newEquations
     * @return the number of duplicate equations from the input list.
     */
    long addEquations(const list<Equation>& newEquations)
    {
        long duplicates = 0;
        for(auto& newEquation : newEquations)
        {
            std::pair<std::set<Equation>::iterator, bool> res = this->uniqueEquations.insert(newEquation);
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

    void setPrimes(long n)
    {
        this->primes.SetLength(n);

        for(long i = 0; i < min(n,300); ++i)
            this->primes[i] = _primes[i];

        for(long i = 301; i <= n; i++)
            this->primes(i) = NextPrime(this->primes(i-1)+2);
    }

public:

    /**
     * Starts the program
     *
     * @param N This is going to be factorized
     * @param n The dimension
     * @param c The c of the prime number lattice
     * @param max_level The maximum pruning level
     * @param A_start_factor The upper bound of the distance is reduced by this factor before starting NewEnum
     * @param restart_ratio If the ratio of distances of a old close vector and a new close vector is lower
     * than this value, NewEnum will e restarted with a new maximal distance.
     * @param accuracy_factor Factor for the accuracy of the Basis
     * @param strong_bkz Block size for the string BKZ reduction
     * @param slight_bkz Block size for the slight BKZ reduction
     */
    explicit Factoring(const FactoringSettings& settings)
        : settings(settings),
          N(settings.N),
          c(settings.c),
          scaled_primes(vector<bool>(settings.n, false))
    {
        this->timer = Timer();

        this->setPrimes(settings.n);

        // setup random number generator
        unsigned int seed;
        if(settings.seed_type <= -2)
            seed = (unsigned int) time(nullptr);
        else if(settings.seed_type == -1)
        {
            random_device rd;
            seed = rd();
        }
        else
            seed = (unsigned int) settings.seed_type;

        this->rgen = mt19937(seed);

        // write the current settings
        this->file.writeSettings(this->settings, this->primes(settings.n), seed);

        this->setBasis(settings.accuracy_factor);
        this->reduceBasis(settings.strong_bkz);
        this->setTargetCoordinates();
        this->search();
    }
};

/**
 * Generates a N containing exact two prime factors of size approximately \f$10^{e/2}\f$
 * @param e
 * @return N
 */
inline ZZ getN(long e)
{
    ZZ p1 = NextPrime(SqrRoot(power(conv<ZZ>(10), e)));
    ZZ p2 = NextPrime(p1 + 2);
    return p1 * p2;
}

#endif	/* FACTORINGINTEGERS_H */