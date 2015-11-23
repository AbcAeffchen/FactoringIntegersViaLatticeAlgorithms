
#ifndef FACTORINGINTEGERS_H
#define	FACTORINGINTEGERS_H

#define __NUM_THREADS__ 2

#include "NewEnum.h"
#include "Equation.h"
#include "Timer.h"
#include "FileOutput.h"
#include "Statistics.h"
#include "FactoringSettings.h"

#include <NTL/LLL.h>

#include <random>
#include <omp.h>

using namespace std;
using namespace NTL;

class Factoring
{
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

    // Output
    Timer timer;
    FileOutput file;
    std::set<Equation> uniqueEquations;         /**< contains the equations that will be found */


    Statistics stats;                           /**< all the statistics */
    long eqnDuplicates = 0;

    /**
     * Scales every row with probability 1/2
     */
    void randomScale1(Mat<ZZ> &basis);
    /**
     * Scales every row with probability 1/4
     */
    void randomScale2(Mat<ZZ> &basis);

    /**
     * scales the lattice basis, BKZ-reduces it slightly, and also converting the
     * target coordinates
     */
    void setScaledAndReducedBasis(Mat<ZZ> &B_scaled_transposed,Mat<ZZ> &U_scaled,Vec<RR> &target_scaled_coordinates);

    /**
     * Generates the prime lattice basis and make a strong BKZ reduction with
     * block size strong_bkz. The entries of the basis have to be (big) integers
     * so the basis is multiplied by accuracy_factor before rounding.
     * @param accuracy_factor
     */
    void setBasis(long accuracy_factor);

    /**
     * Computes the coordinates of the target vector projected orthogonally into
     * the lattice plane. This method runs fast since the coordinates of the
     * projection are all the same an can be computed by a formula. NOTICE:
     * This method works only for the standard target vector.
     * @return Coordinates of the projected target vector.
     */
    Vec<RR> orthogonalProjection_pl();

    /**
     * Computes the coordinates of the target vector respecting the strong reduced
     * basis and shifted to the ground mesh. Requires the the transition matrix
     * U_inv to be set.
     */
    void setTargetCoordinates();

    /**
     * Performs a strong BKZ reduction and sets this->U
     * @param strong_bkz Block size of the BKZ reduction
     */
    void reduceBasis(long strong_bkz);

    /**
     * Runs the search for relations. This method runs a loop until 3*n relations
     * were found
     */
    void search();

    /**
     * Adds new equations to the equation set and count duplicate equations.
     * @param newEquations
     * @return the number of duplicate equations from the input list.
     */
    long addEquations(list<Equation>& newEquations);

    void setPrimes(long n);

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
    Factoring(const FactoringSettings &settings);

};

/**
 * Generates a N containing exact two prime factors of size approximately \f$10^{e/2}\f$
 * @param e
 * @return N
 */
ZZ getN(long e);

#endif	/* FACTORINGINTEGERS_H */