
#ifndef FACTORINGINTEGERS_H
#define	FACTORINGINTEGERS_H

#include "NewEnum.h"
#include "Equation.h"
#include "Timer.h"
#include "FileOutput.h"
#include "Statistics.h"

#include <NTL/LLL.h>

#include <random>

using namespace std;
using namespace NTL;

struct FactoringSettings
{
    ZZ N;
    unsigned long n;
    RR c;
    int s_max = 14;
    double A_start_factor = 0.2;
    double restart_ratio = 0;
    double reduce_ratio = 0.8;
    long accuracy_factor = 10000;
    long strong_bkz = 32;
    long slight_bkz = 20;
    unsigned long min_eqns;
    long long int seed_type = -2;

    FactoringSettings(const ZZ &N, unsigned long n, const RR &c) : N(N), n(n), c(c), min_eqns(n+1)
    { }

    FactoringSettings(const ZZ &N, unsigned long n, const RR &c, unsigned long min_eqns) : N(N), n(n), c(c), min_eqns(min_eqns)
    { }

    FactoringSettings(const ZZ &N, unsigned long n, const RR &c, int s_max, double A_start_factor,
                      double restart_ratio, double reduce_ratio, long accuracy_factor,
                      long strong_bkz, long slight_bkz, unsigned long min_eqns) : N(N), n(
            n), c(c), s_max(s_max), A_start_factor(A_start_factor), restart_ratio(
            restart_ratio), reduce_ratio(reduce_ratio), accuracy_factor(
            accuracy_factor), strong_bkz(strong_bkz), slight_bkz(slight_bkz), min_eqns(
            min_eqns)
    { }

    FactoringSettings(const ZZ &N, unsigned long n, const RR &c, int s_max, double A_start_factor,
                      double restart_ratio, double reduce_ratio, long accuracy_factor,
                      long strong_bkz, long slight_bkz, unsigned long min_eqns,
                      long long int seed_type) : N(N), n(n), c(c), s_max(
            s_max), A_start_factor(A_start_factor), restart_ratio(restart_ratio), reduce_ratio(
            reduce_ratio), accuracy_factor(accuracy_factor), strong_bkz(
            strong_bkz), slight_bkz(slight_bkz), min_eqns(min_eqns), seed_type(seed_type)
    { }
};

class Factoring
{
private:

    ZZ N;                               /**< the number which is going to be factorized */
    RR c;                               /**< */
    double A_start_factor;              /**< */
    double restart_ratio;               /**< */
    double reduce_ratio;                /**< */

    long slight_bkz;                    /**< */
    Vec<long> primes;                   /**< the primes used for factoring */

    unsigned long min_eqns;             /**< the program stops after finding at least that much equations*/
    mt19937 rgen;

    // Parameters used in NewEnum
    int s_max;                          /**< The maximum pruning level */
    // Mat<ZZ> B_scaled;                /**< The scaled lattice basis */

    // Data often used
    Mat<ZZ> B;                          /**< The used strong reduced prim lattice basis */
    Mat<ZZ> U;                          /**< The transition matrix with B_old*U = BKZ */
    Mat<ZZ> U_inv;                      /**< The inverse of U */
    Vec<RR> target_coordinates;         /**< The coordinates of th target vector (reduced and shifted) */
    Vec<RR> shift;                      /**< The shift */

    // Buffer
    Mat<ZZ> B_scaled;                   /**< The scaled BKZ-basis */
    Mat<ZZ> U_scaled;                   /**< The transition matrix to the scaled BKZ basis */
    Mat<ZZ> U_scaled_inv;               /**< The inverse of U_scale */
    Vec<RR> target_scaled_coordinates;  /**< The scaled shifted target vector coordinates */

    // Output
    Timer timer;
    FileOutput file;
    std::set<Equation> uniqueEquations;       /**< contains the equations that will be found */


    Statistics stats;                   /**< all the statistics */
    long eqnDuplicates = 0;

    /**
     * scales the lattice basis, BKZ-reduces it slightly, and also converting the
     * target coordinates
     */
    void randomScale();

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
     * @param s_max The maximum pruning level
     * @param A_start_factor The upper bound of the distance is reduced by this factor before starting NewEnum
     * @param restart_ratio If the ratio of distances of a old close vector and a new close vector is lower
     * than this value, NewEnum will e restarted with a new maximal distance.
     * @param accuracy_factor Factor for the accuracy of the Basis
     * @param strong_bkz Block size for the string BKZ reduction
     * @param slight_bkz Block size for the slight BKZ reduction
     */
    Factoring(ZZ N, long n, RR c, int s_max, double A_start_factor, double restart_ratio, double reduce_ratio, long accuracy_factor, long strong_bkz, long slight_bkz, unsigned long min_eqns, long long int seed_type = -2);

};

/**
 * Generates a N containing exact two prime factors of size approximately \f$10^{e/2}\f$
 * @param e
 * @return N
 */
ZZ getN(long e);

#endif	/* FACTORINGINTEGERS_H */

