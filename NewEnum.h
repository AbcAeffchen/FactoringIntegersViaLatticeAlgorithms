/*
 * File:   NewEnum.h
 * Author: AbcAeffchen
 *
 * Created on 21. September 2014, 17:13
 */

#ifndef NEWENUM_H
#define	NEWENUM_H

#include "Equation.h"
#include "Timer.h"
#include "Statistics.h"
#include "FileOutput.h"

#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <NTL/LLL.h>
#include <math.h>
#include <list>
#include <queue>
#include <vector>
#include <iostream>


using namespace std;
using namespace NTL;

/**
 * struct NewEnumStage This is a data structure to save a stage of NewEnum
 */
struct NewEnumStage
{
    RR y_t;                     // equals y(t)
    RR c_t;                     // equals c(t)
    RR c_tp1;                   // equals c(t+1)
    Vec<double> u;              // equals Vec<RR> u but need less momory. u contains only integers, so its not importend if they are stored in RR or double
    unsigned long t;            // current coordinate
    Vec<RR> target;             // the coordinates of the target vector
    RR v_adjust;                // the v value used to adjust the target vector

    NewEnumStage(const RR& y_t, const RR& c_t, const RR& c_tp1, const Vec<RR>& u, unsigned long t, Vec<RR> target, RR v_adjust);

    /**
     * Convertrs the stored Vec<double> u to a Vec<RR> u, that is required by NewEnum
     * @return
     */
    Vec<RR> get_u();

};

class NewEnum
{

private:
    // Output
    Timer& timer;
    FileOutput& file;
    Statistics& stats;
    long round;
    // Lattice Data

    Mat<RR> mu;                     /**< Contains the Gram-Schmidt coefficients */
    Vec<RR> R_ii_squared;           /**< Contains \f$\|\hat{b}_i\|^2 = r_{ii}^2\f$ */

    Vec<RR> tau;                    /**< The shifted target vector coordinates in the lattice basis */
    Vec<RR> shift;                  /**< The shift. This Vector contains only integral entrys */
    bool decrease_max_distance;     /**< if false, the minimal distance will not decrease any more */
    double min_restart_ratio;       /**< the algorithm starts from the beginning, if the new close vector is at least this factor closer.
                                         min_restart_factor = 0 -> never restart, min_restart_factor = 1 -> always restart */
    double min_reduce_ratio;        /**< the algorithm reduces A_Curr, if the reduce ratio is lower (greater reduction) than this value,
                                         also after a equation was found and the algorithm is supposed not to reduce A_curr */

    const Mat<ZZ>& B;               /**< The reduced, scacled and again reduced Basis */
    const Mat<ZZ>& U;               /**< Transitionmatrix of the strong BKZ reduction */
    const Mat<ZZ>& U_scaled;        /**< Transitionmatrix of the slight BKZ reduction */
    long n;                         /**< Lattice dimension and the number of prims */

    // Factoring extension
    ZZ N;
    const Vec<long>& primes;
    list<Equation> equations;

    // statistics
    long stages = 0;                    /**< Counts the stages in total */

    vector<queue<NewEnumStage> > L;     /**< array of lists of delayd stages */
    queue<NewEnumStage> restart_list;   /**< queue of stages that require to be performed from start */
    // precomputed
    Vec<RR> V;                          /**< Contains the values \f$V_t / (r_11 * ... * r_tt)\f$ */
    Vec<RR> level_probabilitys;         /**< the probability a stage must have to belong to the levels */

    int s_max;                          /**< maximum level */
    int current_s = 10;                 /**< current level */
    RR A_curr;
    RR theoreticalMaxDistance;
    vector<long> delayedStagesCounter;

    // working space for getEquation
    Vec<RR> close_vec, temp_vec;
    Vec<long> raw_equation, equation;
    ZZ u,threshold,vN,h_n,k_n,h_nm1,k_nm1,h_nm2,k_nm2,a_nm1,left_side,right_side,ride_side_factor;
    RR v,d,alpha_nm1;


    void precomputeLevelProbabilitys();

    /**
     * returns the closest interger to x -> \f$\lceil x \rfloor\f$
     * runs in \f$\mathcal{O}(1)\f$
     * @param x realnumber which to round
     * @return closest integer
     */
    RR closest_RR (RR x);

    /**
     * Returns the smalest integer to y with \f$|u-y| < |next(u,y) - y|\f$
     * runs in \f$\mathcal{O}(1)\f$
     * @param u integer with minimal distance
     * @param y the 'center'
     * @return the next integer
     */
    RR next(RR u, RR y);

    /**
     * clears the list of delayed stages
     */
    void clearL();

    /**
     *
     * @param u
     * @param v
     * @param z
     */
    void computeUV(const Vec<RR>& input, ZZ& u, ZZ& v);

    /**
     * starts the performing by setting the first stage to the values that are
     * described under point 1 in the NewEnum algorithm.
     */
    void perform();

    /**
     * performs all stages starting with a start stage
     * @param start The stage we are starting with
     * @param perform_delayed_stages If true, the method will also perform delayed
     * stages. Only the first call of this method should set this parameter to true.
     */
    void perform(NewEnumStage& start, bool perform_delayed_stages);

    /**
     * Precomputes some data used for the volumeheuristic
     * runs in \f$\mathcal{O}(n)\f$
     * @param dim The dimension n of the lattice
     */
    void precomputeVolumes(long dim);

    /**
     * computes (if possible) an equation from a given close vector
     * @param input Coordinates of a close vector
     * @return true if an equation was found, else false.
     */
    bool getEquation(const Vec<RR>& input, RR& c_1 );

    /**
     * Computes the next continued fraction $h_n / k_n$ with $h_n = a_n * h_nm1 + h_nm2$
     * and $k_n = a_n * k_nm1 + k_nm2$. $a_n$ is computed as $a_n = floor( (alpha_nm1 - a_nm1)^-1 )$
     * It also sets $h_nm2 = h_nm1$ and $h_nm1 = k_n$ and $k_nm2 = k_nm1$ and $k_nm1 = k_n$ and
     * $a_nm1 = a_n$ and $alpha_nm1 = (alpha_nm1 - a_nm1)^-1$, so the function can be called
     * multiple times. To work correct, the function has to be called with correct values of the
     * n-1th and n-2th path of the continued fraction. If the computation starts at all, it has
     * to be called with h_nm1 = 1, k_nm1 = 0, h_nm2 = 0, k_nm2 = 1, a_nm1 = floor(x) and alpha_nm1 = x
     * @param [out] h_n
     * @param [out] k_n
     * @param [in,out] h_nm1
     * @param [in,out] k_nm1
     * @param [in] h_nm2
     * @param [in] k_nm2
     * @param [in] a_nm1
     * @param [in,out] alpha_nm1
     */
    void nextContinuedFraction(ZZ& h_n, ZZ& k_n, ZZ& h_nm1, ZZ& k_nm1, ZZ& h_nm2, ZZ& k_nm2, ZZ& a_nm1, RR& alpha_nm1);

    /**
     * Checks if k_n and left_side are p_n-smooth and stores the prime factors in equation.
     * @return If k_n and the left side are smooth true is returned, else false.
     */
    bool isSmooth(Vec<long>& equation, ZZ& k_n, ZZ& left_side, ZZ& right_side);

    Vec<RR> getCurrentTargetCoordinates(RR v_old, RR v_new);

public:

    /**
     * Starts the NewEnum Algorithm to find close vectors and extract equations
     * @param file                  Reference to the file object, that handls th
     * data output.
     * @param N                     The N that is going to be factorized
     * @param primes                The vector of prime numbers used in the lattice
     * @param basis                 The strong reduced, random scaled and slight
     * reduced lattice basis. Required for the gram-schmidt-coefficents and the
     * length of the orthogonal basis vectors.
     * @param U                     The transition matrix, that does the strong
     * BKZ reduction. Required to get the coordinates of the close vector respecting
     * the prime number lattice.
     * @param U_scaled              The transition matrix, that does the slight
     * BKZ reduction. Also required to get the coordinats of the close vector.
     * @param target                The coordinates of the (shifted) target vector
     * @param target_shift          The shift that was done. required to shift
     * the close vector back where it should be.
     * @param s_max                 The maximal pruning level.
     * @param min_restart_ratio     The minimal ratio that is required to restart
     * NewEnum.
     * @param start_factor_A        The factor to set the minimal distance at the
     * beginning.
     */
    NewEnum(Timer& timer, FileOutput& file, Statistics& stats, long round, ZZ N, const Vec<long>& primes, const Mat<ZZ>& basis,
        const Mat<ZZ>& U, const Mat<ZZ>& U_scaled, Vec<RR> target_coordinates, Vec<RR> target_shift,
        int s_max, double min_restart_ratio, double min_reduce_ratio, double start_factor_A );

    /**
     * Returns a list of equations that were found in this round.
     */
    list<Equation> getEquations();

};

#endif	/* NEWENUM_H */

