
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
#include <algorithm>

using namespace std;
using namespace NTL;

/**
 * struct NewEnumStage This is a data structure to save a stage of NewEnum
 */
struct NewEnumStage
{
    RR y_t;                     // equals y(t)
    RR c_tp1;                   // equals c(t+1)
    Vec<double> u;              // equals Vec<RR> u but need less memory. u contains only integers, so it's not important if they are stored in RR or double
    long t;                     // current coordinate

    NewEnumStage() {}

    NewEnumStage(const RR& y_t, const RR& c_tp1, const Vec<RR>& u, long t)
            : y_t(y_t), c_tp1(c_tp1), u(conv<Vec<double>>(u)), t(t)
    {}

    void set(const RR& y_t, const RR& c_tp1, const Vec<RR>& u, long t)
    {
        this->y_t = y_t;
        this->c_tp1 = c_tp1;
        conv(this->u, u);
        this->t = t;
    }

    /**
     * Converts the stored Vec<double> u to a Vec<RR> u, that is required by NewEnum
     * @return
     */
    void get_u(Vec<RR> &u) const;

};

class StageStorage
{
public:

    // statistics
    /**
     * Contains the number of stages stored for the level
     * at the moment this level gets performed.
     * Organized as [alpha_2_indicator][t_indicator][level-11]
     */
    std::vector<std::vector<std::vector<unsigned long long>>> maxDelayedAndPerformedStages;
    unsigned long long totalDelayedAndPerformedStages = 0;
    /**
     * Contains the number of stages stored for the level (performed and deleted)
     * Organized as [alpha_2_indicator][t_indicator][level-11]
     */
    std::vector<std::vector<std::vector<unsigned long long>>> delayedStages;

    std::vector<std::vector<double>> alpha_2_min = {{2.0,2.0,2.0},{2.0,2.0,2.0},{2.0,2.0,2.0}};     /**< organized as alpha_2_min[alpha_2_indicator][t_indicator] */

    StageStorage(unsigned long dim, unsigned long pruningLevel)
            : dim(dim), min_level(10), pruningLevel(pruningLevel),
              stageCounterByLevel(vector<unsigned long long>(pruningLevel - min_level,0))
    {
        this->storage = std::vector<std::vector<std::vector<std::list<NewEnumStage*>>>>(3,std::vector<std::vector<std::list<NewEnumStage*>>>(3, std::vector<std::list<NewEnumStage*>>(pruningLevel - min_level)));
        this->maxDelayedAndPerformedStages = std::vector<std::vector<std::vector<unsigned long long>>>(3, std::vector<std::vector<unsigned long long>>(3, std::vector<unsigned long long>(pruningLevel - min_level,0)));
        this->delayedStages = std::vector<std::vector<std::vector<unsigned long long>>>(3, std::vector<std::vector<unsigned long long>>(3, std::vector<unsigned long long>(pruningLevel - min_level,0)));
        // allocate 1000 stages for the beginning
        for(int i = 0; i < 1000; ++i)
            this->pool.push_back(new NewEnumStage);
    }

    bool getNext(NewEnumStage* &stage);

    void incrementCurrentLevel();

    void storeStage(const RR& y_t, const RR& c_t, const RR& c_tp1, const Vec<RR>& u, long t, long level);

    void updateMaxDistance(const RR &distance);

    void resetStorage(const RR &A)
    {
        // reset counters/statistics
        this->totalDelayedAndPerformedStages = 0;
        for(long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
            for(long t_indicator = 0; t_indicator < 3; t_indicator++)
                for(int level = 0; level < this->pruningLevel - this->min_level; level++)
                {
                    this->maxDelayedAndPerformedStages[alpha_2_indicator][t_indicator][level] = 0;
                    this->delayedStages[alpha_2_indicator][t_indicator][level] = 0;
                }

        this->alpha_2_min = {{2.0,2.0,2.0},{2.0,2.0,2.0},{2.0,2.0,2.0}};
        this->maxDistance = A;
        this->currentLevel = this->min_level;
    }

    /**
     * returns the stage to the pool of available stages. This prevents to much allocating
     */
    void returnStage(NewEnumStage* stage)
    {
        this->pool.push_back(stage);
    }

private:
    const unsigned long dim;

    const double alpha_1_threshold = 0.99;      /**< A_new/A_old = alpha_1 < alpha_1_threshold.
                                                     Used to prevent to much useless level recalculations. */

    const long pruningLevel;
    const long min_level = 10;                  /**< minimum level */

    unsigned long long stageCounterTotal = 0;
    vector<unsigned long long> stageCounterByLevel;

    long currentLevel = 10;

    RR maxDistance;

    list<NewEnumStage*> pool;
    std::vector<std::vector<std::vector<std::list<NewEnumStage*>>>> storage;        /**< A queue of stages for every level s and projection t
                                                                                         organized as storage[t_indicator][alpha_2_indicator][s-11]. */

    inline long alpha_2_indicator(const double &alpha_2)
    {
        return (alpha_2 > 0.65) ? 2 : (alpha_2 > 0.4 ? 1 : 0);
    }

    inline long t_indicator(const long &t)
    {
        return (t >= 40) ? 2 : (t >= 18 ? 1 : 0);
    }

    NewEnumStage* getStage()
    {
        if(this->pool.empty())
        {
            return new NewEnumStage;
        }

        NewEnumStage* stage = this->pool.front();
        this->pool.pop_front();
        return stage;
    }

    void recalculateLevels(const double &newDistance);

    inline unsigned long levelChange(const double &alpha_1, const double &alpha_2, long t_indicator);

    inline void returnStages(list<NewEnumStage*> &stages)
    {
        this->pool.splice(this->pool.end(),stages);
    }
};

class NewEnum
{

private:
    const FactoringSettings settings;

    // Output
    Timer& timer;
    long round;

    // Lattice Data
    Mat<RR> mu;                     /**< Contains the Gram-Schmidt coefficients */
    Vec<RR> R_ii_squared;           /**< Contains \f$\|\hat{b}_i\|^2 = r_{ii}^2\f$ */

    Vec<RR> tau;                    /**< The shifted target vector coordinates in the lattice basis */
    const Vec<RR> shift;            /**< The shift. This Vector contains only integral entries */
    bool decrease_max_distance;     /**< if false, the minimal distance will not decrease any more */
    double min_reduce_ratio;        /**< the algorithm reduces A_Curr, if the reduce ratio is lower (greater reduction) than this value,
                                         also after a equation was found and the algorithm is supposed not to reduce A_curr */

    const Mat<RR> U_RR;
    Mat<RR> U_scaled_RR;
    const long n;                   /**< Lattice dimension and the number of primes */

    // Factoring extension
    const ZZ N;
    const RR N_RR;
    const Vec<long>& primes;
    list<Equation> equations;


    // precomputed
    const Vec<double> log_V;                /**< Contains the values \f$V_t\f$ */
    Vec<double> log_V_minus_log_R_prod;     /**< Contains the values \f$log(V_t / (r_{1,1} * ... * r_{t,t}))\f$ */
    const Vec<double> log_t;                /**< contains the values log(t) for t = 1...n */

    const int min_level = 10;               /**< minimum level */
    const int max_level;                    /**< maximum level */
    int current_level = 10;                 /**< current level */
    RR A_curr;
    RR theoreticalMaxDistance;
    RR heuristicMaxDistance;

    // working space of checkForEquation()
    Vec<RR> close_vec, temp_vec;
    Vec<long> raw_equation, equation, close_vec_long;
    const ZZ threshold;
    ZZ u, vN, h_n, k_n, h_nm1, k_nm1, h_nm2,
       k_nm2, a_nm1, left_side, right_side, ride_side_factor;
    RR v, d, alpha_nm1;


    void prepare(unsigned long round, const Mat<ZZ> &newBasis_transposed, const Mat<ZZ> &newU_scaled,
                 const Vec<RR> &new_target_coordinates);

    /**
     * returns the closest integer to x -> \f$\lceil x \rfloor\f$
     * runs in \f$\mathcal{O}(1)\f$
     * @param x real number which to round
     * @return closest integer
     */
    inline void closest_RR (RR &out, const RR &x);

    /**
     * Returns the smallest integer to y with \f$|u-y| < |next(u,y) - y|\f$
     * runs in \f$\mathcal{O}(1)\f$
     * @param u integer with minimal distance
     * @param y the 'center'
     * @return the next integer
     */
    inline void next(RR &out, const RR &u, const RR &y);

    /**
     * performs all stages starting with a start stage
     * @param start The stage we are starting with
     * @param perform_delayed_stages If true, the method will also perform delayed
     * stages. Only the first call of this method should set this parameter to true.
     */
    inline void perform(NewEnumStage* start);

    /**
     * Precomputes some data used for the volume heuristic
     * runs in \f$\mathcal{O}(n)\f$
     * @param dim The dimension n of the lattice
     */
    static Vec<double> precomputeVolumes(long n);
    void precomputeLogV();

    static Vec<double> precomputeLogT(long n);

    /**
     * computes (if possible) an equation from a given close vector
     * @param [in] input Coordinates of a close vector
     * @param [in] c_1 Distance of this vector to the target vector
     * @return true if an equation was found, else false.
     */
    bool checkForEquation(const Vec<RR> &input, const RR &c_1);

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
     * @param [in,out] h_nm2
     * @param [in,out] k_nm2
     * @param [in,out] a_nm1
     * @param [in,out] alpha_nm1
     */
    void nextContinuedFraction(ZZ& h_n, ZZ& k_n, ZZ& h_nm1, ZZ& k_nm1, ZZ& h_nm2, ZZ& k_nm2, ZZ& a_nm1, RR& alpha_nm1);

    /**
     * Checks if k_n and left_side are p_n-smooth and stores the prime factors in equation.
     * @return If k_n and the left side are smooth true is returned, else false.
     */
    bool isSmooth(Vec<long>& equation, ZZ& k_n, ZZ& left_side, ZZ& right_side);

    inline long calculateLevel(long t, const RR &c_t)
    {
        return conv<long>(ceil(-((t-1)/2.0*log(conv<double>(this->A_curr - c_t)) + this->log_V_minus_log_R_prod(t - 1) - this->log_t(t)) / log(2)));
    }

public:
    /**
     * This is for performance reasons public. Never use this from outside of NewEnum
     */
    StageStorage L;
    unsigned long long stagesCheckedForEquations = 0;
    /**
     * Starts the NewEnum Algorithm to find close vectors and extract equations
     * @param primes        The vector of prime numbers used in the lattice
     * @param basis         The strong reduced, random scaled and slight
     *                      reduced lattice basis. Required for the Gram-Schmidt
     *                      coefficients and the length of the orthogonal basis vectors.
     * @param U             The transition matrix, that does the strong
     *                      BKZ reduction. Required to get the coordinates
     *                      of the close vector respecting  the prime number lattice.
     * @param target_shift  The shift that was done. required to shift
     *                      the close vector back where it should be.
     */
    NewEnum(const FactoringSettings &settings, Timer &timer, const Vec<long> &primes,
            const Mat<ZZ> &U, const Vec<RR> &target_shift);

    /**
     * Returns a list of equations that were found in this round.
     */
    list<Equation> getEquations();

    void getDistances(RR &theoretical,RR &heuristic,RR &reduced);
    /**
     * starts the performing by setting the first stage to the values that are
     * described under point 1 in the NewEnum algorithm.
     */
    void run(unsigned long round, const Mat<ZZ> &newBasis_transposed, const Mat<ZZ> &new_U_scaled,
             const Vec<RR> &new_target_coordinates);

};

#endif	/* NEWENUM_H */