
#ifndef NEWENUM_H
#define    NEWENUM_H

#include "Equation.h"
#include "Timer.h"
#include "Statistics.h"
#include "FileOutput.h"

#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <NTL/LLL.h>
#include <cmath>
#include <list>
#include <queue>
#include <utility>
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

    NewEnumStage() : y_t(conv<RR>(0)), c_tp1(conv<RR>(0)), t(0)
    {}

    NewEnumStage(RR y_t, RR c_tp1, const Vec<RR>& u, long t)
        : y_t(std::move(y_t)), c_tp1(std::move(c_tp1)), u(conv<Vec<double>>(u)), t(t)
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
    void get_u(Vec<RR>& u) const
    {
        conv(u, this->u);
    }


};

class StageStorage
{
public:
    // todo improve this storage by flattening. This should reduce the access times.
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

    explicit StageStorage(unsigned long pruningLevel)
        : pruningLevel(pruningLevel), min_level(10),
          stageCounterByLevel(vector<unsigned long long>(pruningLevel - min_level, 0))
    {
        storage = std::vector<std::vector<std::vector<std::list<NewEnumStage*>>>>(3,
                                                                                  std::vector<std::vector<std::list<NewEnumStage*>>>(
                                                                                      3,
                                                                                      std::vector<std::list<NewEnumStage*>>(
                                                                                          pruningLevel - min_level)));
        maxDelayedAndPerformedStages = std::vector<std::vector<std::vector<unsigned long long>>>(3,
                                                                                                 std::vector<std::vector<unsigned long long>>(
                                                                                                     3,
                                                                                                     std::vector<unsigned long long>(
                                                                                                         pruningLevel -
                                                                                                         min_level,
                                                                                                         0)));
        delayedStages = std::vector<std::vector<std::vector<unsigned long long>>>(3,
                                                                                  std::vector<std::vector<unsigned long long>>(
                                                                                      3,
                                                                                      std::vector<unsigned long long>(
                                                                                          pruningLevel - min_level,
                                                                                          0)));
        alpha_2_min = std::vector<std::vector<std::vector<double>>>(3, std::vector<std::vector<double>>(3,
                                                                                                        std::vector<double>(
                                                                                                            pruningLevel -
                                                                                                            min_level,
                                                                                                            2.0)));
        // allocate 1000 stages for the beginning
        for(int i = 0; i < 200000; ++i)
            pool.push_back(new NewEnumStage);
    }

    ~StageStorage()
    {
        resetStorage();   // return all stages to the pool
        NewEnumStage* temp;
        // delete all stages
        while(!pool.empty())
        {
            temp = pool.front();
            pool.pop_front();
            delete temp;
        }
    }

    bool getNext(NewEnumStage*& stage)
    {
        if(stageCounterByLevel[currentLevel - min_level - 1] <= 0 || stageCounterTotal >= stageCounterThreshold)
            return false;

        for(unsigned long t_indicator = 0; t_indicator < 3; ++t_indicator)
        {
            for(unsigned long alpha_2_indicator = 0; alpha_2_indicator < 3; ++alpha_2_indicator)
            {
                if(!storage[t_indicator][alpha_2_indicator][currentLevel - min_level - 1].empty())
                {
                    stage = storage[t_indicator][alpha_2_indicator][currentLevel - min_level - 1].front();
                    storage[t_indicator][alpha_2_indicator][currentLevel - min_level - 1].pop_front();
                    stageCounterByLevel[currentLevel - min_level - 1]--;
                    stageCounterTotal--;
                    return true;
                }
            }
        }

        // It never should come this far, but who knows...
        return false;
    }

    void incrementCurrentLevel()
    {
        currentLevel++;

        for(unsigned long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
            for(unsigned long t_indicator = 0; t_indicator < 3; t_indicator++)
            {
                maxDelayedAndPerformedStages[alpha_2_indicator][t_indicator][currentLevel - min_level -
                                                                             1] = storage[t_indicator][alpha_2_indicator][
                    currentLevel - min_level - 1].size();
                totalDelayedAndPerformedStages += storage[t_indicator][alpha_2_indicator][currentLevel - min_level -
                                                                                          1].size();
            }
    }

    bool storeStage(const RR& y_t, const RR& c_t, const RR& c_tp1, const Vec<RR>& u, long t, long level)
    {
        NewEnumStage* stage = getStage();
        stage->set(y_t, c_tp1, u, t);

        double alpha_2;
        conv(alpha_2, c_t / maxDistance);
        long alpha_2_indicator = this->alpha_2_indicator(alpha_2);
        long t_indicator = this->t_indicator(stage->t);
        storage[t_indicator][alpha_2_indicator][level - min_level - 1].push_back(stage);

        if(alpha_2_min[alpha_2_indicator][t_indicator][level - min_level - 1] > alpha_2)
            alpha_2_min[alpha_2_indicator][t_indicator][level - min_level - 1] = alpha_2;

        stageCounterTotal++;
        stageCounterByLevel[level - min_level - 1]++;

        // statistics
        delayedStages[alpha_2_indicator][t_indicator][level - min_level - 1]++;

        return stageCounterTotal >= stageCounterThreshold;      // about 6 GB of RAM depending on the dimension
    }

    void updateMaxDistance(const RR& distance)
    {
        double alpha_1;
        conv(alpha_1, distance / maxDistance);

        if(alpha_1 < alpha_1_threshold && stageCounterTotal > 100)
            recalculateLevels(alpha_1);

        for(long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
            for(long t_indicator = 0; t_indicator < 3; t_indicator++)
                for(long s = std::max(0L, currentLevel - min_level - 1); s < pruningLevel - min_level; s++)
                    alpha_2_min[alpha_2_indicator][t_indicator][s] /= alpha_1;

        maxDistance = distance;
    }

    void resetStorage(const RR& A = conv<RR>(0))
    {
        // reset counters/statistics
        totalDelayedAndPerformedStages = 0;
        stageCounterTotal = 0;
        for(int level = 0; level < pruningLevel - min_level; level++)
            stageCounterByLevel[level] = 0;

        for(long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
            for(long t_indicator = 0; t_indicator < 3; t_indicator++)
                for(int level = 0; level < pruningLevel - min_level; level++)
                {
                    maxDelayedAndPerformedStages[alpha_2_indicator][t_indicator][level] = 0;
                    delayedStages[alpha_2_indicator][t_indicator][level] = 0;
                    alpha_2_min[alpha_2_indicator][t_indicator][level] = 2.0;
                    returnStages(storage[t_indicator][alpha_2_indicator][level]);
                }

        maxDistance = A;
        currentLevel = min_level;
    }

    /**
     * returns the stage to the pool of available stages. This prevents to much allocating
     */
    void returnStage(NewEnumStage* stage)
    {
        pool.push_back(stage);
    }

    unsigned long long getTotalDelayedStages()
    {
        return stageCounterTotal;
    }

private:
    const double alpha_1_threshold = 0.99;      /**< A_new/A_old = alpha_1 < alpha_1_threshold.
                                                     Used to prevent to much useless level recalculations. */
    const unsigned long long stageCounterThreshold = 2500000;   // About 6 GB of RAM depending on the dimension
    const long pruningLevel;
    const long min_level = 10;                  /**< minimum level */

    unsigned long long stageCounterTotal = 0;
    vector<unsigned long long> stageCounterByLevel;

    long currentLevel = 10;

    RR maxDistance;

    std::vector<std::vector<std::vector<double>>> alpha_2_min;     /**< organized as alpha_2_min[alpha_2_indicator][t_indicator][level] */

    list<NewEnumStage*> pool;
    std::vector<std::vector<std::vector<std::list<NewEnumStage*>>>> storage;        /**< A queue of stages for every level s and projection t
                                                                                         organized as storage[t_indicator][alpha_2_indicator][s-11]. */

    inline long alpha_2_indicator(const double alpha_2) const
    {
//        const bool a = alpha_2 > 0.65;
        // return static_cast<long>(a) * 2 + static_cast<long>(!a) * static_cast<long>(alpha_2 > 0.4);
        return (alpha_2 > 0.65) ? 2 : (alpha_2 > 0.4 ? 1 : 0);
    }

    inline long t_indicator(const long t) const
    {
//        const bool a = t >= 40;
        // return static_cast<long>(a) * 2 + static_cast<long>(!a) * static_cast<long>(t >= 18);

        return (t >= 40) ? 2 : (t >= 18 ? 1 : 0);
    }

    NewEnumStage* getStage()
    {
        if(pool.empty())
        {
            for(long i = 0; i < 10000; i++)
                pool.push_back(new NewEnumStage);
        }

        NewEnumStage* stage = pool.front();
        pool.pop_front();
        return stage;
    }

    void recalculateLevels(const double alpha_1)
    {
        double new_alpha_2;
        long new_alpha_2_indicator;
        long level_change;
        long s_max = pruningLevel - min_level - 1;
        long s_min = currentLevel - min_level - 1;

        /**
         * skip t_indicator = 0, since the stages here have t between 4 and 20.
         * So in almost all cases the will be no change.
         */
        for(int t_indicator = 2; t_indicator >= 1; t_indicator--)
        {
            for(int alpha_2_indicator = 2; alpha_2_indicator >= 0; alpha_2_indicator--)
            {
                for(long s = s_max; s > s_min; s--)
                {
                    if(storage[t_indicator][alpha_2_indicator][s].size() < 10)
                        continue;

                    level_change = levelChange(alpha_1, alpha_2_min[alpha_2_indicator][t_indicator][s], t_indicator);
                    if(level_change <= 0)
                        continue;

                    new_alpha_2 = alpha_2_min[alpha_2_indicator][t_indicator][s] / alpha_1;
                    new_alpha_2_indicator = this->alpha_2_indicator(new_alpha_2);

                    if(s + level_change > pruningLevel - 11)
                    {   // return stages
                        stageCounterByLevel[s] -= storage[t_indicator][alpha_2_indicator][s].size();
                        stageCounterTotal -= storage[t_indicator][alpha_2_indicator][s].size();
                        returnStages(storage[t_indicator][alpha_2_indicator][s]);
                        alpha_2_min[t_indicator][alpha_2_indicator][s] = 2.0;
                    }
                    else
                    {   // move stages
                        stageCounterByLevel[s] -= storage[t_indicator][alpha_2_indicator][s].size();
                        stageCounterByLevel[s + level_change] += storage[t_indicator][alpha_2_indicator][s].size();
                        storage[t_indicator][new_alpha_2_indicator][s + level_change].splice(
                            storage[t_indicator][new_alpha_2_indicator][s + level_change].end(),
                            storage[t_indicator][alpha_2_indicator][s]);
                        alpha_2_min[t_indicator][alpha_2_indicator][s + level_change] = std::min(
                            alpha_2_min[t_indicator][alpha_2_indicator][s + level_change],
                            alpha_2_min[t_indicator][alpha_2_indicator][s]);
                        alpha_2_min[t_indicator][alpha_2_indicator][s] = 2.0;
                    }
                }
            }
        }
    }

    inline unsigned long levelChange(const double alpha_1, const double alpha_2, const long t_indicator) const
    {
        return static_cast<unsigned long>(t_indicator == 0
                                          ? 0
                                          : (alpha_1 <= alpha_2
                                             ? pruningLevel - min_level -
                                               1        // 2*log(2)= 1.386294361119890..., 39=40-1
                                             : ceil(
                    (t_indicator == 2 ? 39.0 : 17.0) * 1.38629436111989 / log((1.0 - alpha_2) / (alpha_1 - alpha_2))) -
                                               1));
    }

    inline void returnStages(list<NewEnumStage*>& stages)
    {
        pool.splice(pool.end(), stages);
    }
};

class NewEnum
{

private:
    const FactoringSettings settings;

    // Output
    Timer& timer;
    unsigned long round;

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
    Vec<double> log_V_minus_log_R_prod_minus_log_t;     /**< Contains the values \f$log(V_t / (r_{1,1} * ... * r_{t,t}))\f$ */
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


    void prepare(unsigned long round, const Mat<ZZ>& newBasis_transposed, const Mat<ZZ>& newU_scaled,
                 const Vec<RR>& new_target_coordinates)
    {
        this->round = round;
        conv(U_scaled_RR, newU_scaled);
        tau = new_target_coordinates;

        current_level = min_level;
        stagesCheckedForEquations = 0;

        ComputeGS(newBasis_transposed, mu, R_ii_squared);

        equations.clear();

        // setting the decreasing behavior
        decrease_max_distance = true;

        precomputeLogV();

        // start algorithm with a start parameter A
        A_curr = 0;
        for(long i = 1; i <= R_ii_squared.length(); i++)
        {
            A_curr += R_ii_squared(i);
        }

        A_curr *= 0.25;
        theoreticalMaxDistance = A_curr;
        A_curr *= settings.A_start_factor;
        heuristicMaxDistance = A_curr;

        L.resetStorage(A_curr);
    }

    /**
     * returns the closest integer to x -> \f$\lceil x \rfloor\f$
     * runs in \f$\mathcal{O}(1)\f$
     * @param x real number which to round
     * @return closest integer
     */
    inline void closest_RR(RR& out, const RR& x) const
    {
        sub(out, x, 0.5);
        ceil(out, out);
    }

    /**
     * Returns the smallest integer to y with \f$|u-y| < |next(u,y) - y|\f$
     * runs in \f$\mathcal{O}(1)\f$
     * @param u integer with minimal distance
     * @param y the 'center'
     * @return the next integer
     */
    inline void next(RR& out, const RR& u, const RR& y) const
    {
        RR temp;
        closest_RR(temp, y);

        auto side1 = static_cast<float>(temp >= y) * 2 - 1.0f;
//        float side1 = temp >= y ? 1.0f : -1.0f;
        auto side2 = static_cast<float>(u >= y) * 2 - 1.0f;
//        float side2 = u >= y ? 1.0f : -1.0f;

        // out = 2 * closest_y - u;
        mul(temp, temp, 2.0);
        sub(out, temp, u);

        out -= side2 * static_cast<float>(side1 == side2);
//        if(side1 == side2)
//            out -= side2;
    }

    /**
     * performs all stages starting with a start stage
     * @param start The stage we are starting with
     * @param perform_delayed_stages If true, the method will also perform delayed
     * stages. Only the first call of this method should set this parameter to true.
     */
    inline void perform(NewEnumStage* current_stage)
    {
        decrease_max_distance = true;

        long max_t = current_stage->t;
        long t = current_stage->t;          // projection to the n-t-1 last coordinates
        Vec<RR> c,                          // Length in the projection
            u,                          // coordinates of a close vector
            y;                          // uncovered parts of the target vector
        c.SetLength(t + 1);
        c(t + 1) = current_stage->c_tp1;
        current_stage->get_u(u);
        y.SetLength(t);
        y(t) = current_stage->y_t;

        long level;

        bool success;

        RR temp;

        // perform stages with s = current_level
        while(t <= max_t)
        {
            // c(t) = c(t+1) + power(abs(u(t) - y(t)), 2) * R_ii_squared(t);
            sub(c(t), u(t), y(t));
            abs(c(t), c(t));
            power(c(t), c(t), 2);
            c(t) *= R_ii_squared(t);
            c(t) += c(t + 1);

            if(c(t) >= A_curr)
            {
                // 2.1
                goto cleanup;
            }

            if(t == 1)
            {
                success = checkForEquation(u, c(1));

                if(c(1) < A_curr)
                {
                    if(decrease_max_distance || c(1) / A_curr < min_reduce_ratio)
                    {
                        A_curr = c(1);        // reduce max distance
                        L.updateMaxDistance(A_curr);
                    }
                }

                decrease_max_distance = decrease_max_distance && !success;

                // 2.1
                goto cleanup;
            }

            // compute the probability beta
            if(t < 4 || (level = calculateLevel(t, c(t))) <= current_level)
            {
                t--;

                clear(y(t));        // y(t) = 0
                for(long i = t + 1; i <= n; i++)    // 1/r_tt * sum_{i=t+1}^n (\tau_i - u_i) r_ti
                {
                    // y(t) += (tau(i) - u(i)) * mu(i, t);
                    sub(temp, tau(i), u(i));
                    mul(temp, temp, mu(i, t));
                    y(t) += temp;
                }
                y(t) += tau(t);

                closest_RR(u(t), y(t));
                continue;
            }

            // the program comes only this far, if level > current_s holds
            if(level <= max_level)
            {
                if(L.storeStage(y(t), c(t), c(t + 1), u, t, level))      // store the stage for later
                    return;
            }

            // 2.1
            cleanup:
            t++;
            if(t > max_t)
                break;
            next(u(t), u(t), y(t));
        }
    }


    /**
     * Precomputes some data used for the volume heuristic
     * runs in \f$\mathcal{O}(n)\f$
     * @param dim The dimension n of the lattice
     */
    static Vec<double> precomputeVolumes(const long n)
    {
        Vec<RR> V;
        Vec<double> log_V;
        V.SetLength(n);
        log_V.SetLength(n);

        RR pi = ComputePi_RR();
        V.SetLength(n);
        Vec<RR> R_diag_prod;
        R_diag_prod.SetLength(n);

        if(n >= 1)
            V(1) = conv<RR>(2);
        if(n >= 2)
            V(2) = pi;

        for(long i = 3; i <= n; i++)
        {
            V(i) = V(i - 2) * pi / (i / 2.0);
            log_V(i) = log(conv<double>(V(i)));
        }

        return log_V;
    }

    void precomputeLogV()
    {
        Vec<double> log_R_diag_prod;
        conv(log_R_diag_prod, R_ii_squared);

        for(long i = 1; i <= n; i++)
        {
            log_R_diag_prod(i) = log(log_R_diag_prod(i)) / 2.0;
        }

        for(long i = 2; i <= n; i++)
        {
            log_R_diag_prod(i) += log_R_diag_prod(i - 1);
        }

        log_V_minus_log_R_prod_minus_log_t = log_V;

        for(long i = 1; i <= n; i++)
        {
            log_V_minus_log_R_prod_minus_log_t(i) -= log_R_diag_prod(i);
            log_V_minus_log_R_prod_minus_log_t(i) -= log_t(i + 1);
        }
    }

    static Vec<double> precomputeLogT(const long n)
    {
        Vec<double> log_t;
        log_t.SetLength(n + 1);
        for(long t = 1; t <= n + 1; t++)
            log_t(t) = log(t);

        return log_t;
    }

    /**
     * computes (if possible) an equation from a given close vector
     * @param [in] input Coordinates of a close vector
     * @param [in] c_1 Distance of this vector to the target vector
     * @return true if an equation was found, else false.
     */
    bool checkForEquation(const Vec<RR>& input, const RR& c_1)
    {
        stagesCheckedForEquations++;
        mul(temp_vec, U_scaled_RR, input);
        temp_vec += shift;
        mul(close_vec, U_RR, temp_vec);

        conv(close_vec_long, close_vec);

        raw_equation.SetLength(n + 1);     // exponents of the first primes and -1
        equation.SetLength(n + 1);         // exponents of the first primes and -1

        NTL::set(u);       // u = 1

        ZZ temp_ZZ, temp_ZZ_2;
        RR temp_RR;

        for(long i = 1; i <= n; i++)
        {
            if(close_vec(i) > 0)
            {
                power(temp_ZZ, primes(i), close_vec_long(i));
                u *= temp_ZZ;
                raw_equation(i) = close_vec_long(i);
            }
            else
            {
                raw_equation(i) = 0;
            }
        }

        conv(v, u);
        conv(temp_RR, N);
        v /= temp_RR;

        closest_RR(temp_RR, v);

        sub(d, v, temp_RR);

        NTL::abs(alpha_nm1, d);

        closest_RR(temp_RR, v);
        conv(vN, temp_RR);
        vN *= N;

        NTL::clear(h_n);
        NTL::set(h_nm1);
        NTL::clear(h_nm2);
        NTL::set(k_n);
        NTL::clear(k_nm1);
        NTL::set(k_nm2);
#ifndef FS_CCF
        NTL::clear(a_nm1);
#else
        if(alpha_nm1 >= 0.5)
        NTL::set(a_nm1);
    else
        NTL::clear(a_nm1);

    ZZ h_n_abs, k_n_abs;
#endif

        long sign = NTL::sign(d),
            equation_counter = 0;
        bool cf_equation = false;

        do
        {
            nextContinuedFraction(h_n, k_n, h_nm1, k_nm1, h_nm2, k_nm2,
                                  a_nm1, alpha_nm1);

            equation = raw_equation;

#ifndef FS_CCF
//        left_side = u * k_n;
            mul(left_side, u, k_n);
//        abs(right_side,left_side - vN * k_n - sign * h_n * N);
            mul(temp_ZZ, vN, k_n);
            sub(temp_ZZ, left_side, temp_ZZ);
            mul(temp_ZZ_2, sign, h_n);
            mul(temp_ZZ_2, temp_ZZ_2, N);
            sub(right_side, temp_ZZ, temp_ZZ_2);
            abs(right_side, right_side);
#else
            //        left_side = u * k_n;
        abs(k_n_abs,k_n);
        abs(h_n_abs,h_n);
        mul(left_side,u,k_n_abs);
//        abs(right_side,left_side - vN * k_n - sign * h_n * N);
        mul(temp_ZZ,vN,k_n_abs);
        sub(temp_ZZ,left_side,temp_ZZ);
        mul(temp_ZZ_2,sign,h_n_abs);
        mul(temp_ZZ_2,temp_ZZ_2,N);
        sub(right_side,temp_ZZ,temp_ZZ_2);
        abs(right_side,right_side);
#endif

#ifndef FS_CCF
            if(isSmooth(equation, k_n, left_side, right_side))
#else
                if(isSmooth(equation, k_n_abs, left_side, right_side))
#endif
            {
                // ride_side_factor = conv<ZZ>(closest_RR(conv<RR>(left_side) / N_RR));
                equation_counter++;
                conv(temp_RR, left_side);
                div(temp_RR, temp_RR, N_RR);
                closest_RR(temp_RR, temp_RR);
                conv(ride_side_factor, temp_RR);
                // temp_RR = c_1 / A_max
                div(temp_RR, c_1, theoreticalMaxDistance);
                equations.emplace_back(equation, ride_side_factor, current_level, conv<double>(temp_RR),
                                       round, timer.step(), cf_equation);
            }

            cf_equation = true;
        }
#ifndef FS_CCF
        while(k_n < threshold && settings.useContinuedFractions);
#else
        while(k_n_abs < threshold && settings.useContinuedFractions);
#endif

        return equation_counter > 0;
    }

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
    void
    nextContinuedFraction(ZZ& h_n, ZZ& k_n, ZZ& h_nm1, ZZ& k_nm1, ZZ& h_nm2, ZZ& k_nm2, ZZ& a_nm1, RR& alpha_nm1) const
    {
//    h_n = a_nm1 * h_nm1 + h_nm2;
        mul(h_n, a_nm1, h_nm1);
        add(h_n, h_n, h_nm2);

//    k_n = a_nm1 * k_nm1 + k_nm2;
        mul(k_n, a_nm1, k_nm1);
        add(k_n, k_n, k_nm2);

        swap(h_nm2, h_nm1);
        h_nm1 = h_n;

        swap(k_nm2, k_nm1);
        k_nm1 = k_n;

#ifndef FS_CCF
        sub(alpha_nm1, alpha_nm1, floor(alpha_nm1));
#else
        sub(alpha_nm1,alpha_nm1,NTL::round(alpha_nm1));
#endif
        inv(alpha_nm1, alpha_nm1);                           // this is actually alpha_n
#ifndef FS_CCF
        FloorToZZ(a_nm1, alpha_nm1);                         // this is actually a_n
#else
        RoundToZZ(a_nm1,alpha_nm1);                         // this is actually a_n
#endif
    }

    /**
     * Checks if k_n and left_side are p_n-smooth and stores the prime factors in equation.
     * @return If k_n and the left side are smooth true is returned, else false.
     */
    bool isSmooth(Vec<long>& equation, ZZ& k_n, ZZ& left_side, ZZ& right_side) const
    {
        // check smoothness of k_n
        for(long i = 1; i <= n; i++)
        {
            while(k_n % primes(i) == 0)
            {
                equation(i) += 1;
                k_n /= primes(i);
            }
        }

        if(k_n > 1)        // if not smooth
            return false;

        if(right_side < 0)
            equation(n + 1) = 1;
        else
            equation(n + 1) = 0;

        ZZ right_side_abs = abs(right_side);

        // check smoothness of the right side
        for(long i = 1; i <= n; i++)
        {
            while(right_side_abs % primes(i) == 0)
            {
                if(equation(i) > 0)       // cancel
                    left_side /= primes(i);

                equation(i) -= 1;
                right_side_abs /= primes(i);
            }
        }

        return !(right_side_abs > 1 || left_side == 1);     // if not smooth or the equation is 1 = 1
    }

    inline long calculateLevel(long t, const RR& c_t) const
    {
        return conv<long>(floor(
            -((t - 1) / 2.0 * log(conv<double>(A_curr - c_t)) + log_V_minus_log_R_prod_minus_log_t(t - 1)) /
            0.693147180559945));    // log(2) = 0.69314...
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
    NewEnum(const FactoringSettings& settings, Timer& timer, const Vec<long>& primes,
            const Mat<ZZ>& U, Vec<RR> target_shift)
        : settings(settings), timer(timer), round(0), shift(std::move(target_shift)),
          decrease_max_distance(true), min_reduce_ratio(settings.reduce_ratio),
          U_RR(conv<Mat<RR>>(U)), n(primes.length()), N(settings.N), N_RR(conv<RR>(settings.N)),
          primes(primes), log_V(NewEnum::precomputeVolumes(n)), log_t(NewEnum::precomputeLogT(n)),
          max_level(settings.max_level), threshold(power_ZZ(primes(n), 3)),
          L(StageStorage(static_cast<unsigned long>(settings.max_level)))
    {
        // setting up the checkForEquation workspace
        raw_equation.SetLength(n + 1);     // exponents of the first primes and -1
        equation.SetLength(n + 1);         // exponents of the first primes and -1
    }

    /**
     * Returns a list of equations that were found in this round.
     */
    const list<Equation>& getEquations() const
    {
        return equations;
    }

    void getDistances(RR& theoretical, RR& heuristic, RR& reduced)
    {
        theoretical = theoreticalMaxDistance;
        heuristic = heuristicMaxDistance;
        reduced = A_curr;
    }

    /**
     * starts the performing by setting the first stage to the values that are
     * described under point 1 in the NewEnum algorithm.
     */
    void run(unsigned long round, const Mat<ZZ>& newBasis_transposed, const Mat<ZZ>& new_U_scaled,
             const Vec<RR>& new_target_coordinates)
    {
        prepare(round, newBasis_transposed, new_U_scaled, new_target_coordinates);

        cout << "Performing: ";
        // Reset list of delayed stages
        current_level = min_level;

        Vec<RR> u;                              // coordinates of a close vector
        u.SetLength(n);
        closest_RR(u(n), tau(n));

        NewEnumStage* stage = new NewEnumStage(tau(n), conv<RR>(0), u, n);

        perform(stage);
        L.returnStage(stage);

        for(long l = 0; l <= max_level - min_level - 1; l++)
        {
            L.incrementCurrentLevel();
            current_level++;

            // perform delayed stages
            while(L.getNext(stage))
            {
                perform(stage);
                L.returnStage(stage);
            }

            cout << ".";
        }

        cout << "\n";
    }
};

#endif    /* NEWENUM_H */