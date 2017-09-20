
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

    NewEnumStage() : y_t(conv<RR>(0)), c_tp1(conv<RR>(0)), t(0) {}

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
            : min_level(10), pruningLevel(pruningLevel),
              stageCounterByLevel(vector<unsigned long long>(pruningLevel - min_level,0))
    {
        this->storage = std::vector<std::vector<std::vector<std::list<NewEnumStage*>>>>(3,std::vector<std::vector<std::list<NewEnumStage*>>>(3, std::vector<std::list<NewEnumStage*>>(pruningLevel - min_level)));
        this->maxDelayedAndPerformedStages = std::vector<std::vector<std::vector<unsigned long long>>>(3, std::vector<std::vector<unsigned long long>>(3, std::vector<unsigned long long>(pruningLevel - min_level,0)));
        this->delayedStages = std::vector<std::vector<std::vector<unsigned long long>>>(3, std::vector<std::vector<unsigned long long>>(3, std::vector<unsigned long long>(pruningLevel - min_level,0)));
        this->alpha_2_min = std::vector<std::vector<std::vector<double>>>(3,std::vector<std::vector<double>>(3, std::vector<double>(pruningLevel-min_level,2.0)));
        // allocate 1000 stages for the beginning
        for(int i = 0; i < 200000; ++i)
            this->pool.push_back(new NewEnumStage);
    }

    ~StageStorage()
    {
        this->resetStorage();   // return all stages to the pool
        NewEnumStage* temp;
        // delete all stages
        while(!this->pool.empty())
        {
            temp = this->pool.front();
            this->pool.pop_front();
            delete temp;
        }
    }

    bool getNext(NewEnumStage* &stage)
    {
        if(this->stageCounterByLevel[this->currentLevel - this->min_level - 1] <= 0 || this->stageCounterTotal >= this->stageCounterThreshold)
            return false;

        for(unsigned long t_indicator = 0; t_indicator < 3; ++t_indicator)
        {
            for(unsigned long alpha_2_indicator = 0; alpha_2_indicator < 3; ++alpha_2_indicator)
            {
                if(!this->storage[t_indicator][alpha_2_indicator][this->currentLevel - this->min_level - 1].empty())
                {
                    stage = this->storage[t_indicator][alpha_2_indicator][this->currentLevel - this->min_level - 1].front();
                    this->storage[t_indicator][alpha_2_indicator][this->currentLevel - this->min_level - 1].pop_front();
                    this->stageCounterByLevel[this->currentLevel - this->min_level-1]--;
                    this->stageCounterTotal--;
                    return true;
                }
            }
        }

        // It never should come this far, but who knows...
        return false;
    }

    void incrementCurrentLevel()
    {
        this->currentLevel++;

        for(long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
            for(long t_indicator = 0; t_indicator < 3; t_indicator++)
            {
                this->maxDelayedAndPerformedStages[alpha_2_indicator][t_indicator][this->currentLevel - this->min_level - 1] = this->storage[t_indicator][alpha_2_indicator][this->currentLevel - this->min_level - 1].size();
                this->totalDelayedAndPerformedStages += this->storage[t_indicator][alpha_2_indicator][this->currentLevel - this->min_level - 1].size();
            }
    }


    bool storeStage(const RR& y_t, const RR& c_t, const RR& c_tp1, const Vec<RR>& u, long t, long level)
    {
        NewEnumStage* stage = this->getStage();
        stage->set(y_t, c_tp1, u, t);

        double alpha_2;
        conv(alpha_2, c_t / this->maxDistance);
        long alpha_2_indicator = this->alpha_2_indicator(alpha_2);
        long t_indicator = this->t_indicator(stage->t);
        this->storage[t_indicator][alpha_2_indicator][level - this->min_level - 1].push_back(stage);
        if(this->alpha_2_min[alpha_2_indicator][t_indicator][level - this->min_level - 1] > alpha_2)
            this->alpha_2_min[alpha_2_indicator][t_indicator][level - this->min_level - 1] = alpha_2;
        this->stageCounterTotal++;
        this->stageCounterByLevel[level - this->min_level - 1]++;

        // statistics
        this->delayedStages[alpha_2_indicator][t_indicator][level - this->min_level - 1]++;

        return this->stageCounterTotal >= this->stageCounterThreshold;      // about 6 GB of RAM depending on the dimension
    }

    void updateMaxDistance(const RR& distance)
    {
        double alpha_1;
        conv(alpha_1, distance/this->maxDistance);

        if(alpha_1 < this->alpha_1_threshold && this->stageCounterTotal > 100)
            this->recalculateLevels(alpha_1);

        for(long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
            for(long t_indicator = 0; t_indicator < 3; t_indicator++)
                for(long s = std::max(0L,this->currentLevel-this->min_level-1); s < this->pruningLevel-this->min_level;s++)
                    this->alpha_2_min[alpha_2_indicator][t_indicator][s] /= alpha_1;

        this->maxDistance = distance;
    }

    void resetStorage(const RR& A = conv<RR>(0))
    {
        // reset counters/statistics
        this->totalDelayedAndPerformedStages = 0;
        this->stageCounterTotal = 0;
        for(int level = 0; level < this->pruningLevel - this->min_level; level++)
            this->stageCounterByLevel[level] = 0;

        for(long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
            for(long t_indicator = 0; t_indicator < 3; t_indicator++)
                for(int level = 0; level < this->pruningLevel - this->min_level; level++)
                {
                    this->maxDelayedAndPerformedStages[alpha_2_indicator][t_indicator][level] = 0;
                    this->delayedStages[alpha_2_indicator][t_indicator][level] = 0;
                    this->alpha_2_min[alpha_2_indicator][t_indicator][level] = 2.0;
                    this->returnStages(this->storage[t_indicator][alpha_2_indicator][level]);
                }

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

    unsigned long long getTotalDelayedStages()
    {
        return this->stageCounterTotal;
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
            for(long i = 0; i < 10000; i++)
                this->pool.push_back(new NewEnumStage);
        }

        NewEnumStage* stage = this->pool.front();
        this->pool.pop_front();
        return stage;
    }

    void recalculateLevels(const double alpha_1)
    {
        double new_alpha_2;
        long new_alpha_2_indicator;
        long level_change;
        long s_max = this->pruningLevel - this->min_level - 1;
        long s_min = this->currentLevel-this->min_level - 1;

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
                    if (this->storage[t_indicator][alpha_2_indicator][s].size() < 10)
                        continue;

                    level_change = this->levelChange(alpha_1,
                                                     this->alpha_2_min[alpha_2_indicator][t_indicator][s],
                                                     t_indicator);
                    if (level_change <= 0)
                        continue;

                    new_alpha_2 = this->alpha_2_min[alpha_2_indicator][t_indicator][s] / alpha_1;
                    new_alpha_2_indicator = this->alpha_2_indicator(new_alpha_2);

                    if(s+level_change > this->pruningLevel-11)
                    {   // return stages
                        this->stageCounterByLevel[s] -= this->storage[t_indicator][alpha_2_indicator][s].size();
                        this->stageCounterTotal -= this->storage[t_indicator][alpha_2_indicator][s].size();
                        this->returnStages(this->storage[t_indicator][alpha_2_indicator][s]);
                        this->alpha_2_min[t_indicator][alpha_2_indicator][s] = 2.0;
                    }
                    else
                    {   // move stages
                        this->stageCounterByLevel[s] -= this->storage[t_indicator][alpha_2_indicator][s].size();
                        this->stageCounterByLevel[s + level_change] += this->storage[t_indicator][alpha_2_indicator][s].size();
                        this->storage[t_indicator][new_alpha_2_indicator][s + level_change].splice(
                                this->storage[t_indicator][new_alpha_2_indicator][s + level_change].end(),
                                this->storage[t_indicator][alpha_2_indicator][s]);
                        this->alpha_2_min[t_indicator][alpha_2_indicator][s + level_change] = std::min(this->alpha_2_min[t_indicator][alpha_2_indicator][s + level_change],
                                                                                                       this->alpha_2_min[t_indicator][alpha_2_indicator][s]);
                        this->alpha_2_min[t_indicator][alpha_2_indicator][s] = 2.0;
                    }
                }
            }
        }
    }

    inline unsigned long levelChange(const double alpha_1, const double alpha_2, const long t_indicator) const
    {
        return t_indicator == 0
               ? 0
               : (alpha_1 <= alpha_2
                  ? (unsigned long) this->pruningLevel - this->min_level - 1        // 2*log(2)= 1.386294361119890..., 39=40-1
                  : (unsigned long) ceil((t_indicator == 2 ? 39.0 : 17.0)  * 1.38629436111989 / log((1.0-alpha_2) / (alpha_1 - alpha_2))) - 1);
    }


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


    void prepare(unsigned long round, const Mat<ZZ>& newBasis_transposed, const Mat<ZZ>& newU_scaled, const Vec<RR>& new_target_coordinates)
    {
        this->round = round;
        conv(this->U_scaled_RR,newU_scaled);
        this->tau = new_target_coordinates;

        this->current_level = this->min_level;
        this->stagesCheckedForEquations = 0;

        ComputeGS(newBasis_transposed, this->mu, this->R_ii_squared);

        this->equations.clear();

        // setting the decreasing behavior
        this->decrease_max_distance = true;

        this->precomputeLogV();

        // start algorithm with a start parameter A
        this->A_curr = 0;
        for(long i = 1; i <= this->R_ii_squared.length(); i++)
        {
            this->A_curr += this->R_ii_squared(i);
        }

        this->A_curr *= 0.25;
        this->theoreticalMaxDistance = this->A_curr;
        this->A_curr *= this->settings.A_start_factor;
        this->heuristicMaxDistance = this->A_curr;

        this->L.resetStorage(this->A_curr);
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
        this->closest_RR(temp, y);

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
        this->decrease_max_distance = true;

        long max_t = current_stage->t;
        long t = current_stage->t;          // projection to the n-t-1 last coordinates
        Vec<RR> c,                          // Length in the projection
                u,                          // coordinates of a close vector
                y;                          // uncovered parts of the target vector
        c.SetLength(t+1);
        c(t+1) = current_stage->c_tp1;
        current_stage->get_u(u);
        y.SetLength(t);
        y(t) = current_stage->y_t;

        long level;

        bool success;

        RR temp;

        // perform stages with s = current_level
        while(t <= max_t)
        {
            // c(t) = c(t+1) + power(abs(u(t) - y(t)), 2) * this->R_ii_squared(t);
            sub(c(t),u(t),y(t));
            abs(c(t),c(t));
            power(c(t),c(t),2);
            c(t) *= this->R_ii_squared(t);
            c(t) += c(t+1);

            if(c(t) >= this->A_curr)
            {
                // 2.1
                goto cleanup;
            }

            if(t == 1)
            {
                success = this->checkForEquation(u, c(1));

                if(c(1) < this->A_curr)
                {
                    if(this->decrease_max_distance || c(1) / this->A_curr < this->min_reduce_ratio)
                    {
                        this->A_curr = c(1);        // reduce max distance
                        this->L.updateMaxDistance(this->A_curr);
                    }
                }

                this->decrease_max_distance = this->decrease_max_distance && !success;

                // 2.1
                goto cleanup;
            }

            // compute the probability beta
            if(t < 4 || (level = calculateLevel(t,c(t))) <= this->current_level)
            {
                t--;

                clear(y(t));        // y(t) = 0
                for(long i = t + 1; i <= this->n; i++)    // 1/r_tt * sum_{i=t+1}^n (\tau_i - u_i) r_ti
                {
                    // y(t) += (this->tau(i) - u(i)) * this->mu(i, t);
                    sub(temp,this->tau(i),u(i));
                    mul(temp,temp,this->mu(i, t));
                    y(t) += temp;
                }
                y(t) += this->tau(t);

                this->closest_RR(u(t),y(t));
                continue;
            }

            // the program comes only this far, if level > current_s holds
            if(level <= this->max_level)
            {
                if(this->L.storeStage(y(t),c(t),c(t+1),u,t,level))      // store the stage for later
                    return;
            }

            // 2.1
            cleanup:
            t++;
            if(t > max_t)
                break;
            this->next(u(t), u(t), y(t));
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
            V(i) = V(i-2) * pi / (i/2.0);
            log_V(i) = log(conv<double>(V(i)));
        }

        return log_V;
    }

    void precomputeLogV()
    {
        Vec<double> log_R_diag_prod;
        conv(log_R_diag_prod, this->R_ii_squared);

        for(long i = 1; i <= this->n; i++)
        {
            log_R_diag_prod(i) = log(log_R_diag_prod(i)) / 2.0;
        }

        for(long i = 2; i <= this->n; i++)
        {
            log_R_diag_prod(i) += log_R_diag_prod(i - 1);
        }

        this->log_V_minus_log_R_prod_minus_log_t = this->log_V;

        for(long i = 1; i <= this->n; i++)
        {
            this->log_V_minus_log_R_prod_minus_log_t(i) -= log_R_diag_prod(i);
            this->log_V_minus_log_R_prod_minus_log_t(i) -= this->log_t(i + 1);
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
        this->stagesCheckedForEquations++;
        mul(this->temp_vec, this->U_scaled_RR, input);
        this->temp_vec += this->shift;
        mul(this->close_vec, this->U_RR, this->temp_vec);

        conv(this->close_vec_long, this->close_vec);

        this->raw_equation.SetLength(this->n + 1);     // exponents of the first primes and -1
        this->equation.SetLength(this->n + 1);         // exponents of the first primes and -1

        NTL::set(this->u);       // u = 1

        ZZ temp_ZZ, temp_ZZ_2;
        RR temp_RR;

        for(long i = 1; i <= this->n; i++)
        {
            if(this->close_vec(i) > 0)
            {
                power(temp_ZZ, this->primes(i), this->close_vec_long(i));
                this->u *= temp_ZZ;
                this->raw_equation(i) = this->close_vec_long(i);
            }
            else
            {
                this->raw_equation(i) = 0;
            }
        }

        conv(this->v, this->u);
        conv(temp_RR, this->N);
        this->v /= temp_RR;

        this->closest_RR(temp_RR, this->v);

        sub(this->d, this->v, temp_RR);

        NTL::abs(this->alpha_nm1, this->d);

        this->closest_RR(temp_RR, this->v);
        conv(this->vN, temp_RR);
        this->vN *= this->N;

        NTL::clear(this->h_n);
        NTL::set(this->h_nm1);
        NTL::clear(this->h_nm2);
        NTL::set(this->k_n);
        NTL::clear(this->k_nm1);
        NTL::set(this->k_nm2);
#ifndef FS_CCF
        NTL::clear(this->a_nm1);
#else
        if(this->alpha_nm1 >= 0.5)
        NTL::set(this->a_nm1);
    else
        NTL::clear(this->a_nm1);

    ZZ h_n_abs, k_n_abs;
#endif

        long sign = NTL::sign(this->d),
                equation_counter = 0;
        bool cf_equation = false;

        do
        {
            this->nextContinuedFraction(this->h_n, this->k_n, this->h_nm1, this->k_nm1, this->h_nm2, this->k_nm2,
                                        this->a_nm1, this->alpha_nm1);

            this->equation = this->raw_equation;

#ifndef FS_CCF
//        this->left_side = this->u * this->k_n;
            mul(this->left_side, this->u, this->k_n);
//        abs(this->right_side,this->left_side - this->vN * this->k_n - sign * this->h_n * this->N);
            mul(temp_ZZ, this->vN, this->k_n);
            sub(temp_ZZ, this->left_side, temp_ZZ);
            mul(temp_ZZ_2, sign, this->h_n);
            mul(temp_ZZ_2, temp_ZZ_2, this->N);
            sub(this->right_side, temp_ZZ, temp_ZZ_2);
            abs(this->right_side, this->right_side);
#else
            //        this->left_side = this->u * this->k_n;
        abs(k_n_abs,this->k_n);
        abs(h_n_abs,this->h_n);
        mul(this->left_side,this->u,k_n_abs);
//        abs(this->right_side,this->left_side - this->vN * this->k_n - sign * this->h_n * this->N);
        mul(temp_ZZ,this->vN,k_n_abs);
        sub(temp_ZZ,this->left_side,temp_ZZ);
        mul(temp_ZZ_2,sign,h_n_abs);
        mul(temp_ZZ_2,temp_ZZ_2,this->N);
        sub(this->right_side,temp_ZZ,temp_ZZ_2);
        abs(this->right_side,this->right_side);
#endif

#ifndef FS_CCF
            if(this->isSmooth(this->equation, this->k_n, this->left_side, this->right_side))
#else
                if(this->isSmooth(this->equation, k_n_abs, this->left_side, this->right_side))
#endif
            {
                // this->ride_side_factor = conv<ZZ>(this->closest_RR(conv<RR>(left_side) / this->N_RR));
                equation_counter++;
                conv(temp_RR, left_side);
                div(temp_RR, temp_RR, this->N_RR);
                this->closest_RR(temp_RR, temp_RR);
                conv(this->ride_side_factor, temp_RR);
                // temp_RR = c_1 / A_max
                div(temp_RR, c_1, this->theoreticalMaxDistance);
                this->equations.emplace_back(this->equation, this->ride_side_factor, this->current_level, conv<double>(temp_RR),
                                  this->round, this->timer.step(), cf_equation);
            }

            cf_equation = true;
        }
#ifndef FS_CCF
        while(this->k_n < this->threshold && this->settings.useContinuedFractions);
#else
        while(k_n_abs < this->threshold && this->settings.useContinuedFractions);
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
    void nextContinuedFraction(ZZ& h_n, ZZ& k_n, ZZ& h_nm1, ZZ& k_nm1, ZZ& h_nm2, ZZ& k_nm2, ZZ& a_nm1, RR& alpha_nm1) const
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
        for(long i = 1; i <= this->n; i++)
        {
            while(k_n % this->primes(i) == 0)
            {
                equation(i) += 1;
                k_n /= this->primes(i);
            }
        }

        if(k_n > 1)        // if not smooth
            return false;

        if(right_side < 0)
            equation(this->n + 1) = 1;
        else
            equation(this->n + 1) = 0;

        ZZ right_side_abs = abs(right_side);

        // check smoothness of the right side
        for(long i = 1; i <= this->n; i++)
        {
            while(right_side_abs % this->primes(i) == 0)
            {
                if(equation(i) > 0)       // cancel
                    left_side /= this->primes(i);

                equation(i) -= 1;
                right_side_abs /= this->primes(i);
            }
        }

        return !(right_side_abs > 1 || left_side == 1);     // if not smooth or the equation is 1 = 1
    }

    inline long calculateLevel(long t, const RR& c_t) const
    {
        return conv<long>(floor(-((t-1)/2.0*log(conv<double>(this->A_curr - c_t)) + this->log_V_minus_log_R_prod_minus_log_t(t - 1)) / 0.693147180559945));    // log(2) = 0.69314...
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
        : settings(settings), timer(timer), N(settings.N), N_RR(conv<RR>(settings.N)), U_RR(conv<Mat<RR>>(U)),
          primes(primes), shift(std::move(target_shift)), n(primes.length()), min_reduce_ratio(settings.reduce_ratio),
          max_level(settings.max_level), log_t(NewEnum::precomputeLogT(n)), threshold(power_ZZ(primes(n),3)),
          log_V(NewEnum::precomputeVolumes(n)), L(StageStorage((unsigned long) settings.max_level))
    {
        // setting up the checkForEquation workspace
        this->raw_equation.SetLength(this->n + 1);     // exponents of the first primes and -1
        this->equation.SetLength(this->n + 1);         // exponents of the first primes and -1
    }

    /**
     * Returns a list of equations that were found in this round.
     */
    const list<Equation>& getEquations() const
    {
        return this->equations;
    }

    void getDistances(RR& theoretical, RR& heuristic, RR& reduced)
    {
        theoretical = this->theoreticalMaxDistance;
        heuristic = this->heuristicMaxDistance;
        reduced = this->A_curr;
    }

    /**
     * starts the performing by setting the first stage to the values that are
     * described under point 1 in the NewEnum algorithm.
     */
    void run(unsigned long round, const Mat<ZZ>& newBasis_transposed, const Mat<ZZ>& new_U_scaled, const Vec<RR>& new_target_coordinates)
    {
        this->prepare(round, newBasis_transposed, new_U_scaled, new_target_coordinates);

        cout << "Performing: ";
        // Reset list of delayed stages
        this->current_level = this->min_level;

        Vec<RR> u;                              // coordinates of a close vector
        u.SetLength(this->n);
        this->closest_RR(u(this->n),this->tau(this->n));

        NewEnumStage* stage = new NewEnumStage(this->tau(this->n), conv<RR>(0), u, this->n);

        this->perform(stage);
        this->L.returnStage(stage);

        for(long l = 0; l <= this->max_level - this->min_level - 1; l++)
        {
            this->L.incrementCurrentLevel();
            this->current_level++;

            // perform delayed stages
            while(this->L.getNext(stage))
            {
                this->perform(stage);
                this->L.returnStage(stage);
            }

            cout << ".";
        }

        cout << endl;
    }

};

#endif	/* NEWENUM_H */