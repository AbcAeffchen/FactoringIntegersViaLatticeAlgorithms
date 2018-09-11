//
// Created by Alex Schickedanz <alex@ae.cs.uni-frankfurt.de> on 11.09.18.
//

#ifndef FACTORINGINTEGERSVIALATTICEALGORITHMS_STAGESTORAGE_H
#define FACTORINGINTEGERSVIALATTICEALGORITHMS_STAGESTORAGE_H

#include "Equation.h"
#include "Statistics.h"
#include "NewEnumStage.h"

#include <vector>
#include <list>

class StageStorage
{
private:
    const double alpha_1_threshold = 0.99;             /**< A_new/A_old = alpha_1 < alpha_1_threshold.
                                                            Used to prevent to much useless level recalculations. */
    const unsigned long long stageCounterThreshold = 2500000;   // About 6 GB of RAM depending on the dimension
    const long unsigned minLevel = 10;                  /**< minimum level */
    const long unsigned numLevel;                       /**< number of levels stages are stored for */

    unsigned long long stageCounterTotal = 0;
    std::vector<unsigned long long> stageCounterByLevel;

    unsigned long currentLevel = 10;

    RR maxDistance;

    std::vector<double> alpha_2_min;     /**< organized as alpha_2_min[alpha_2_indicator][t_indicator][level] */

    std::list<NewEnumStage*> pool;
    std::vector<std::list<NewEnumStage*>> storage;        /**< A queue of stages for every level s and projection t
                                                               organized as storage[t_indicator][alpha_2_indicator][s-11]. */

public:
    // statistics
    /**
     * Contains the number of stages stored for the level
     * at the moment this level gets performed.
     * Organized as [alpha_2_indicator][t_indicator][level-11]
     */
    std::vector<unsigned long long> maxDelayedAndPerformedStages;
    unsigned long long totalDelayedAndPerformedStages = 0;
    /**
     * Contains the number of stages stored for the level (performed and deleted)
     * Organized as [alpha_2_indicator][t_indicator][level-11]
     */
    std::vector<unsigned long long> delayedStages;

    template<bool precalculatedLevel = false>
    inline size_t coordinatesToIndex(unsigned long t_indicator, unsigned long alpha_2_indicator, unsigned long level) const noexcept
    {
        if(precalculatedLevel)
            return alpha_2_indicator + 3 * t_indicator + 9 * level;
        else
            return alpha_2_indicator + 3 * t_indicator + 9 * (level - minLevel - 1);
    }

    explicit StageStorage(unsigned long pruningLevel)
        : numLevel(pruningLevel - minLevel),
          stageCounterByLevel(std::vector<unsigned long long>(numLevel, 0)),
          alpha_2_min(std::vector<double>(9 * numLevel, 2.0)),
          storage(std::vector<std::list<NewEnumStage*>> (9 * numLevel)),
          maxDelayedAndPerformedStages(std::vector<unsigned long long>(9 * numLevel, 0)),
          delayedStages(std::vector<unsigned long long> (9 * numLevel, 0))
    {
        // allocate 200000 stages for the beginning
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
        if(stageCounterByLevel[currentLevel - minLevel - 1] <= 0 || stageCounterTotal >= stageCounterThreshold)
            return false;

        for(unsigned long t_indicator = 0; t_indicator < 3; ++t_indicator)
        {
            for(unsigned long alpha_2_indicator = 0; alpha_2_indicator < 3; ++alpha_2_indicator)
            {
                if(!storage[coordinatesToIndex(t_indicator, alpha_2_indicator, currentLevel)].empty())
                {
                    stage = storage[coordinatesToIndex(t_indicator, alpha_2_indicator, currentLevel)].front();
                    storage[coordinatesToIndex(t_indicator, alpha_2_indicator, currentLevel)].pop_front();
                    stageCounterByLevel[currentLevel - minLevel - 1]--;
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
                maxDelayedAndPerformedStages[coordinatesToIndex(t_indicator, alpha_2_indicator, currentLevel)]
                    = storage[coordinatesToIndex(t_indicator, alpha_2_indicator, currentLevel)].size();
                totalDelayedAndPerformedStages += storage[coordinatesToIndex(t_indicator, alpha_2_indicator, currentLevel)].size();
            }
    }

    bool storeStage(const RR& y_t, const RR& c_t, const RR& c_tp1, const Vec<RR>& u, long t, unsigned long level)
    {
        NewEnumStage* stage = getStage();
        stage->set(y_t, c_tp1, u, t);

        double alpha_2;
        conv(alpha_2, c_t / maxDistance);
        unsigned long alpha_2_indicator = this->alpha_2_indicator(alpha_2);
        unsigned long t_indicator = this->t_indicator(stage->t);
        storage[coordinatesToIndex(t_indicator, alpha_2_indicator, level)].push_back(stage);

        if(alpha_2_min[coordinatesToIndex(t_indicator, alpha_2_indicator, level)] > alpha_2)
            alpha_2_min[coordinatesToIndex(t_indicator, alpha_2_indicator, level)] = alpha_2;

        stageCounterTotal++;
        stageCounterByLevel[level - minLevel - 1]++;

        // statistics
        delayedStages[coordinatesToIndex(t_indicator, alpha_2_indicator, level)]++;

        return stageCounterTotal >= stageCounterThreshold;      // about 6 GB of RAM depending on the dimension
    }

    void updateMaxDistance(const RR& distance)
    {
        double alpha_1;
        conv(alpha_1, distance / maxDistance);

        if(alpha_1 < alpha_1_threshold && stageCounterTotal > 100)
            recalculateLevels(alpha_1);

        for(unsigned long alpha_2_indicator = 0; alpha_2_indicator < 3; alpha_2_indicator++)
            for(unsigned long t_indicator = 0; t_indicator < 3; t_indicator++)
                for(unsigned long s = currentLevel - minLevel - 1; s < numLevel; s++)
                    alpha_2_min[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)] /= alpha_1;

        maxDistance = distance;
    }

    void resetStorage(const RR& A = conv<RR>(0))
    {
        // reset counters/statistics
        totalDelayedAndPerformedStages = 0;
        stageCounterTotal = 0;

        std::fill(stageCounterByLevel.begin(), stageCounterByLevel.end(), 0);
        std::fill(maxDelayedAndPerformedStages.begin(), maxDelayedAndPerformedStages.end(), 0);
        std::fill(delayedStages.begin(), delayedStages.end(), 0);
        std::fill(alpha_2_min.begin(), alpha_2_min.end(), 2.0);

        for(auto& stageList : storage)
        {
            returnStages(stageList);
        }

        maxDistance = A;
        currentLevel = minLevel;
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
    inline unsigned long alpha_2_indicator(const double alpha_2) const
    {
//        const bool a = alpha_2 > 0.65;
        // return static_cast<long>(a) * 2 + static_cast<long>(!a) * static_cast<long>(alpha_2 > 0.4);
        return (alpha_2 > 0.65) ? 2 : (alpha_2 > 0.4 ? 1 : 0);
    }

    inline unsigned long t_indicator(const long t) const
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
        unsigned long new_alpha_2_indicator;
        unsigned long s_max = numLevel - 1;
        unsigned long s_min = currentLevel - minLevel - 1;

        /**
         * skip t_indicator = 0, since the stages here have t between 4 and 20.
         * So in almost all cases there will be no change.
         */
        for(unsigned long t_indicator = 2; t_indicator >= 1; t_indicator--)
        {
            for(unsigned long alpha_2_indicator = 3; alpha_2_indicator-- > 0; /*alpha_2_indicator--*/)
            {
                for(unsigned long s = s_max; s > s_min; s--)
                {
                    if(storage[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)].size() < 10)
                        continue;

                    const auto level_change = levelChange(alpha_1, alpha_2_min[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)], t_indicator);
                    if(level_change <= 0)
                        continue;

                    new_alpha_2 = alpha_2_min[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)] / alpha_1;
                    new_alpha_2_indicator = this->alpha_2_indicator(new_alpha_2);

                    if(s + level_change > numLevel - 1)
                    {   // return stages
                        stageCounterByLevel[s] -= storage[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)].size();
                        stageCounterTotal -= storage[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)].size();
                        returnStages(storage[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)]);
                        alpha_2_min[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)] = 2.0;
                    }
                    else
                    {   // move stages
                        stageCounterByLevel[s] -= storage[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)].size();
                        stageCounterByLevel[s + level_change] += storage[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)].size();
                        storage[coordinatesToIndex<true>(t_indicator, new_alpha_2_indicator, s + level_change)].splice(
                            storage[coordinatesToIndex<true>(t_indicator, new_alpha_2_indicator, s + level_change)].end(),
                            storage[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)]);
                        alpha_2_min[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s + level_change)] = std::min(
                            alpha_2_min[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s + level_change)],
                            alpha_2_min[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)]);
                        alpha_2_min[coordinatesToIndex<true>(t_indicator, alpha_2_indicator, s)] = 2.0;
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
                                             ? numLevel - 1        // 2*log(2)= 1.386294361119890..., 39=40-1
                                             : ceil(
                    (t_indicator == 2 ? 39.0 : 17.0) * 1.38629436111989 / log((1.0 - alpha_2) / (alpha_1 - alpha_2))) -
                                               1));
    }

    inline void returnStages(std::list<NewEnumStage*>& stages)
    {
        pool.splice(pool.end(), stages);
    }
};

#endif //FACTORINGINTEGERSVIALATTICEALGORITHMS_STAGESTORAGE_H
