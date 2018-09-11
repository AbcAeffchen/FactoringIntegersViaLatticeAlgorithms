
#ifndef FACTORINGINTEGERSVIALATTICEALGORITHMS_FACTORINGSETTINGS_H
#define FACTORINGINTEGERSVIALATTICEALGORITHMS_FACTORINGSETTINGS_H

#include <NTL/ZZ.h>
#include <NTL/RR.h>

#include <cmath>
#include <utility>

#define FS_SCALING_ONE_HALF 1
#define FS_SCALING_ONE_FOURTH 2
#define FS_SCALING_THREE_FORTH 3
#define FS_SCALING_ONE_FOURTH_ONE_HALF 4
#define FS_SCALING_ONE_HALF_ONE_FOURTH 5
#define FS_SCALING_MIXED 0

using namespace NTL;

struct FactoringSettings
{
    const ZZ N;
    const unsigned long n;
    const RR c;
    const unsigned int max_level = 14;
    const double A_start_factor = 0.2;
    const double reduce_ratio = 0.8;
    const long accuracy_factor = 10000;
    const long strong_bkz = 32;
    const long slight_bkz = 20;
    const unsigned long min_eqns;
    const long long int seed_type = -2;
    const bool useContinuedFractions = true;
    const int scalingType = 0;

    FactoringSettings(ZZ N, unsigned long n, RR c)
        : N(std::move(N)), n(n), c(std::move(c)), min_eqns(n+1)
    { }

    FactoringSettings(ZZ N, unsigned long n, RR c, unsigned long min_eqns)
        : N(std::move(N)), n(n), c(std::move(c)), min_eqns(min_eqns)
    { }

    FactoringSettings(ZZ N, unsigned long n, RR c, unsigned int max_level,
                      double A_start_factor, double reduce_ratio, long accuracy_factor,
                      long strong_bkz, long slight_bkz, unsigned long min_eqns,
                      long long int seed_type = -2, bool useContinuedFractions = true, int scalingType = 0)
        : N(std::move(N)), n(n), c(std::move(c)), max_level(std::max(10u, max_level)), A_start_factor(A_start_factor),
          reduce_ratio(reduce_ratio > 0 ? reduce_ratio : 0), accuracy_factor(accuracy_factor),
          strong_bkz(strong_bkz), slight_bkz(slight_bkz), min_eqns(min_eqns), seed_type(seed_type),
          useContinuedFractions(useContinuedFractions),scalingType(scalingType)
    { }
};

#endif //FACTORINGINTEGERSVIALATTICEALGORITHMS_FACTORINGSETTINGS_H