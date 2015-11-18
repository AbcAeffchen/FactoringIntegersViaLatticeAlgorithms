
#ifndef FACTORINGINTEGERSVIALATTICEALGORITHMS_FACTORINGSETTINGS_H
#define FACTORINGINTEGERSVIALATTICEALGORITHMS_FACTORINGSETTINGS_H

#include <NTL/ZZ.h>
#include <NTL/RR.h>

struct FactoringSettings
{
    const ZZ N;
    const unsigned long n;
    const RR c;
    const int max_level = 14;
    const double A_start_factor = 0.2;
    const double reduce_ratio = 0.8;
    const long accuracy_factor = 10000;
    const long strong_bkz = 32;
    const long slight_bkz = 20;
    const unsigned long min_eqns;
    const long long int seed_type = -2;
    const int scalingType = 1;

    FactoringSettings(const ZZ &N, unsigned long n, const RR &c)
            : N(N), n(n), c(c), min_eqns(n+1)
    { }

    FactoringSettings(const ZZ &N, unsigned long n, const RR &c, unsigned long min_eqns)
            : N(N), n(n), c(c), min_eqns(min_eqns)
    { }

    FactoringSettings(const ZZ &N, unsigned long n, const RR &c, int max_level,
                      double A_start_factor, double reduce_ratio, long accuracy_factor,
                      long strong_bkz, long slight_bkz, unsigned long min_eqns)
            : N(N), n(n), c(c), max_level(max(10, max_level)), A_start_factor(A_start_factor),
              reduce_ratio(reduce_ratio > 0 ? reduce_ratio : 0), accuracy_factor(accuracy_factor), strong_bkz(strong_bkz),
              slight_bkz(slight_bkz), min_eqns(min_eqns)
    { }

    FactoringSettings(const ZZ &N, unsigned long n, const RR &c, int max_level,
                      double A_start_factor, double reduce_ratio, long accuracy_factor,
                      long strong_bkz, long slight_bkz, unsigned long min_eqns,
                      long long int seed_type)
            : N(N), n(n), c(c), max_level(max(10, max_level)), A_start_factor(A_start_factor),
              reduce_ratio(reduce_ratio > 0 ? reduce_ratio : 0), accuracy_factor(accuracy_factor), strong_bkz(strong_bkz),
              slight_bkz(slight_bkz), min_eqns(min_eqns), seed_type(seed_type)
    { }

    FactoringSettings(const ZZ &N, unsigned long n, const RR &c, int max_level,
                      double A_start_factor, double reduce_ratio, long accuracy_factor,
                      long strong_bkz, long slight_bkz, unsigned long min_eqns,
                      long long int seed_type, int scalingType)
            : N(N), n(n), c(c), max_level(max(10, max_level)), A_start_factor(A_start_factor),
              reduce_ratio(reduce_ratio > 0 ? reduce_ratio : 0), accuracy_factor(accuracy_factor), strong_bkz(strong_bkz),
              slight_bkz(slight_bkz), min_eqns(min_eqns), seed_type(seed_type), scalingType(scalingType)
    { }
};

#endif //FACTORINGINTEGERSVIALATTICEALGORITHMS_FACTORINGSETTINGS_H