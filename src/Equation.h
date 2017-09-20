
#ifndef EQUATION_H
#define	EQUATION_H

#include <NTL/ZZ.h>
#include <NTL/vector.h>

/**
 * Class Equation Data structure to save the equation data
 */
struct Equation
{
    NTL::Vec<long> e;
    NTL::ZZ v;
    long level;
    double reduced;
    long round;
    double time;
    long counter;
    bool fromContinuedFraction;

    Equation(NTL::Vec<long> e, NTL::ZZ v, long level, double reduced, long round, double time, bool fromContinuedFraction = false)
        : e(std::move(e)),
          v(std::move(v)),
          level(level),
          reduced(reduced),
          round(round),
          time(time),
          counter(1),
          fromContinuedFraction(fromContinuedFraction)
    {}

};

/**
 * @param left
 * @param right
 * @return bool
 */
inline bool operator==(const Equation& left, const Equation& right)
{
    return static_cast<bool>(left.e == right.e);    // NTLs vector comparator returns a long with 0 and 1 as only possible value...
}


/**
 * @param left
 * @param right
 * @return bool
 */
inline bool operator<(const Equation& left, const Equation& right)
{
    long n = left.e.length();

    // if left and right have not the same length
    if(n < right.e.length())
        return true;
    else if(n > right.e.length())
        return false;

    // compare the exponents
    for(long i = 0; i < n; ++i)
    {
        if(left.e[i] < right.e[i])
            return true;
        else if(left.e[i] > right.e[i])
            return false;
    }

    // they are equal
    return false;
}

/**
 * Sorts equations by round and time
 * @param left
 * @param right
 * @return
 */
inline bool sort_equations(const Equation& left, const Equation& right)
{
    if(left.round < right.round)
        return true;
    else if(left.round > right.round)
        return false;

    return left.time < right.time;
}

/**
 * Sorts equations by the value of v
 * @param left
 * @param right
 * @return
 */
inline bool sort_equations_by_v(const Equation& left, const Equation& right)
{
    return left.v < right.v != 0;
}

#endif	/* EQUATION_H */