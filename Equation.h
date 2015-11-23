
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
    long threadId;
    double time;
    long counter;
    bool fromContinuedFraction;

    Equation(NTL::Vec<long> e, NTL::ZZ v, long level, double reduced, long round, long threadId, double time, bool fromContinuedFraction = false);
};

/**
 * @param left
 * @param right
 * @return bool
 */
bool operator== (const Equation &left, const Equation &right);

/**
 * @param left
 * @param right
 * @return bool
 */
bool operator< (const Equation &left, const Equation &right);

bool sort_equations(const Equation &left, const Equation &right);

bool sort_equations_by_v(const Equation &left, const Equation &right);

#endif	/* EQUATION_H */