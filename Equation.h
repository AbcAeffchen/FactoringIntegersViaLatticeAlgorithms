
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

    Equation(NTL::Vec<long> e, NTL::ZZ v, long level, double reduced, long round, double time);
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

bool sort_equations(Equation left, Equation right);

#endif	/* EQUATION_H */

