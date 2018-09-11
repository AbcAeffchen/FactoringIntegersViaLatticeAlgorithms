//
// Created by Alex Schickedanz <alex@ae.cs.uni-frankfurt.de> on 11.09.18.
//

#ifndef FACTORINGINTEGERSVIALATTICEALGORITHMS_NEWENUMSTAGE_H
#define FACTORINGINTEGERSVIALATTICEALGORITHMS_NEWENUMSTAGE_H

#include <NTL/LLL.h>
#include <queue>
#include "Equation.h"
#include "Statistics.h"

/**
 * This is a data structure to save a stage of NewEnum
 */
struct NewEnumStage
{
    RR y_t;                     // equals y(t)
    RR c_tp1;                   // equals c(t+1)
    Vec<double> u;              // equals Vec<RR> u but need less memory. u contains only integers, so it's not important if they are stored in RR or double
    long t;                     // current coordinate

    NewEnumStage()
        : y_t(conv<RR>(0)), c_tp1(conv<RR>(0)), t(0)
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

#endif //FACTORINGINTEGERSVIALATTICEALGORITHMS_NEWENUMSTAGE_H
