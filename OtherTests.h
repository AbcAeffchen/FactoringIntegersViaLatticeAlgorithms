
#ifndef FACTORINGINTEGERSVIALATTICEALGORITHMS_OTHERTESTS_H
#define FACTORINGINTEGERSVIALATTICEALGORITHMS_OTHERTESTS_H

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/LLL.h>
#include <NTL/vector.h>

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <set>
#include <vector>
#include <time.h>
#include <cstdlib>

#include <unistd.h>
#include <stdlib.h>
#include <cstdio>
using namespace NTL;
using namespace std;

const vector<long> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
                             61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
                             131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
                             197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269,
                             271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
                             353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431,
                             433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
                             509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599,
                             601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673,
                             677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761,
                             769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857,
                             859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947,
                             953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031,
                             1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097,
                             1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187,
                             1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277,
                             1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327,
                             1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439,
                             1447, 1451, 1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499,
                             1511, 1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583,
                             1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663,
                             1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747,
                             1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 1823, 1831, 1847,
                             1861, 1867, 1871, 1873, 1877, 1879, 1889, 1901, 1907, 1913, 1931,
                             1933, 1949, 1951, 1973, 1979, 1987};

void BasisTests(ZZ N, RR c, long dim, long blockSize)
{
    Mat<ZZ> B,U;

    B.SetDims(dim, dim + 1);    // transposed
    // Setting the basis
    for(long i = 1; i <= dim; i++)
    {
        B(i,i) = conv<ZZ>(10000 * SqrRoot(log(conv<RR>(primes[i-1]))));
        if(c < 1)
            B(i, dim + 1) = conv<ZZ>(pow(conv<RR>(N),c) * 10000 * log(conv<RR>(primes[i-1])));
        else
            B(i, dim + 1) = conv<ZZ>(conv<RR>(N * 10000) * log(conv<RR>(primes[i-1])));
    }

    BKZ_FP(B, U, 0.99, blockSize);                 // strong reducing

    transpose(U,U);

    // open file

    system("mkdir output");

    time_t now = time(0);
    struct tm  tstruct;
    tstruct = *localtime(&now);
    char buf[80];
    strftime(buf, sizeof(buf), "./output/%Y-%m-%d_%H-%M-%S_BasisTests.tex", &tstruct);
    fstream basisTest;
    basisTest.open(buf,ios::out);

    basisTest << "\\documentclass[a4paper,twoside,10pt]{report}" << endl << endl
              << "\\usepackage[paperheight=" << 30.0/90.0 * dim + 2 << "cm,paperwidth=" << 62.0/90.0 * dim << "cm, left=1cm, right=1cm, top=2cm,bottom=2cm]{geometry}" << endl
              << "\\usepackage{amsmath}" << endl
              << "\\usepackage{setspace}" << endl
              << "\\pagestyle{empty}" << endl << endl
              << "\\setcounter{MaxMatrixCols}{300}" << endl
              << "\\begin{document}" << endl
              << "\\onehalfspacing" << endl << endl;

    // print transformation
    basisTest << "\\subsection*{Transformation Matrix}" << endl
              << "\\tiny" << endl
              << "\\[ \\begin{bmatrix}" << endl;

    for(long i = 1; i <= dim; i++)
    {
        if(U(i,1) != 0)
            basisTest << U(i,1);
        for(long j = 2; j < dim; j++)
        {
            basisTest << " & ";
            if(U(i,j) != 0)
                basisTest << U(i,j);
        }
        basisTest << " \\\\" << endl;
    }

    basisTest << "\\end{bmatrix} \\]" << endl;

    Mat<RR> mu;
    Vec<RR> r_ii_square;
    ComputeGS(B,mu,r_ii_square);

    RR quot_min, quot_max, quot_avg;

    quot_min = r_ii_square[1]/r_ii_square[0];
    quot_max = r_ii_square[1]/r_ii_square[0];
    quot_avg = r_ii_square[1]/r_ii_square[0];

    for(long i = 2; i < dim; i++)
    {
        if(r_ii_square[i]/r_ii_square[i-1] < quot_min)
            quot_min = r_ii_square[i]/r_ii_square[i-1];
        if(r_ii_square[i]/r_ii_square[i-1] > quot_max)
            quot_max = r_ii_square[i]/r_ii_square[i-1];
        quot_avg += r_ii_square[i]/r_ii_square[i-1];
    }
    quot_avg /= (dim-1);

    basisTest << "$\\min\\frac{r_{i,i}^2}{r_{i-1,i-1}^2}$ " << quot_min << endl
              << "\\\\$\\max\\frac{r_{i,i}^2}{r_{i-1,i-1}^2}$ " << quot_max << endl
              << "\\\\avg $\\frac{r_{i,i}^2}{r_{i-1,i-1}^2}$ " << quot_avg << endl;

    basisTest << "\\end{document}" << endl;

    basisTest.close();
}

#endif //FACTORINGINTEGERSVIALATTICEALGORITHMS_OTHERTESTS_H