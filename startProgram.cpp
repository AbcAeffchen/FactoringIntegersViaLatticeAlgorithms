
#include "Factoring.h"
#include "OtherTests.h"
#include <NTL/version.h>

using namespace std;

void manual_settings_input()
{
    long N_size;
    cout << endl << "N ~ 10^";
    cin >> N_size;

    unsigned long n;
    cout << endl << "n is the dimension of the lattice (default is n = 90)" << endl << endl;
    cout << "n = ";
    cin >> n;

    RR c;
    cout << endl << "c is a parameter of the lattice (default c = 0.714)" << endl << endl;
    cout << "c = ";
    cin >> c;

    long s;
    cout << endl << "s is the pruning level (default s = 14)" << endl << endl;
    cout << "s = ";
    cin >> s;

    double reduce_ratio;
    cout << endl << "reduce_ratio: bypass not-reducing the minimal distance if it is lower than reduce_ratio*current minimal distance  (0 = never bypass, 1 = always bypass, default = 0.85)" << endl << endl;
    cout << "reduce_ratio = ";
    cin >> reduce_ratio;

    long bkz_strong;
    long bkz_slight;
    cout << endl << "Block sizes of the strong (default 32) and slight (default 20) BKZ-reduce:" << endl << endl;
    cout << "strong BKZ block size: ";
    cin >> bkz_strong;
    cout << "slight BKZ block size: ";
    cin >> bkz_slight;

    double A_factor;
    cout << endl << "A_factor controls the maximum start value of A = A_factor * 0.25 * sum(r_ii^2) (default alpha = 0.2)" << endl << endl;
    cout << "A_factor = ";
    cin >> A_factor;

    unsigned long min_eqns;
    cout << endl << "How many equations should be found before stopping the programm and calculating statistics" << endl << endl;
    cout << "Min. equations = ";
    cin >> min_eqns;

    long long int seed_type;
    cout << endl << "Choose a seed for the random number generator:" << endl
         << "-2 = timestamp" << endl
         << "-1 = random device (sometimes buggy)" << endl
         << "any non negativ number: this number is used as seed" << endl << endl
         << "Seed type: ";
    cin >> seed_type;

    cout << endl << endl;

    Factoring(FactoringSettings(getN(N_size), n, c, s, A_factor, reduce_ratio, 10000, bkz_strong, bkz_slight, min_eqns, seed_type));
}

void test_series_n()
{
    Factoring(FactoringSettings(getN(20), 120, conv<RR>(0.6), 17, 0.2, 0.82, 10000, 32, 20, 25));

    Factoring(FactoringSettings(getN(20), 130, conv<RR>(0.6), 17, 0.2, 0.82, 10000, 32, 20, 25));

    Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.6), 17, 0.2, 0.82, 10000, 32, 20, 25));

    Factoring(FactoringSettings(getN(20), 150, conv<RR>(0.6), 17, 0.2, 0.82, 10000, 32, 20, 25));

    Factoring(FactoringSettings(getN(20), 160, conv<RR>(0.6), 17, 0.2, 0.82, 10000, 32, 20, 25));
}

void test_series_c()
{
    Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.4), 17, 0.2, 0.82, 10000, 32, 20, 25));

    Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.5), 17, 0.2, 0.82, 10000, 32, 20, 25));

    Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.6), 17, 0.2, 0.82, 10000, 32, 20, 25));

    Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.7), 17, 0.2, 0.82, 10000, 32, 20, 25));

    Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.8), 17, 0.2, 0.82, 10000, 32, 20, 25));
}

void test_series_prune()
{
    Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.6), 15, 0.2, 0.82, 10000, 32, 20, 25));

    Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.6), 16, 0.2, 0.82, 10000, 32, 20, 25));

    Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.6), 17, 0.2,0.82, 10000, 32, 20, 25));

    Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.6), 18, 0.2, 0.82, 10000, 32, 20, 25));

    Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.6), 19, 0.2, 0.82, 10000, 32, 20, 25));
}

void default_settings()
{
    Factoring(FactoringSettings(getN(14), 90, conv<RR>(5.0/7.0), 14, 0.2, 0.8, 10000, 32, 20, 91));
}

void speedTestSmall()
{
    Factoring(FactoringSettings(getN(14), 90, conv<RR>(5.0/7.0), 14, 0.2, 0.8, 10000, 32, 20, 91,1446721645,true,0));
}

void speedTestBig()
{
    Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 50,-2,true,0));
}

void fast_test()
{
    Factoring(FactoringSettings(getN(14), 90, conv<RR>(0.5), 14, 0.2, 0.8, 10000, 32, 20, 40,34992116121LL));
}

void default_settings_big_N()
{
    Factoring(FactoringSettings(getN(20), 150, conv<RR>(1.0/2.0), 19, 0.2, 0.75, 10000, 32, 20, 20));
}

void test_series_1()
{
    Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 50,-2,true,0));
    Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 50,-2,true,0));

    Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 50,-2,true,4));
    Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 50,-2,true,5));
    Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 50,-2,true,4));
    Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 50,-2,true,5));

    Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 50,-2,true,3));
}

void menu()
{
    short choice;
    cout << "Choose what to do:" << endl
         << "(1) Run with default settings " << endl
         << "    -> N = 10^14, n = 90, c = 5/7, prune = 14, min. equations = 91" << endl
         << "(2) Run n test series" << endl
         << "    -> N = 10^20, n = 120 to 170 by steps of 10," << endl
         << "       c = 0.6, prune = 17, min. equations = 50" << endl
         << "(3) Run c test series" << endl
         << "    -> N = 10^20, n = 140, c = 0.4 to 0.8 by steps of 0.1," << endl
         << "       prune = 17, min. equations = 50" << endl
         << "(4) Run prune test series" << endl
         << "    -> N = 10^20, n = 140, c = 0.6, prune = 15 to 19 by steps of 1," << endl
         << "       min. equations = 50" << endl << endl
         << "(9) Run with custom settings" << endl << endl
         << "(14) Run the speedTest (N~10^14)" << endl << endl
         << "(20) Run the speedTest (N~10^20)" << endl << endl
         << "(100) Prime Lattice Basis tests" << endl << endl
         << "(200) Test Series 1" << endl << endl
         << "Quit by [Ctrl]+[C]" << endl
         << "Selection: ";
        cin >> choice;

    switch(choice)
    {
    case 1:
        default_settings();
        break;
    case 2:
        test_series_n();
        break;
    case 3:
        test_series_c();
        break;
    case 4:
        test_series_prune();
        break;
    case 9:
        manual_settings_input();
        break;
    case 14:
        speedTestSmall();
        break;
    case 20:
        speedTestBig();
        break;
    case 100:
    {
        long e, dim, blockSize;
        RR c;
        cout << "N ~ 10^";
        cin >> e;
        cout << "Dim (<= 300): ";
        cin >> dim;
        cout << "c = ";
        cin >> c;
        cout << "Block size: ";
        cin >> blockSize;
        BasisTests(getN(e), c, dim, blockSize);
        break;
    }
    case 200:
        test_series_1();
    default:
        menu();
    }
}

int main()
{
    RR::SetPrecision(150);
    RR::SetOutputPrecision(20);

    cout << "GCC: " << __GNUC__ << "." <<__GNUC_MINOR__ << "."  << __GNUC_PATCHLEVEL__ << endl;
    cout << "NTL: " << NTL_VERSION << endl << endl;

    menu();

    return 0;
}