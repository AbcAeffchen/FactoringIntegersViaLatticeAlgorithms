
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
    cout << endl << "How many equations should be found before stopping the program and calculating statistics" << endl << endl;
    cout << "Min. equations = ";
    cin >> min_eqns;

    long long int seed_type;
    cout << endl << "Choose a seed for the random number generator:" << endl
         << "-2 = timestamp (default)" << endl
         << "-1 = random device (sometimes buggy)" << endl
         << "any non negative number: this number is used as seed" << endl << endl
         << "Seed type: ";
    cin >> seed_type;

    char cf;
    cout << endl << "Use continued fractions (y/n) (defaults to yes): ";
    cin >> cf;

    bool scalingType;
    cout << endl << "Choose the scaling type:" << endl
         << "0 = mixed (default)" << endl
         << "1 = every row with probability 1/2" << endl
         << "2 = every row with probability 1/4" << endl
         << "3 = every row with probability 3/4" << endl
         << "4 = first n/2 rows with prop. 1/4, other with prop. 1/2" << endl
         << "5 = first n/2 rows with prop. 1/2, other with prop. 1/4" << endl
         << "Scaling type: ";
    cin >> scalingType;

    cout << endl << endl;

    Factoring(FactoringSettings(getN(N_size), n, c, s, A_factor, reduce_ratio, 10000, bkz_strong, bkz_slight, min_eqns, seed_type,cf =='y',scalingType));
}

void test_series_n()
{
//    for(long i = 0; i < 5; i++)
//    {
//        for(long dim = 170; dim <= 240; dim += 10)
//        {
//            Factoring(FactoringSettings(getN(20), dim, conv<RR>(0.5), 17, 0.2, 0.82, 10000, 32, 20, 25));
//        }
//    }
}

void test_series_c()
{
    // done
}

void test_series_prune()
{
//    for(long i = 0; i < 2; i++)
//    {
//        for(long level = 15; level <= 19; level++)
//            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.6), level, 0.2, 0.82, 10000, 32, 20, 25));
//    }
}

void test_series_start_factor()
{
    // all is done
}

void test_series_scaling()
{
    for (long i = 0; i < 6; i++)
    {
        Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 25,-2,true,FS_SCALING_MIXED));
        Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 25,-2,true,FS_SCALING_ONE_HALF));
        Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 25,-2,true,FS_SCALING_ONE_FOURTH));
        Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 25,-2,true,FS_SCALING_THREE_FORTH));
        Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 25,-2,true,FS_SCALING_ONE_FOURTH_ONE_HALF));
        Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 25,-2,true,FS_SCALING_ONE_HALF_ONE_FOURTH));
    }
}

void speedTestSmall()
{
    Factoring(FactoringSettings(getN(14), 90, conv<RR>(5.0/7.0), 14, 0.2, 0.8, 10000, 32, 20, 91,1448726806,true,0));
}

void speedTestBig()
{
    Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 25,1448682782,true,0));
}

void menu()
{
    short choice;
    cout << "Choose what to do:" << endl
         << "(  0) Run with custom settings" << endl
         << "(  1) Run the speedTest (N~10^14)" << endl
         << "(  2) Run the speedTest (N~10^20)" << endl
         << endl
         << "Test series" << endl
         << "( 10) n test series" << endl
         << "( 11) c test series" << endl
         << "( 12) Prune test series" << endl
         << "( 13) Scaling test series" << endl
         << "( 14) Start Factor test series" << endl
         << endl
         << "(100) Prime Lattice Basis tests" << endl << endl
         << "Quit by [Ctrl]+[C]" << endl
         << "Selection: ";
        cin >> choice;

    switch(choice)
    {
    case 0:
        manual_settings_input();
        break;
    case 1:
        speedTestSmall();
        break;
    case 2:
        speedTestBig();
        break;
    case 10:
        test_series_n();
        break;
    case 11:
        test_series_c();
        break;
    case 12:
        test_series_prune();
        break;
    case 13:
        test_series_scaling();
        break;
    case 14:
        test_series_start_factor();
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
        break;
    }
    case 1000:
    {
        // cf on
        for(int i = 0; i < 50; i++)
            Factoring(FactoringSettings(getN(14), 90, conv<RR>(5.0/7.0), 14, 0.2, 0.8, 10000, 32, 20, 92, -2, true, FS_SCALING_MIXED));

        break;
    }
    case 2000:
    {
        Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.8), 17, 0.2, 0.82, 10000, 32, 20, 25,-2,true,0));
        // c (5*4 = 20)
        for(int i = 0; i < 2; i++)
        {
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.4), 17, 0.2, 0.82, 10000, 32, 20, 25,-2,true,0));
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.5), 17, 0.2, 0.82, 10000, 32, 20, 25,-2,true,0));
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.6), 17, 0.2, 0.82, 10000, 32, 20, 25,-2,true,0));
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.7), 17, 0.2, 0.82, 10000, 32, 20, 25,-2,true,0));
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.8), 17, 0.2, 0.82, 10000, 32, 20, 25,-2,true,0));
        }

        // n (13*4 = 52)
        for(int i = 0; i < 4; i++)
        {
            for(long dim = 120; dim <= 240; dim += 10)
                Factoring(FactoringSettings(getN(20), dim, conv<RR>(0.5), 17, 0.2, 0.82, 10000, 32, 20, 25,-2,true,0));
        }

        // prune (5*4 = 20)
        for(int i = 0; i < 4; i++)
        {
            for(long prune = 15; prune <= 19; prune++)
                Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.5), prune, 0.2, 0.82, 10000, 32, 20, 25, -2, true, 0));
        }

        // scaling (6*4 = 24)
        for(int i = 0; i < 4; i++)
        {
            Factoring(FactoringSettings(getN(20), 150, conv<RR>(0.5), 17, 0.2, 0.8, 10000, 32, 20, 25, -2, true, FS_SCALING_MIXED));
            Factoring(FactoringSettings(getN(20), 150, conv<RR>(0.5), 17, 0.2, 0.8, 10000, 32, 20, 25, -2, true, FS_SCALING_ONE_HALF));
            Factoring(FactoringSettings(getN(20), 150, conv<RR>(0.5), 17, 0.2, 0.8, 10000, 32, 20, 25, -2, true, FS_SCALING_ONE_FOURTH));
            Factoring(FactoringSettings(getN(20), 150, conv<RR>(0.5), 17, 0.2, 0.8, 10000, 32, 20, 25, -2, true, FS_SCALING_THREE_FORTH));
            Factoring(FactoringSettings(getN(20), 150, conv<RR>(0.5), 17, 0.2, 0.8, 10000, 32, 20, 25, -2, true, FS_SCALING_ONE_FOURTH_ONE_HALF));
            Factoring(FactoringSettings(getN(20), 150, conv<RR>(0.5), 17, 0.2, 0.8, 10000, 32, 20, 25, -2, true, FS_SCALING_ONE_HALF_ONE_FOURTH));
        }

        // min reduction (8*4=32)
        for(int i = 0; i < 4; i++)
        {
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.5), 17, 0.18, 0.82, 10000, 32, 20, 25, -2, true, FS_SCALING_MIXED));
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.5), 17, 0.19, 0.82, 10000, 32, 20, 25, -2, true, FS_SCALING_MIXED));
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.5), 17, 0.20, 0.82, 10000, 32, 20, 25, -2, true, FS_SCALING_MIXED));
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.5), 17, 0.21, 0.82, 10000, 32, 20, 25, -2, true, FS_SCALING_MIXED));
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.5), 17, 0.22, 0.82, 10000, 32, 20, 25, -2, true, FS_SCALING_MIXED));
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.5), 17, 0.23, 0.82, 10000, 32, 20, 25, -2, true, FS_SCALING_MIXED));
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.5), 17, 0.24, 0.82, 10000, 32, 20, 25, -2, true, FS_SCALING_MIXED));
            Factoring(FactoringSettings(getN(20), 140, conv<RR>(0.5), 17, 0.25, 0.82, 10000, 32, 20, 25, -2, true, FS_SCALING_MIXED));
        }



        break;
    }
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