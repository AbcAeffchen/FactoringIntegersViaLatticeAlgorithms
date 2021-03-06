
#include "Factoring.h"
#include "OtherTests.h"
#include <NTL/version.h>

using namespace std;

void manual_settings_input()
{
    long N_size;
    cout << "\nN ~ 10^";
    cin >> N_size;

    unsigned long n;
    cout << "\nn is the dimension of the lattice (default is n = 90)\n\n";
    cout << "n = ";
    cin >> n;

    RR c;
    cout << "\nc is a parameter of the lattice (default c = 0.714)\n\n";
    cout << "c = ";
    cin >> c;

    long s;
    cout << "\ns is the pruning level (default s = 14)\n\n";
    cout << "s = ";
    cin >> s;

    double reduce_ratio;
    cout << "\nreduce_ratio: bypass not-reducing the minimal distance if it is lower than reduce_ratio*current minimal distance  (0 = never bypass, 1 = always bypass, default = 0.85)\n\n";
    cout << "reduce_ratio = ";
    cin >> reduce_ratio;

    long bkz_strong;
    long bkz_slight;
    cout << "\nBlock sizes of the strong (default 32) and slight (default 20) BKZ-reduce:\n\n";
    cout << "strong BKZ block size: ";
    cin >> bkz_strong;
    cout << "slight BKZ block size: ";
    cin >> bkz_slight;

    double A_factor;
    cout << "\nA_factor controls the maximum start value of A = A_factor * 0.25 * sum(r_ii^2) (default alpha = 0.2)\n\n";
    cout << "A_factor = ";
    cin >> A_factor;

    unsigned long min_eqns;
    cout << "\nHow many equations should be found before stopping the program and calculating statistics\n\n"
            "Min. equations = ";
    cin >> min_eqns;

    long long int seed_type;
    cout << "\nChoose a seed for the random number generator:\n"
            "-2 = timestamp (default)\n"
            "-1 = random device (sometimes buggy)\n"
            "any non negative number: this number is used as seed\n\n"
            "Seed type: ";
    cin >> seed_type;

    char cf;
    cout << "\nUse continued fractions (y/n) (defaults to yes): ";
    cin >> cf;

    int scalingType;
    cout << "\nChoose the scaling type:\n"
            "0 = mixed (default)\n"
            "1 = every row with probability 1/2\n"
            "2 = every row with probability 1/4\n"
            "3 = every row with probability 3/4\n"
            "4 = first n/2 rows with prop. 1/4, other with prop. 1/2\n"
            "5 = first n/2 rows with prop. 1/2, other with prop. 1/4\n"
            "Scaling type: ";
    cin >> scalingType;

    cout << "\n\n";

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

void randomSpeedTestsSmall()
{
    unsigned int k;
    cout << "Number of runs to make: ";
    cin >> k;

    for(unsigned int i = 0; i < k; i++)
        Factoring(FactoringSettings(getN(14), 90, conv<RR>(5.0/7.0), 14, 0.2, 0.8, 10000, 32, 20, 91,-2,true,0));
}

void randomSpeedTestsBig()
{
    unsigned int k;
    cout << "Number of runs to make: ";
    cin >> k;

    for(unsigned int i = 0; i < k; i++)
        Factoring(FactoringSettings(getN(20), 150, conv<RR>(1/2.0), 17, 0.2, 0.8, 10000, 32, 20, 25,-2,true,0));
}

void menu()
{
    short choice;
    cout << "Choose what to do:\n"
            "(  0) Run with custom settings\n"
            "(  1) Run the speedTest (N~10^14)\n"
            "(  2) Run the speedTest (N~10^20)\n\n"
            "Performance Tests\n"
            "(  5) Run multiple speedTests with random seeds (N~10^14)\n"
            "(  6) Run multiple speedTests with random seeds (N~10^20)\n\n"
            "Test series\n"
            "( 10) n test series\n"
            "( 11) c test series\n"
            "( 12) Prune test series\n"
            "( 13) Scaling test series\n"
            "( 14) Start Factor test series\n\n"
            "(100) Custom Prime Lattice Basis Test\n"
            "(101) Prime Lattice Basis Testseries\n\n"
            "Quit by [Ctrl]+[C]\n"
            "Selection: ";
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
    case 5:
        randomSpeedTestsSmall();
        break;
    case 6:
        randomSpeedTestsBig();
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
        BasisTests(getN(e),c,dim,blockSize);
        break;
    }
    case 101:
    {
        BasisTests(getN(14),conv<RR>(5.0/7.0),100,20);
        BasisTests(getN(14),conv<RR>(5.0/7.0),200,20);
        BasisTests(getN(14),conv<RR>(5.0/7.0),300,20);
        BasisTests(getN(14),conv<RR>(5.0/7.0),100,32);
        BasisTests(getN(14),conv<RR>(5.0/7.0),200,32);
        BasisTests(getN(14),conv<RR>(5.0/7.0),300,32);
        BasisTests(getN(14),conv<RR>(5.0/7.0),100,50);
        BasisTests(getN(14),conv<RR>(5.0/7.0),200,50);
        BasisTests(getN(14),conv<RR>(5.0/7.0),300,50);
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

    cout << "GCC: " << __GNUC__ << "." <<__GNUC_MINOR__ << "."  << __GNUC_PATCHLEVEL__
         << "\nNTL: " << NTL_VERSION << "\n" << endl;

    menu();

    return 0;
}