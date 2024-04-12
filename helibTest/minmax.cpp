#include <iostream>
#include <time.h>
#include <random>

#include <helib/helib.h>
#include <helib/debugging.h>
#include <helib/Context.h>
#include <helib/polyEval.h>
#include "tools.h"
#include "comparator.h"

using namespace std;
using namespace NTL;
using namespace helib;
using namespace he_cmp;

void minmax_function(char *variableList[]);

// some parameters for quick testing
// 7 1 75 90 1 4 1 10 y
// 7 1 300 90 1 6 2 10 y
// 17 1 145 120 1 7 2 10 y
void minmax_function(char *variableList[]) {
    unsigned long p = atol(variableList[0]);
    // Field extension degree
    unsigned long d = atol(variableList[1]);
    // Cyclotomic polynomial - defines phi(m)
    unsigned long m = atol(variableList[2]);
    // Number of ciphertext prime bits in the modulus chain
    unsigned long nb_primes = atol(variableList[3]);
    // Number of columns of Key-Switching matrix (default = 2 or 3)
    unsigned long c = 2;
    cout << "Initialising context object..." << endl;
    // Intialise context
    auto context = ContextBuilder<BGV>()
            .m(m)
            .p(p)
            .r(1)
            .bits(nb_primes)
            .c(c)
            .scale(6)
            .build();

    context.printout();
    cout << endl;

    std::cout << "Security: " << context.securityLevel() << std::endl;

    //maximal number of digits in a number
    unsigned long expansion_len = atol(variableList[4]);

    std::cout << "Creating secret key..." << std::endl;
    // Create a secret key associated with the context
    helib::SecKey secret_key(context);
    // Generate the secret key
    secret_key.GenSecKey();
    std::cout << "Generating key-switching matrices..." << std::endl;
    // Compute key-switching matrices that we need
    if (expansion_len > 1)
    {
        if (context.getZMStar().numOfGens() == 1)
        {
            std::set<long> automVals;
            long e = 1;
            long ord = context.getZMStar().OrderOf(0);
            bool native = context.getZMStar().SameOrd(0);
            if(!native)
                automVals.insert(context.getZMStar().genToPow(0, -ord));
            while (e < expansion_len){
                long atm = context.getZMStar().genToPow(0, ord-e);
                //cout << "Automorphism " << -e << " is " << atm << endl;
                automVals.insert(atm);
                e <<=1;
            }
            addTheseMatrices(secret_key, automVals);
        }
        else
        {
            addSome1DMatrices(secret_key);
        }
    }

    if (d > 1)
        addFrbMatrices(secret_key); //might be useful only when d > 1

    // create Comparator (initialize after buildModChain)
    Comparator comparator(context, UNI, d, expansion_len, secret_key, true);

    // number of values to be sorted
    int input_len = atoi(variableList[5]);

    // levels of consecutive comparisons
    long depth = atol(variableList[8]);

    //repeat experiments 'runs' times
    int runs = atoi(variableList[6]);

    //test sorting
    comparator.test_array_min(input_len, depth, runs);

    printAllTimers(cout);

}
