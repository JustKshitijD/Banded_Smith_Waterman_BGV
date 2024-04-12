/* Copyright (C) 2019 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */
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

// void compare_function(char *argv[]);

// some parameters for quick testing
// B 7 1 75 90 1 10 y
// B 7 1 300 90 1 10 y
// U 17 1 145 120 1 10 y
Comparator compare_function(char *argv[]) {

    CircuitType type = UNI;
    if (!strcmp(argv[7], "B"))
    {
        type = BI;
    }
    else if (!strcmp(argv[7], "T"))
    {
        type = TAN;
    }
    else if (!strcmp(argv[7], "U"))
    {
        type = UNI;
    }
    else
    {
        throw invalid_argument("Choose a valid circuit type (U for univariate, B for bivariate and T for Tan et al.\n");
    }

    //////////PARAMETER SET UP////////////////
    // Plaintext prime modulus
    unsigned long p = atol(argv[0]);
    // Field extension degree
    unsigned long d = atol(argv[1]);
    // Cyclotomic polynomial - defines phi(m)
    unsigned long m = atol(argv[2]);
    // Number of ciphertext prime bits in the modulus chain
    unsigned long nb_primes = atol(argv[3]);
    // Number of columns of Key-Switching matix (default = 2 or 3)
    unsigned long c = 3;
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

    // Print the security level
    cout << "Q size: " << context.logOfProduct(context.getCtxtPrimes())/log(2.0) << endl;
    cout << "Q*P size: " << context.logOfProduct(context.fullPrimes())/log(2.0) << endl;
    cout << "Security: " << context.securityLevel() << endl;

    
  // Get the EncryptedArray of the context
  const helib::EncryptedArray& ea = context.getEA();

  // Get the number of slot (phi(m))
  long nslots = ea.size();
  std::cout << "Number of slots: " << nslots << std::endl;

    // Print the context
    context.getZMStar().printout();
    cout << endl;

    //maximal number of digits in a number
    unsigned long expansion_len = atol(argv[4]);

    // Secret key management
    cout << "Creating secret key..." << endl;
    // Create a secret key associated with the context
    SecKey secret_key(context);
    // Generate the secret key
    secret_key.GenSecKey();
    cout << "Generating key-switching matrices..." << endl;
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

    // Simple case- expansion_len=number of digits in number=1; d=1- Both numbers to be compared are <p
    Comparator comparator(context, type, d, expansion_len, secret_key, true);

    //repeat experiments several times
    int runs = atoi(argv[6]);

    // //test comparison circuit
    // comparator.test_compare(runs);

    // cout<<"Done with test_compare"<<endl;

    return comparator;

    //printAllTimers(cout);
}

