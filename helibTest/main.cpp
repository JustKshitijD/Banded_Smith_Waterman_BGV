#include <iostream>
#include <cstdlib> 
#include <vector>
#include <sys/time.h>
#include <time.h>
#include <getopt.h>
#include <string>
#include <filesystem>
#include "omp.h"
#include <random>

#include <helib/helib.h>
#include <helib/debugging.h>
#include <helib/Context.h>
#include <helib/polyEval.h>
#include "tools.h"
#include "comparator.h"
#include "sorting.cpp"
#include "compare.cpp"
#include "minmax.cpp"
#include "host_data_io.h"
#include "host_data.h"
#include "host_kernel.h"
#include "comparator.h"
#include <iostream>
#include <fstream>
#include <helib/Context.h>

#include <map> 
#include <NTL/ZZ_pE.h>
#include <NTL/mat_ZZ_pE.h>
#include <helib/Ptxt.h>
#include <fstream>

#include <chrono>
using namespace std::chrono;

using namespace std;
using namespace NTL;
using namespace helib;
using namespace he_cmp;

// the main function that takes  arguments (type in Terminal: ./sorting_circuit argv[1] argv[2] argv[3] argv[4] argv[5] argv[6] argv[7] argv[8])
// variable[0] - the plaintext modulus
// variable[1] - the dimension of a vector space over a finite field
// variable[2] - the order of the cyclotomic ring
// variable[3] - the bitsize of the ciphertext modulus in ciphertexts (HElib increases it to fit the moduli chain). The modulus used for public-key generation
// variable[4] - the length of vectors to be compared
// variable[5] - the number of values to be sorted
// variable[6] - the number of experiment repetitions
// variable[7] - circuit type (U, B or T)
// variable[8] - the number of tournament stages

Ptxt<helib::BGV> create_plaintext_from_constant(long k, int count, Comparator comparator)
{
    helib::Ptxt<helib::BGV> head(comparator.m_context);
    for(int j=0;j<count;j++)
    {
        head[j]=k;
    }
    return head;
}

long compare_cnt=0;

Ctxt compare_ctxt(Ctxt c, Ctxt c2, Comparator comparator,const helib::PubKey& m_pk, SecKey m_sk, unsigned long ord_p)
{
    Ctxt ctxt_res(m_pk);
    auto start = std::chrono::system_clock::now();
    comparator.compare(ctxt_res, c, c2);
    auto end = std::chrono::system_clock::now();
    if(compare_cnt<5)
    {
     auto elapsed = end - start;
     std::cout << "Time for compare: "<<elapsed.count() << '\n';
    }
    
    compare_cnt++;
    return ctxt_res;
}


// Possible optimization - Implement comparison between Ciphertext and Plaintext instead of CT-CT comparison
// returns ctxt = (c<k)
Ctxt compare_ctxt_num(long int k, Ctxt c, Comparator comparator,const helib::EncryptedArray& ea,long nslots,const helib::PubKey& m_pk,SecKey m_sk, unsigned long ord_p)
{
    vector<ZZX> pol_x(nslots);

	for(int kk=0;kk<nslots;kk++)
	{
		pol_x[kk]=ZZX(INIT_MONO, 0, k);
	}

    long nSlots = ea.size();
    
    Ctxt ctxt_x(m_pk);
    ea.encrypt(ctxt_x, m_pk, pol_x);
    
    Ctxt ctxt_res(m_pk);

	auto start = high_resolution_clock::now();
    comparator.compare(ctxt_res,ctxt_x,c);
	auto stop = high_resolution_clock::now();
   
	return ctxt_res;

}

// returns ctxt = (c<k)
Ctxt compare_ctxt_num(Ctxt c, long int k, Comparator comparator,const helib::EncryptedArray& ea,long nslots,const helib::PubKey& m_pk,SecKey m_sk, unsigned long ord_p)
{
    vector<ZZX> pol_x(nslots);

	for(int kk=0;kk<nslots;kk++)
	{
		pol_x[kk]=ZZX(INIT_MONO, 0, k);
	}

    long nSlots = ea.size();
    
    Ctxt ctxt_x(m_pk);
    ea.encrypt(ctxt_x, m_pk, pol_x);
    
    Ctxt ctxt_res(m_pk);

	auto start = high_resolution_clock::now();
    comparator.compare(ctxt_res, c, ctxt_x);
	auto stop = high_resolution_clock::now();
   
	return ctxt_res;

}

Ptxt<helib::BGV> compare_ptxt(Ptxt<helib::BGV> a, Ptxt<helib::BGV> b, Comparator comparator, unsigned long p)
{
    Ptxt<helib::BGV> cmp(comparator.m_context);  Ptxt<helib::BGV> diff(comparator.m_context);
    cmp=a; cmp-=b; diff=cmp; 

    cmp.power(p-1); cmp*=(p+1)/2;

    for(long i=1;i<=p-2;i+=2)
    {
        long c=0;
        for(long a=1;a<=(p-1)/2;a++)
        {
            long k=1;
            for(long j=0;j<p-1-i;j++)
            {
                k*=a; k=k%p;
            }
    
            c+=k; c=c%p;
        }

        Ptxt<helib::BGV> temp(comparator.m_context); temp=diff;
        temp.power(i); 
        temp*=c;
        cmp+=temp;
    }

    return cmp;
}

Ctxt equal_ctxt(Ctxt a, Ctxt b, long p, Comparator comparator,const helib::PubKey& m_pk, SecKey m_sk)
{
    Ctxt res(m_pk); res=a; res-=b;
    res.power(p-1);
    res*=-1l; res+=1l;
    return res;
}

// returns 1 if a==b
Ctxt equal_ctxt(Ctxt a, Ptxt<helib::BGV> b, long p, Comparator comparator,const helib::PubKey& m_pk, SecKey m_sk)
{
    Ctxt res(m_pk); res=a; res-=b;
    res.power(p-1);res*=-1l; res+=1l;
    return res;
}

Ptxt<helib::BGV> equal_ptxt( Ptxt<helib::BGV> a, Ptxt<helib::BGV> b, long p, Comparator comparator,const helib::PubKey& m_pk, SecKey m_sk)
{
    Ptxt<helib::BGV> res(m_pk); res=a; res-=b;
    res.power(p-1); res*=-1l; res+=1l;
    return res;
}

Ctxt assign_comparison_result(Ctxt a, Ptxt<helib::BGV> b, Ctxt cmp_result,const helib::PubKey& m_pk )
{
    Ctxt res(m_pk);
    res=cmp_result; res*=b;
    Ctxt temp(m_pk); temp=cmp_result; temp*=-1l; temp+=1l; temp*=a;
    res+=temp;
    return res;
}

Ctxt assign_comparison_result(Ctxt a, Ptxt<helib::BGV> b, Ptxt<helib::BGV> cmp_result,const helib::PubKey& m_pk )
{
    Ctxt res(m_pk);
    b*=cmp_result; res+=b;
    Ctxt temp(m_pk); temp=a; cmp_result*=-1l; cmp_result+=1l; temp*=cmp_result;
    res+=temp;
    return res;
}

Ctxt assign_comparison_result(Ctxt a, Ctxt b, Ptxt<helib::BGV> cmp_result,const helib::PubKey& m_pk )
{
    Ctxt res(m_pk);
    res=b; res*=cmp_result;
    cmp_result*=-1l; cmp_result+=1l; a*=cmp_result;// *=a;
    res+=a;
    return res;
}

// if cmp_result=0, return a; else, return b
Ctxt assign_comparison_result(Ctxt a, Ctxt b, Ctxt cmp_result,const helib::PubKey& m_pk)
{
    Ctxt res(m_pk);
    res=b; 
    // cout<<"x"<<endl;
    // secret_key.Decrypt(new_plaintext_result, res);
    // cout << "res: " << new_plaintext_result << endl;
    // secret_key.Decrypt(new_plaintext_result, cmp_result);
    // cout << "cmp_result: " << new_plaintext_result << endl;
    res*=cmp_result;
    // cout<<"y"<<endl;
    Ctxt temp(m_pk); temp=cmp_result; 
    // cout<<"z"<<endl;
    temp*=-1l; 
    // cout<<"z1"<<endl;
    temp+=1l; 
    // cout<<"z2"<<endl;
    temp*=a;
    // cout<<"z3"<<endl;
    res+=temp;
    // cout<<"z4"<<endl;
    return res;
}

// if cmp_result=0, return a; else, return b
Ctxt assign_comparison_result(Ctxt a, Ctxt b, Ctxt cmp_result,const helib::PubKey& m_pk,SecKey secret_key,Comparator comparator)
{
    Ptxt<helib::BGV> new_plaintext_result(comparator.m_context);
    cout<<"new_plaintext_result: "<<new_plaintext_result<<endl;

    Ctxt res(m_pk);
    res=b; 
    cout<<"x"<<endl;
    secret_key.Decrypt(new_plaintext_result, res);
    cout << "res: " << new_plaintext_result << endl;
    secret_key.Decrypt(new_plaintext_result, cmp_result);
    cout << "cmp_result: " << new_plaintext_result << endl;
    res*=cmp_result;  // error; a multiplication is giving division by zero error
    cout<<"y"<<endl;
    Ctxt temp(m_pk); temp=cmp_result; 
    cout<<"z"<<endl;
    temp*=-1l; 
    cout<<"z1"<<endl;
    temp+=1l; 
    cout<<"z2"<<endl;
    temp*=a;
    cout<<"z3"<<endl;
    res+=temp;
    cout<<"z4"<<endl;
    return res;
}

// Ctxt get_match_score(Ctxt a, Ctxt b,Ctxt w_match_ctxt,Ctxt w_mismatch_ctxt,Ctxt w_ambig_ctxt, long p, Comparator comparator,const helib::EncryptedArray& ea,long nslots,const helib::PubKey& m_pk,SecKey secret_key, unsigned long ord_p, int count)
// {
//     Ptxt<helib::BGV> new_plaintext_result(comparator.m_context);

//     // secret_key.Decrypt(new_plaintext_result, a);
//     // cout << "a: " << new_plaintext_result << endl;
//     // secret_key.Decrypt(new_plaintext_result, b);
//     // cout << "b: " << new_plaintext_result << endl;

//     Ctxt res(m_pk);
//     new_plaintext_result=create_plaintext_from_constant(4,count,comparator);
//     Ctxt c_ambig_1=equal_ctxt(a,new_plaintext_result,p,comparator,m_pk,secret_key);
//     Ctxt c_ambig_2=equal_ctxt(b, new_plaintext_result, p,comparator,m_pk,secret_key);
//     Ctxt c_ambig_3(m_pk); c_ambig_3=c_ambig_2; c_ambig_3*=-1l; c_ambig_3+=1l;
//     c_ambig_1*=c_ambig_3; c_ambig_1+=c_ambig_2;        
//     c_ambig_1*=-1l; c_ambig_1+=1l;                   // a==4 OR b==4
//     // secret_key.Decrypt(new_plaintext_result, c_ambig_1);
//     // cout << "c_ambig_1: " << new_plaintext_result << endl;

//     Ctxt c_match=equal_ctxt(a, b,p,comparator,m_pk,secret_key);
//     // secret_key.Decrypt(new_plaintext_result, c_match);
//     // cout << "c_match: " << new_plaintext_result << endl;

//     res = assign_comparison_result(w_match_ctxt,w_mismatch_ctxt,c_match,m_pk);
//     // secret_key.Decrypt(new_plaintext_result, res);
//     // cout << "res1: " << new_plaintext_result << endl;

//     res = assign_comparison_result(res,w_ambig_ctxt,c_ambig_1,m_pk);

//     secret_key.Decrypt(new_plaintext_result, res);
//     // cout << "res: " << new_plaintext_result << endl;
//     return res;
// }

// Function to split a string into individual digits
vector<int> splitDigits(const string& str) {
    vector<int> digits;
    for (char digit : str) {
        if (isdigit(digit)) {
            digits.push_back(digit - '0');
        }
    }
    return digits;
}

int main(int argc, char* argv[])
{

    char *variableList[9];

    for(int i=0; i<sizeof(variableList); i++){
        variableList[i] = new char[1000];
    }

    cout << "Please set arguments to all 8 variables:" << endl;
    cout << "plaintext modulus (prime number):" << endl;
    cin >> variableList[0];
    // cout << "dimension of a vector space over slot finite field:" << endl;
    // cin >> variableList[1];
    cout << "order of cyclotomic ring:" << endl;
    cin >> variableList[1];
    cout << "min bitsize of the ciphertext modulus in ciphertext:" << endl;
    cin >> variableList[2];
    // cout << "length of field vector:" << endl;
    // cin >> variableList[4];
    // cout << "number of values to be sorted:" << endl;
    // cin >> variableList[5];
    // cout << "number of experiments:" << endl;
    // cin >> variableList[6];
    cout << "type of circuit(B,T,U):" << endl;
    cin >> variableList[3];
    // cout << "level of consecutive comparison (depth):" << endl;
    // cin >> variableList[8];

    // int input;

    // cout << "Choose the function you want to run:" << endl;
    // cout << "1) comparison " << endl;
    // cout << "2) sorting " << endl;
    // cout << "3) minmax " << endl;
    // cin >> input;

    // Comparator c = compare_function(variableList);
    
    CircuitType type = UNI;
    if (!strcmp(variableList[3], "B"))
    {
        type = BI;
    }
    else if (!strcmp(variableList[3], "T"))
    {
        type = TAN;
    }
    else if (!strcmp(variableList[3], "U"))
    {
        type = UNI;
    }
    else
    {
        throw invalid_argument("Choose a valid circuit type (U for univariate, B for bivariate and T for Tan et al.\n");
    }

    int d=1; int expansion_len=1;

    //////////PARAMETER SET UP////////////////
    // Plaintext prime modulus
    unsigned long p = atol(variableList[0]);
    // Field extension degree
    // unsigned long d = atol(variableList[1]);
    // Cyclotomic polynomial - defines phi(m)
    unsigned long m = atol(variableList[1]);
    // Number of ciphertext prime bits in the modulus chain
    unsigned long nb_primes = atol(variableList[2]);
    // Number of columns of Key-Switching matix (default = 2 or 3)
    unsigned long cc = 3;
    
    cout << "Initialising context object..." << endl;
    // Intialise context
    auto context = ContextBuilder<BGV>()
            .m(m)
            .p(p)
            .r(1)
            .bits(nb_primes)
            .c(cc)
            .scale(6)
            .build();

    std::filebuf fb;
    fb.open ("context.txt",std::ios::out);
    std::ostream os(&fb);
    context.writeTo(os);
    fb.close();
    ///////////////////////////////////////////////////
    // READ STORED CONTEXT
    // std::ifstream ifs ("context_p_43_m_91342_nb_3000.txt", std::ifstream::in);
    // auto context=Context::readFrom(ifs);
    // unsigned long p=context.getP();
    // unsigned long m=context.getM();

    cout<<"p: "<<p<<endl;
    cout<<"m: "<<m<<endl;
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

    // //maximal number of digits in a number
    // unsigned long expansion_len = atol(variableList[4]);

    // Secret key management
    cout << "Creating secret key..." << endl;
    // Create a secret key associated with the context
    SecKey secret_key(context);
    // Generate the secret key
    secret_key.GenSecKey();
    cout << "Generating key-switching matrices..." << endl;
    helib::PubKey&m_pk = secret_key;
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
    Comparator comparator(context, type, d, expansion_len, secret_key, true);

    //   // get EncryptedArray
    //   const EncryptedArray& ea = m_context.getEA();
    //   //extract number of slots
    //   long nslots = ea.size();
    //   //get p
    //   unsigned long p = m_context.getP();
    //order of p
    unsigned long ord_p = context.getOrdP();


//////////////////////////////////////////////////////
    // t1 = {0,3,3,1,2,3,1,2,3,1,2,1,3,1,2,3}
    // q1 = {1,3,3,2,1,3,2,1,1,2,3,1,2,3,1,3}
    // t2 = {1,0,1,2,0,1,0}
    // q2 = {0,3,0,1,3,2,0,3,3}

    // A-0, G-1, T-2, C-3
    // (ATCA,ACAT); (CATT,ACT)
    // Arrays of 16 values; maxLen1, maxLen2 < 16
    // for max tlen,qlen=4 and w_match=5, macimum score possible=4*5=20. So,(p-1)/2>20=> p=43.
    
    // int t_arr[count][16]={{0,2,3,0,5,5,5,5,5,5,5,5,5,5,5,5},{3,0,2,2,5,5,5,5,5,5,5,5,5,5,5,5}};
    // int q_arr[count][16]={{0,3,0,2,6,6,6,6,6,6,6,6,6,6,6,6},{0,3,2,6,6,6,6,6,6,6,6,6,6,6,6}};
    // int tlen[count]={4,4};
    // int qlen[count]={4,3};
    // int h0_arr[count]={19,17};


    /////////////////////////////////////////////
    int count=nslots;
    int maxLen1=4; int maxLen2=4;         // unencrypted
    // int t_arr[count][maxLen1]={{0,2,3},{1,1,3}};  // {t1,t2}
    // int q_arr[count][maxLen2]={{0,3,2},{1,1,6}};  // {q1,q2}
    // int tlen[count]={3,3};
    // int qlen[count]={3,2};
    /////////////////////////////////////////////

    // Modify
    // int h0_arr[count]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    int cnt=count;  
    int t_arr[count][maxLen1]; // Final array to store arrays from the first section
    int q_arr[count][maxLen2]; // Final array to store arrays from the second section
    int tlen[count]; 
    int qlen[count];

    // string inputFile="number_nucl_file";
    ifstream inputFile("number_nucl_file");
    int kk=0;
    bool flg=false;

    if(inputFile.is_open()){
        string line;
        while (count > 0 && getline(inputFile, line)) {
            if (line.find("-----") != string::npos) {
              cout<<"tlen found ----------------"<<endl;
              flg=true;  
              break; // Stop reading when encountering a line containing "-----"
            }
            vector<int> digits = splitDigits(line);

            tlen[cnt-count]=digits.size();

            for(int i=0;i<maxLen1;i++)
            {
                if(i<digits.size())
                    t_arr[cnt-count][i]=digits[i];
                else
                    t_arr[cnt-count][i]=5;
            }
            // t_arr.push_back(digits);
            count--;
            kk++;
        }
        
       while(count>0)
       {
         for(int i=0;i<maxLen1;i++)
         {
            t_arr[cnt-count][i]=rand()%4;
         }
         // t_arr[cnt-count]=5;
         tlen[cnt-count]=maxLen1;
         count--;
       }

      if(flg==false){
        cout<<"flg==false"<<endl;
        while(getline(inputFile, line))
        {
            kk++;
            if (line.find("-----") != string::npos) {
                break; // Stop reading when encountering a line containing "-----"
            }
        }
      }
        cout<<"kk: "<<kk<<endl;

        count = cnt; // Reset count for the second section
        while (count > 0 && getline(inputFile, line)) {
            vector<int> digits = splitDigits(line);

            qlen[cnt-count]=digits.size();

            for(int i=0;i<maxLen2;i++)
            {
                if(i<digits.size())
                    q_arr[cnt-count][i]=digits[i];
                else
                    q_arr[cnt-count][i]=6;
            }
            // t_arr.push_back(digits);
            count--;
        }

       while(count>0)
       {
         for(int i=0;i<maxLen2;i++)
         {
            q_arr[cnt-count][i]=rand()%4;
         }
         // t_arr[cnt-count]=5;
         qlen[cnt-count]=maxLen2;
         count--;
       }


/*
       cout<<"t_arr: "<<endl;
        for(int i=0;i<cnt;i++)
        {
            for(int j=0;j<maxLen1;j++)
            {
                cout<<t_arr[i][j]<<" ";
            }
            cout<<endl;
        }

        cout << endl;
        cout<<"q_arr: "<<endl;

        for(int i=0;i<cnt;i++)
        {
            for(int j=0;j<maxLen2;j++)
            {
                cout<<q_arr[i][j]<<" ";
            }
            cout<<endl;
        }
*/

    inputFile.close();
    }

    count=cnt;

    Ctxt init_c(m_pk);
    // Modify
    Ctxt t[]={Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)};
    Ctxt q[]={Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)};
   
    cout<<"Initialized arrays!"<<endl;

    for(int i=0;i<maxLen1;i++)
    {
        // cout<<"i: "<<i<<endl;
        vector<ZZX> pol_t(nslots);
        for(int j=0;j<count;j++)
        {
            //cout<<"t_arr["<<j<<"]["<<i<<"]: "<<t_arr[j][i]<<endl;
            pol_t[j]=ZZX(INIT_MONO, 0, t_arr[j][i]);
        }
        //cout<<"pol_t: "<<pol_t<<endl;
        // t[i]=init_c;
        ea.encrypt(t[i], m_pk, pol_t); 

        helib::Ptxt<helib::BGV> temp_ptxt(context); 
        secret_key.Decrypt(temp_ptxt, t[i]);
        cout<<"t["<<i<<"]: "<<temp_ptxt<<endl;
        cout<<endl;
    }

    for(int i=0;i<maxLen2;i++)
    {
        vector<ZZX> pol_q(nslots);
        for(int j=0;j<count;j++)
        {
            pol_q[j]=ZZX(INIT_MONO, 0, q_arr[j][i]);
        }
        
        // q[i]=Ctxt(m_pk);
        ea.encrypt(q[i], m_pk, pol_q); 

        helib::Ptxt<helib::BGV> temp_ptxt(context); 
        secret_key.Decrypt(temp_ptxt, q[i]);
        cout<<"q["<<i<<"]: "<<temp_ptxt<<endl;
        cout<<endl;
    }

    Ctxt h0(m_pk);
    // Create the plaintext polynomials for the text and for the pattern
    vector<ZZX> pol_h0(nslots);

    // for(int i=0;i<nslots;i++)
    //     pol_h0[i]=h0_arr[i];

    // pol_h0[0]=h0_arr[0];
    // pol_h0[1]=h0_arr[1];
    // for(int j=2;j<count;j++)
    // {
    //     pol_h0[j]=ZZX(INIT_MONO, 0, h0_arr[j]);
    // }
    // ea.encrypt(h0, m_pk, pol_h0);          // h0
    
    Ctxt w_match_ctxt(m_pk); Ctxt w_mismatch_ctxt(m_pk); Ctxt w_ambig_ctxt(m_pk); Ctxt max_s_ctxt(m_pk);
    vector<ZZX> pol_w_match_ctxt(nslots); vector<ZZX> pol_w_mismatch_ctxt(nslots); vector<ZZX> pol_w_ambig_ctxt(nslots);
    for(int j=0;j<nslots;j++)
    {
        pol_w_match_ctxt[j]=ZZX(INIT_MONO, 0, 5);
        pol_w_mismatch_ctxt[j]=ZZX(INIT_MONO, 0, 3);
        pol_w_ambig_ctxt[j]=ZZX(INIT_MONO, 0, 1);     // Negative values not considered
    }
    ea.encrypt(w_match_ctxt, m_pk, pol_w_match_ctxt);
    ea.encrypt(w_mismatch_ctxt, m_pk, pol_w_mismatch_ctxt);
    ea.encrypt(w_ambig_ctxt, m_pk, pol_w_ambig_ctxt);
    cout<<"Got encrypted scores"<<endl;

    long int eb=5l, w_open=0l, w_extend=2l;
    long w=3;
    // Constants to be directly multiplied/added to ciphertext need to be of type 
    // long int, not int.

    ////////////////////// Score initialization ////////////////////
    auto start = high_resolution_clock::now();
    Ptxt<helib::BGV> new_plaintext_result(context);
    // Modify

    Ctxt H1[]={Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)};
    Ctxt F[]={Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)};
    Ctxt ctxt_cmp_result(m_pk); Ctxt ctxt_cmp_result_max_s(m_pk);

    ctxt_cmp_result_max_s = compare_ctxt(w_match_ctxt, w_mismatch_ctxt, comparator, m_pk, secret_key, ord_p);
    max_s_ctxt = assign_comparison_result(w_match_ctxt,w_mismatch_ctxt,ctxt_cmp_result_max_s,m_pk);
    ctxt_cmp_result_max_s = compare_ctxt(max_s_ctxt, w_ambig_ctxt, comparator, m_pk, secret_key, ord_p);
    max_s_ctxt = assign_comparison_result(max_s_ctxt,w_ambig_ctxt,ctxt_cmp_result_max_s,m_pk);
    secret_key.Decrypt(new_plaintext_result, max_s_ctxt);
    cout << "max_s_ctxt: " << new_plaintext_result << endl;

    /* 
    //// Score initialization
    H1[0]=h0;
    // H1[1]=max(H1[0]-w_open-w_extend,0)
    H1[1]=h0;
    H1[1]-=w_open; H1[1]-=w_extend; 
    // both CTs are fresh
    ctxt_cmp_result = compare_ctxt_num(h0,w_open+w_extend,comparator,ea,nslots,m_pk,secret_key,ord_p);
    // H1[1] is of depth (1+depth of comparison)
    ctxt_cmp_result*=-1l; ctxt_cmp_result+=1l; ctxt_cmp_result*=H1[1]; H1[1]=ctxt_cmp_result;
    secret_key.Decrypt(new_plaintext_result, H1[0]);
    cout << "H1[0] = " << new_plaintext_result << endl;
    secret_key.Decrypt(new_plaintext_result, H1[1]);
    cout << "H1[1] = " << new_plaintext_result << endl;

    for(int j=2;j<=maxLen2;j++)
    {
        cout<<"j: "<<j<<endl;
        // H1[j]=H1[j-1]; // H1[j]-=w_extend;
        ctxt_cmp_result = compare_ctxt_num(H1[j-1],w_extend,comparator,ea,nslots,m_pk,secret_key,ord_p);
           // Decrypt the modified ciphertext into a new plaintext
        Ptxt<helib::BGV> new_plaintext_result(context);
        secret_key.Decrypt(new_plaintext_result, ctxt_cmp_result);
        cout << "ctxt_cmp_result:  " << new_plaintext_result << endl;

        ctxt_cmp_result*=-1l; ctxt_cmp_result+=1l; 
        H1[j]=H1[j-1]; H1[j]-=w_extend; 
        ctxt_cmp_result*=H1[j]; H1[j]=ctxt_cmp_result;

        // Decrypt the modified ciphertext into a new plaintext
        // Ptxt<helib::BGV> new_plaintext_result(context);
        secret_key.Decrypt(new_plaintext_result, H1[j]);
        cout << "H1["<<j<<"] = " << new_plaintext_result << endl;
    }
    //  Ctxt H_x[][maxLen2+1]={{Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)},{Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)},{Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)},{Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)},{Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)}};
    //  Ctxt H_y[][maxLen2+1]={{Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)},{Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)},{Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)},{Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)},{Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk),Ctxt(m_pk)}};
 

    // for(int i=0;i<5;i++)
    // {
    //     for(int j=0;j<5;j++)
    //     {
    //         if(i==0 || j==0)
    //         {
    //             // H_x[i][j]=i;
    //             // H_y[i][j]=j;

    //             vector<ZZX> pol_i(nslots);
    //             vector<ZZX> pol_j(nslots);

    //             for(int k=0;k<count;k++)
    //             {
    //                 //cout<<"t_arr["<<j<<"]["<<i<<"]: "<<t_arr[j][i]<<endl;
    //                 pol_i[k]=ZZX(INIT_MONO, 0, i);
    //                 pol_j[k]=ZZX(INIT_MONO, 0, j);
    //             }
                
    //             ea.encrypt(H_x[i][j], m_pk, pol_i); 
    //             ea.encrypt(H_y[i][j], m_pk, pol_j);
    //         }
    //         else
    //         {
    //             vector<ZZX> pol_i(nslots);
    //             vector<ZZX> pol_j(nslots);

    //             for(int k=0;k<count;k++)
    //             {
    //                 //cout<<"t_arr["<<j<<"]["<<i<<"]: "<<t_arr[j][i]<<endl;
    //                 pol_i[k]=ZZX(INIT_MONO, 0, -1);
    //                 pol_j[k]=ZZX(INIT_MONO, 0, -1);
    //             }
    //             //cout<<"pol_t: "<<pol_t<<endl;
    //             // t[i]=init_c;
    //             ea.encrypt(H_x[i][j], m_pk, pol_i); 
    //             ea.encrypt(H_y[i][j], m_pk, pol_j);
    //         }
    //     }
    // }
    */  
    auto stop = high_resolution_clock::now();

    // Ptxt<helib::BGV> new_plaintext_result(context);
    secret_key.Decrypt(new_plaintext_result, F[0]);
    cout << "F[0]: " << new_plaintext_result << endl;

    auto duration = duration_cast<microseconds>(stop - start);
	// To get the value of duration use the count()
	// member function on the duration object
	cout <<"Time for score_initialization: "<< duration.count() << endl;

    ////////////////////////////////////////////////////////////////////////////
    // cout<<"qlen[0]: "<<qlen[0]<<endl;
    // cout<<"count: "<<count;

    // band size
    helib::Ptxt<helib::BGV> qlen_ptxt(context);
    cout<<"qlen_ptxt.size(): "<<qlen_ptxt.size()<<endl;
    
    for(int j=0;j<count;j++)
    {
        qlen_ptxt[j]=ZZX(INIT_MONO, 0, qlen[j]);
    }
    
    cout<<"qlen_ptxt: "<<qlen_ptxt<<endl;
    // t[i]=init_c;
    // ea.encrypt(qlen_ctxt, m_pk, qlen_pol_t); 
    // secret_key.Decrypt(new_plaintext_result, qlen_ctxt);
    // cout << "qlen_ctxt: " << new_plaintext_result << endl;

    /////////////////////////////////////////////////////////////////
    auto start_band = high_resolution_clock::now();
    // qlen={1,2,3,4}; max_s=Enc{1,2,3,4,5}
    // To be interpoalted - min(w,max(1,(int)((qlen*max_s-3)/3+1)))
    // [(-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1), (3, 2), (5, 2), (7, 2), (6, 2), (9, 2), (12, 2), (13, 2), (17, 2)]
    // For p=257, interpolation polynomial is -[13, 14, 159, 91, 236, 99, 204, 119, 234, 237, 182, 211, 1]
    // For p=43, interpolation polynomial is - [5, 3, 41, 13, 34, 29, 11, 30, 18, 37, 20, 17, 1]
    
    // Band size w=3
    // [(-2, 1), (-1, 1), (0, 1), (1, 1), (2, 1), (3, 2), (5, 2), (7, 3), (6, 3), (9, 3), (12, 3), (13, 3), (17, )]

    Ctxt myband(m_pk); Ctxt temp(m_pk); temp=max_s_ctxt; temp*=qlen_ptxt; 
    // temp-=3l;
    temp+=1l;

    // secret_key.Decrypt(new_plaintext_result, temp);
    // cout << "temp: " << new_plaintext_result << endl;
    Ctxt m1(m_pk); m1=temp; Ctxt m2(m_pk); m2=temp; Ctxt m3(m_pk); 
    Ctxt m4(m_pk);Ctxt m5(m_pk);Ctxt m6(m_pk);Ctxt m7(m_pk);Ctxt m8(m_pk);
    Ctxt m9(m_pk);Ctxt m10(m_pk);Ctxt m11(m_pk);Ctxt m12(m_pk);
    m2*=m2;  // m2 => x^2; depth=1
    m3=m2; m3*=m1;      // m3 => x^3; depth=2
    m4=m2; m4*=m4;      // m4 => x^4; depth=2
    m5=m3; m5*=m2;      // m5 => x^5; depth=3
    m6=m4; m6*=m2;      // m6 => x^6; depth=3
    m7=m4; m7*=m3;      // m7 => x^7; depth=3
    m8=m4; m8*=m4;      // m8 => x^8; depth=3
    m9=m7; m9*=m2;      // m9 => x^9; depth=4
    m10=m8; m10*=m2;    // m10 => x^10; depth=4
    m11=m7; m11*=m4;    // m11 => x^11; depth=4
    m12=m7; m12*=m5;     // m12 => x^12; depth=4

    // m12*=13l; m11*=14l; m10*=159l; m9*=91l; m8*=236l; m7*=99l; m6*=204l; m5*=119l;
    // m4*=234l; m3*=237l; m2*=182l; m1*=211l; 

    // p = 97
    // m12*=12l; m11*=65l; m10*=51l; m9*=35l; m8*=23l; m7*=8l; m6*=31l; m5*=34l;
    // m4*=53l; m3*=33l; m2*=24l; m1*=19l;

    // // p=43, band size=2
    // m12*=5l; m11*=3l; m10*=41l; m9*=13l; m8*=34l; m7*=29l; m6*=11l; m5*=30l;
    // m4*=18l; m3*=37l; m2*=20l; m1*=17l;

    // p=43, band size=3
    // m12*=13l; m11*=9l; m10*=18l; m9*=10l; m8*=9l; m7*=41l; m6*=2l; m5*=28l;
    // m4*=30l; m3*=25l; m2*=14l; m1*=16l;

    // m12*=20l; m11*=11l; m10*=7l; m9*=9l; m8*=18l; m7*=31l; m6*=26l; m5*=25l;
    // m4*=34l; m3*=35l; m2*=20l; m1*=17l;

   // m12*=34l; m11*=31l; m10*=8l; m9*=12l; m8*=29l; m7*=6l; m6*=5l; m5*=37l;
    // m4*=37l; m3*=4l; m2*=16l; m1*=1l;
m12*=6l; m11*=28l; m10*=30l; m9*=27l; m8*=35l; m7*=9l; m6*=40l; m5*=9l;
    m4*=12l; m3*=32l; m2*=19l; m1*=3l;

    // secret_key.Decrypt(new_plaintext_result, m12);
    // cout << "m12: " << new_plaintext_result << endl;
    // secret_key.Decrypt(new_plaintext_result, m11);
    // cout << "m11: " << new_plaintext_result << endl;
    // secret_key.Decrypt(new_plaintext_result, m10);
    // cout << "m10: " << new_plaintext_result << endl;
    // secret_key.Decrypt(new_plaintext_result, m9);
    // cout << "m9: " << new_plaintext_result << endl;
    // secret_key.Decrypt(new_plaintext_result, m8);
    // cout << "m8: " << new_plaintext_result << endl;
    // secret_key.Decrypt(new_plaintext_result, m7);
    // cout << "m7: " << new_plaintext_result << endl;
    // secret_key.Decrypt(new_plaintext_result, m6);
    // cout << "m6: " << new_plaintext_result << endl;
    // secret_key.Decrypt(new_plaintext_result, m5);
    // cout << "m5: " << new_plaintext_result << endl;
    // secret_key.Decrypt(new_plaintext_result, m4);
    // cout << "m4: " << new_plaintext_result << endl;
    // secret_key.Decrypt(new_plaintext_result, m3);
    // cout << "m3: " << new_plaintext_result << endl;
    // secret_key.Decrypt(new_plaintext_result, m2);
    // cout << "m2: " << new_plaintext_result << endl;
    // secret_key.Decrypt(new_plaintext_result, m1);
    // cout << "m1: " << new_plaintext_result << endl;
    // // secret_key.Decrypt(new_plaintext_result, m0);
    // // cout << "m0: " << new_plaintext_result << endl;

    myband=m12; myband+=m11; myband+=m10; myband+=m9; myband+=m8; myband+=m7; 
    myband+=m6; myband+=m5; myband+=m4; myband+=m3; myband+=m2; myband+=m1; myband+=22l; 
    helib::Ptxt<helib::BGV> myband_ptxt(context); 
    secret_key.Decrypt(myband_ptxt, myband);
    cout << "myband_ptxt: " << myband_ptxt << endl;

    for(int i=0;i<count;i++)
    {
        w=max((long)(new_plaintext_result[i]),w);
    }
    cout<<"w: "<<w<<endl;
    auto stop_band = high_resolution_clock::now();
    auto duration_band = duration_cast<microseconds>(stop_band - start_band);
	// To get the value of duration use the count()
	// member function on the duration object
	cout <<"Time for band finding: "<< duration_band.count() << endl;

    ////////////////////////////////////////////////////

    auto start_dp = high_resolution_clock::now();

    Ctxt max_i(m_pk); Ctxt max_j(m_pk); Ctxt gscore(m_pk); Ctxt max_ie(m_pk); Ctxt max_sc(m_pk);
    Ctxt max_off(m_pk);
    vector<ZZX> pol_max_i(nslots); vector<ZZX> pol_max_j(nslots); vector<ZZX> pol_gscore(nslots);
    vector<ZZX> pol_max_ie(nslots); vector<ZZX> pol_max_off(nslots);
    
    for(int j=0;j<count;j++)
    {
        pol_max_i[j]=ZZX(INIT_MONO, 0, -1);     // max_i=Enc(-1)
        pol_max_j[j]=ZZX(INIT_MONO, 0, -1);     // max_j=Enc(-1)
        pol_gscore[j]=ZZX(INIT_MONO, 0, 0);
        pol_max_ie[j]=ZZX(INIT_MONO, 0, -1);
        pol_max_off[j]=ZZX(INIT_MONO, 0, 0);
    }
    ea.encrypt(max_i, m_pk, pol_max_i);
    ea.encrypt(max_j, m_pk, pol_max_j);
    ea.encrypt(gscore, m_pk, pol_gscore);
    ea.encrypt(max_ie, m_pk, pol_max_ie);
    ea.encrypt(max_off, m_pk, pol_max_off);
    max_sc=h0;                      // max = h0

    int minq = 10000000;
    for (int l=0; l<count; l++) {
        if (qlen[l]< minq) minq = qlen[l];
    }
    minq -= 1;                      // minq=min(qlen[k])-1
    
    int beg=0, end=w+1;
    helib::Ptxt<helib::BGV> head(context); helib::Ptxt<helib::BGV> tail(context);
    for(int j=0;j<count;j++)
    {
        head[j]=0;
        // tail[j]=qlen[j];
        tail[j]=myband_ptxt[j]+1;
    }

    for(int i=0;i<maxLen1;i++)
    {
        auto start_i = high_resolution_clock::now();
        cout<<"######################################"<<endl;
        cout<<"i: "<<i<<endl;
        helib::Ptxt<helib::BGV> i_ptxt(context);
        for(int k=0;k<count;k++)
            i_ptxt[k]=i;

        Ctxt e(m_pk); Ctxt m(m_pk); Ctxt mj(m_pk);  Ctxt i_ctxt(m_pk);               // e,m,j=Enc(0)
        vector<ZZX> pol_e(nslots); vector<ZZX> pol_m(nslots); vector<ZZX> pol_mj(nslots); vector<ZZX> pol_i_ctxt(nslots);
        for(int j=0;j<count;j++)
        {
            pol_e[j]=ZZX(INIT_MONO, 0, 0);
            pol_m[j]=ZZX(INIT_MONO, 0, 0);
            pol_mj[j]=ZZX(INIT_MONO, 0, -1);
            pol_i_ctxt[j]=ZZX(INIT_MONO, 0, i);
        }
        ea.encrypt(e, m_pk, pol_e); ea.encrypt(m, m_pk, pol_m);ea.encrypt(mj, m_pk, pol_mj); ea.encrypt(i_ctxt, m_pk, pol_i_ctxt); 

        beg=max(beg,i-w); end=min(max(end,i+w+1),maxLen2);    // beg=max(beg,i-w); end=min(min(end,i+w+1),maxLen2);  
        cout<<"beg: "<<beg<<"; end: "<<end<<endl;

        helib::Ptxt<helib::BGV> comp(context); helib::Ptxt<helib::BGV> comp2(context);
        helib::Ptxt<helib::BGV> temp_ptxt(context); helib::Ptxt<helib::BGV> head_temp(context);
        
        head_temp=head; head_temp+=myband_ptxt;     
        comp = compare_ptxt(head_temp, i_ptxt, comparator, p); // comp=(head+myband<i);
        comp2=comp; comp2*=-1l; comp2+=1l;
        temp_ptxt=i_ptxt; temp_ptxt-=myband_ptxt; temp_ptxt*=comp;  // temp_ptxt=(i-myband)*comp
        comp2*=head;                // comp2 = (1-comp)*head
        head+=comp2; head+=temp_ptxt;   // head=max(head,i-myband)
        cout<<"head: "<<head<<endl;
     
        cout<<"qlen_ptxt: "<<qlen_ptxt<<endl;
        temp_ptxt=i_ptxt; temp_ptxt+=myband_ptxt; temp_ptxt+=1l;      // temp_ptxt=i+myband+1
        cout<<"temp_ptxt: "<<temp_ptxt<<endl;
        comp = compare_ptxt(tail, temp_ptxt, comparator, p);    // comp=(tail<i+myband+1)
        comp2=comp; comp2*=-1l; comp2+=1l;
        comp*=temp_ptxt; tail*=comp2; tail+=comp;                        // tail=max(tail,i+1+myband)
        comp=compare_ptxt(tail, qlen_ptxt, comparator, p);      // comp=(tail<qlen)
        comp2=comp; comp2*=-1l; comp2+=1l;
        tail*=comp; comp2*=qlen_ptxt; tail+=comp2;          // tail=min(max(tail,i+1+myband),qlen)
        cout<<"tail: "<<tail<<endl;
        // tail can be equal to qlen; If we make band as [head,tail)
        // instead of [head,tail],i.e, not consider j==qlen, then,
        // H1[qlen] will not be evaluated appropriately
        // SO, head<=j AND j<=tail is correct comparison and not
        // head<=j AND j<tail

        Ctxt h1(m_pk); Ctxt h_x_temp(m_pk); Ctxt h_x_temp_2(m_pk);

        if(beg==0)
        {
            h1=h0; h1-=w_open; h1-=((w_extend)*(i+1));
            Ctxt ctxt_cmp_result(m_pk);
            ctxt_cmp_result = compare_ctxt_num(w_open+w_extend*(i+1),h0,comparator,ea,nslots,m_pk,secret_key,ord_p);
            h1*= ctxt_cmp_result;
        }

        Ctxt zero_ctxt(m_pk); Ctxt tmp(m_pk); helib::Ptxt<helib::BGV> tmp_ptxt(comparator.m_context);

        // We compute SIMD. [head,tail) are the starting and ending columns
        // respectively for multiple sequence pairs, i.e, head is Plaintext
        // that encodes starting column of multiple sequence pairs.

        // [beg,end) are the minimum head and maximum tail amongst all sequence pairs.
        // We only iterate across i in [beg,end) columns and for each
        // i, check if it is within [head,tail) for each sequence pair
        // and only then compute  

        for(int j=beg; j<end;j++)
        {
            cout<<"-------------------------"<<endl;
            cout<<"j: "<<j<<endl;
            helib::Ptxt<helib::BGV> j_ptxt=create_plaintext_from_constant(j,count,comparator);
            Ctxt j_ctxt(m_pk); j_ctxt+=j_ptxt;
            
            comp=compare_ptxt(j_ptxt,head,comparator,p);     // j<head
            comp2=compare_ptxt(tail,j_ptxt,comparator,p);    // tail<j
            // Let comp=(j<head); comp2=(tail<j)
            helib::Ptxt<helib::BGV> comp3(context);
            comp3=comp2; comp3*=-1l; comp3+=1l;
            comp*=comp3;                                                                  
            comp+=comp2;                                     
            // comp = comp*(1-comp2)+comp2 => (j<head OR tail<j)

            Ctxt f(m_pk); f=F[j];                // f=F[j]

            helib::Ptxt<helib::BGV> temp_ptxt(context); 
            secret_key.Decrypt(temp_ptxt, t[i]);
            cout << "t["<<i<<"]: " << temp_ptxt << endl;
            secret_key.Decrypt(temp_ptxt, q[j]);
            cout << "q["<<j<<"]: " << temp_ptxt << endl;
            // Ctxt w_m(m_pk); w_m=get_match_score(t[i],q[j],w_match_ctxt,w_mismatch_ctxt,w_ambig_ctxt,p,comparator,ea,nslots,m_pk,secret_key,ord_p,count);
            cout<<"old H1["<<j<<"]: "<<endl;
            secret_key.Decrypt(temp_ptxt, H1[j]);
            cout << temp_ptxt << endl;

            Ptxt<helib::BGV> new_plaintext_result(comparator.m_context);

            // secret_key.Decrypt(new_plaintext_result, a);
            // cout << "a: " << new_plaintext_result << endl;
            // secret_key.Decrypt(new_plaintext_result, b);
            // cout << "b: " << new_plaintext_result << endl;

            Ctxt res(m_pk);
            new_plaintext_result=create_plaintext_from_constant(4,count,comparator);
            Ctxt c_ambig_1=equal_ctxt(t[i],new_plaintext_result,p,comparator,m_pk,secret_key);
            Ctxt c_ambig_2=equal_ctxt(q[j], new_plaintext_result, p,comparator,m_pk,secret_key);
            // c_ambig_1*=-1l; c_ambig_1+=1l;
            // c_ambig_2*=-1l; c_ambig_2+=1l;
            Ctxt c_ambig_3(m_pk); c_ambig_3=c_ambig_2; c_ambig_3*=-1l; c_ambig_3+=1l;
            c_ambig_1*=c_ambig_3; c_ambig_1+=c_ambig_2;  
            secret_key.Decrypt(new_plaintext_result, c_ambig_1);
            // cout << "c_ambig_1: " << new_plaintext_result << endl;      
            //c_ambig_1 = (a==4 OR b==4)

            // secret_key.Decrypt(new_plaintext_result, c_ambig_1);
            // cout << "c_ambig_1: " << new_plaintext_result << endl;

            Ctxt c_match(m_pk);
            c_match = equal_ctxt(t[i], q[j],p,comparator,m_pk,secret_key);
            secret_key.Decrypt(new_plaintext_result, c_match);
            cout << "c_match: " << new_plaintext_result << endl;

            // res = assign_comparison_result(w_match_ctxt,w_mismatch_ctxt,c_match,m_pk);
            // // secret_key.Decrypt(new_plaintext_result, res);
            // // cout << "res1: " << new_plaintext_result << endl;

            // res = assign_comparison_result(res,w_ambig_ctxt,c_ambig_1,m_pk);
            
            
            // M_temp_1=max(H1[j]-w_mismatch,0)
            Ctxt M_temp_1(m_pk); 
            M_temp_1=H1[j]; M_temp_1-=w_mismatch_ctxt;  // w_mismatch is subtracted
            ctxt_cmp_result = compare_ctxt(w_mismatch_ctxt,H1[j], comparator, m_pk, secret_key, ord_p);
            // M_temp_1 = assign_comparison_result(zero_ctxt,M_temp_1,ctxt_cmp_result,m_pk);    
            M_temp_1*=ctxt_cmp_result;
            secret_key.Decrypt(new_plaintext_result, M_temp_1);
            // cout << "M_temp_1: " << new_plaintext_result << endl;


            // M_temp_2 =max(H1[j]+w_match,0)
            Ctxt M_temp_2(m_pk); 
            M_temp_2=H1[j]; M_temp_2+=w_match_ctxt;
            ctxt_cmp_result = compare_ctxt(zero_ctxt, M_temp_2, comparator, m_pk, secret_key, ord_p);
            // M_temp_2 = assign_comparison_result(zero_ctxt,M_temp_2,ctxt_cmp_result,m_pk);
            M_temp_2*=ctxt_cmp_result;
            secret_key.Decrypt(new_plaintext_result, M_temp_2);
            // cout << "M_temp_2: " << new_plaintext_result << endl;


            // M_temp_2 = (t[i]==q[j]?max(H1[j]+w_match,0): max(H1[j]-w_mismatch))
            M_temp_2 = assign_comparison_result(M_temp_1,M_temp_2,c_match,m_pk);    
            secret_key.Decrypt(new_plaintext_result, M_temp_2);
            // cout << "new M_temp_2: " << new_plaintext_result << endl;


            // M_temp_3 = max(H1[j]-w_ambig,0)
            Ctxt M_temp_3(m_pk);
            M_temp_3=H1[j];M_temp_3-=w_ambig_ctxt;
            ctxt_cmp_result = compare_ctxt(w_ambig_ctxt,H1[j], comparator, m_pk, secret_key, ord_p);
            //M_temp_3 = assign_comparison_result(zero_ctxt,M_temp_3,ctxt_cmp_result,m_pk);    
            M_temp_3*=ctxt_cmp_result;
            secret_key.Decrypt(new_plaintext_result, M_temp_3);
            // cout << "M_temp_3: " << new_plaintext_result << endl;

            Ctxt M(m_pk); 
            M = assign_comparison_result(M_temp_2,M_temp_3,c_ambig_1,m_pk);    
            secret_key.Decrypt(temp_ptxt, M);
            cout << "M: " << temp_ptxt << endl;

            // H1[j]=comp*H1[j]+(1-comp)*h1
            secret_key.Decrypt(temp_ptxt, h1);
            cout<<"h1 to be assigned to H1[j]: "<<temp_ptxt<<endl;
            cout<<"comp: "<<comp<<endl;
            H1[j] = assign_comparison_result(h1,H1[j],comp,m_pk); 
            secret_key.Decrypt(temp_ptxt, H1[j]);
            cout << "new H1["<<j<<"]: " << temp_ptxt << endl;

            Ctxt h(m_pk);
            // // h =max(M,e,f)
            ctxt_cmp_result = compare_ctxt(M, e, comparator, m_pk, secret_key, ord_p);
            // h_x_temp=ctxt_cmp_result; 
            // h_x_temp_2=ctxt_cmp_result;h_x_temp_2*=-1l; h_x_temp_2+=1l;
            // h_x_temp*=(long)(i+1); // if M<e, choose e, i.e, left cell to H[i+1][j+1]
            // h_x_temp_2*=(long)(i);
            // h_x_temp+=h_x_temp_2;
            // H_x[i+1][j+1]=h_x_temp;
            // H_y[i+1][j+1]=j_ctxt;
            secret_key.Decrypt(temp_ptxt,ctxt_cmp_result);
            cout << "ctxt_cmp_result " << temp_ptxt << endl;
            h = assign_comparison_result(M,e,ctxt_cmp_result,m_pk);
            secret_key.Decrypt(temp_ptxt, h);
            cout << "h=max(M,e): " << temp_ptxt << endl; 
            

            ctxt_cmp_result = compare_ctxt(h, f, comparator, m_pk, secret_key, ord_p);
            // h_x_temp=ctxt_cmp_result; 
            // h_x_temp_2=ctxt_cmp_result;h_x_temp_2*=-1l; h_x_temp_2+=1l;
            // h_x_temp*=(long)(i); // if f is greatest, choose it i.e, cell above to H[i+1][j+1]
            // H_x[i+1][j+1]*=h_x_temp_2;
            // H_x[i+1][j+1]+=h_x_temp;

            // h_x_temp=ctxt_cmp_result; 
            // h_x_temp*=(long)(j+1);
            // H_y[i+1][j+1]*=h_x_temp_2;
            // H_y[i+1][j+1]+=h_x_temp;

            // secret_key.Decrypt(temp_ptxt, H_x[i+1][j+1]);
            // cout << "1 H_x["<<i+1<<"]["<<j+1<<"]: " << temp_ptxt << endl; 
            // secret_key.Decrypt(temp_ptxt, H_y[i+1][j+1]);
            // cout << "1 H_y["<<i+1<<"]["<<j+1<<"]: " << temp_ptxt << endl; 
                  
            h = assign_comparison_result(h,f,ctxt_cmp_result,m_pk);

            // if score is 0, H_x=i+1, H_y=j+1
            Ctxt i_cpy_ctxt(m_pk); Ctxt j_cpy_ctxt(m_pk);
            i_cpy_ctxt=i_ctxt; i_cpy_ctxt+=1l;
            j_cpy_ctxt=j_ctxt; j_cpy_ctxt+=1l;
            ctxt_cmp_result = equal_ctxt(h, zero_ctxt,p,comparator,m_pk,secret_key);
            // H_x[i+1][j+1] = assign_comparison_result(H_x[i+1][j+1],i_cpy_ctxt,ctxt_cmp_result,m_pk);
            // H_y[i+1][j+1] = assign_comparison_result(H_y[i+1][j+1],j_cpy_ctxt,ctxt_cmp_result,m_pk);

            // secret_key.Decrypt(temp_ptxt, H_x[i+1][j+1]);
            // cout << "H_x["<<i+1<<"]["<<j+1<<"]: " << temp_ptxt << endl; 
            // secret_key.Decrypt(temp_ptxt, H_y[i+1][j+1]);
            // cout << "H_y["<<i+1<<"]["<<j+1<<"]: " << temp_ptxt << endl; 
                        
            h1=h;
            secret_key.Decrypt(temp_ptxt, h1);
            cout << "h1=max(M,e,f): " << temp_ptxt << endl; 
            //////// right

            comp2=comp; comp2*=-1l; comp2+=1l;      // comp2=(head<=j AND j<=tail)
            ctxt_cmp_result = compare_ctxt(m, h1, comparator, m_pk, secret_key, ord_p);
            ctxt_cmp_result*=comp2;         // m<h1 AND head<=j AND j<=tail

            // // mj=(m<h1 AND head<=j AND j<=tail)?j:mj
            // max_j = assign_comparison_result(max_j,j_ptxt,ctxt_cmp_result,m_pk); 
            // max_i = assign_comparison_result(max_i,i_ptxt,ctxt_cmp_result,m_pk); 
            
            // secret_key.Decrypt(temp_ptxt, max_i);
            // cout << "row number i with max score m is max_i: " << temp_ptxt << endl; 
            // secret_key.Decrypt(temp_ptxt, max_j);
            // cout << "column number j with max score m is max_j: " << temp_ptxt << endl; 
            
            // m=(m<h1 AND head<=j AND j<=tail)?h1:m
            m = assign_comparison_result(m,h1,ctxt_cmp_result,m_pk,secret_key,comparator); 
            secret_key.Decrypt(temp_ptxt, m);
            cout << "largest score in row "<<i<<" so far is m: " << temp_ptxt << endl; 
            mj = assign_comparison_result(mj,j_ptxt,ctxt_cmp_result,m_pk); 
            secret_key.Decrypt(temp_ptxt, mj);
            cout << "largest score in row "<<i<<" so far is at mj: " << temp_ptxt << endl; 
            
            // temp = max(M-(w_open+w_extend),0)
            Ctxt temp(m_pk); 
            temp=M; temp-=(w_open+w_extend);
            
            // ctxt_cmp_result = compare_ctxt_num(M,w_open+w_extend,comparator,ea,nslots,m_pk,secret_key,ord_p);
            // temp = assign_comparison_result(temp_M_diff,zero_ctxt,ctxt_cmp_result,m_pk); 
            ctxt_cmp_result = compare_ctxt_num(w_open+w_extend,M,comparator,ea,nslots,m_pk,secret_key,ord_p);
            temp*=ctxt_cmp_result;

            // temp2 = max(f-e_del,0)
            Ctxt temp2(m_pk);
            temp2=f; temp2-=(w_extend);
            ctxt_cmp_result = compare_ctxt_num(w_extend,f,comparator,ea,nslots,m_pk,secret_key,ord_p);
            temp2*=ctxt_cmp_result;

            // temp3 = max((M-(w_open+w_extend),f-e_del),0)
            Ctxt temp3(m_pk);
            ctxt_cmp_result = compare_ctxt(temp,temp2, comparator, m_pk, secret_key, ord_p);
            temp3 = assign_comparison_result(temp,temp2,ctxt_cmp_result,m_pk); 
            
            // F[j]=(1-comp)*temp3+comp*F[j]
            F[j] = assign_comparison_result(temp3,F[j],comp,m_pk); 
            secret_key.Decrypt(temp_ptxt, F[j]);
            cout << "F["<<j<<"]: " << temp_ptxt << endl; 

            // temp2 = max(e-w_extend,0)
            ctxt_cmp_result = compare_ctxt_num(w_extend,e,comparator,ea,nslots,m_pk,secret_key,ord_p);
            temp2=e; temp2-=(w_extend);
            temp2*=ctxt_cmp_result;
            // ctxt_cmp_result = compare_ctxt_num(e,w_extend,comparator,ea,nslots,m_pk,secret_key,ord_p);
            // temp2 = assign_comparison_result(temp_e_diff,zero_ctxt,ctxt_cmp_result,m_pk); 

            // temp3 = max((M-(w_open+w_extend),e-e_del),0)
            ctxt_cmp_result = compare_ctxt(temp,temp2, comparator, m_pk, secret_key, ord_p);
            temp3 = assign_comparison_result(temp,temp2,ctxt_cmp_result,m_pk); 
            
            //e=(1-comp)*temp3+comp*e
            e = assign_comparison_result(temp3,e,comp,m_pk); 
            secret_key.Decrypt(temp_ptxt, e);
            cout << "e: " << temp_ptxt << endl; 

		/* 
            if(j>=minq)
            {
                /////gscore calculate
                cout<<"gscore calculate............"<<endl;
                // tail=min(max(tail,i+1+myband),qlen); so, tail<=qlen
                // If j+1==qlen, j<qlen
                comp=compare_ptxt(j_ptxt,head,comparator,p);    // j>head
                comp2=compare_ptxt(tail,j_ptxt,comparator,p);   // tail<j
                comp*=-1l; comp+=1l; comp2*=-1l; comp2+=1l;
                comp*=comp2;        // j<=head AND j<=tail
                helib::Ptxt<helib::BGV> comp3(context);
                helib::Ptxt<helib::BGV> j_2_ptxt(context);
                j_2_ptxt=j_ptxt; j_2_ptxt+=1;
                comp3 = equal_ptxt(j_2_ptxt,qlen_ptxt, p, comparator, m_pk, secret_key); // j+1==qlen
                comp*=comp3;                // j<=head AND j<=tail AND j+1==qlen
                // secret_key.Decrypt(temp_ptxt, comp);
                cout << "comp: " << comp << endl;  
                secret_key.Decrypt(temp_ptxt, h1);
                cout << "h1 init: " << temp_ptxt << endl; 
                secret_key.Decrypt(temp_ptxt, gscore);
                cout << "gscore_init: " << temp_ptxt << endl; 

                Ctxt comp_gscore(m_pk);
                comp_gscore = compare_ctxt(h1,gscore, comparator, m_pk, secret_key, ord_p);
                comp_gscore*=-1l; comp_gscore+=1l; // gscore<=h1
                // comp_gscore*=comp;          // gscore<=h1
                secret_key.Decrypt(temp_ptxt, comp_gscore);
                cout << "comp_gscore1: " << temp_ptxt << endl; 
                comp_gscore*=comp;          // j<=head AND j<=tail AND j+1==qlen AND gscore<=h1
                secret_key.Decrypt(temp_ptxt, comp_gscore);
                cout << "comp_gscore2: " << temp_ptxt << endl; 
                // xdouble division by 0 after this
                
                gscore = assign_comparison_result(gscore,h1,comp_gscore,m_pk,secret_key,comparator);
                // secret_key.Decrypt(temp_ptxt, gscore);
                // cout << "gscore1: " << temp_ptxt << endl;
            
                Ctxt i2_ctxt(m_pk); i2_ctxt=i_ctxt; i2_ctxt+=1l;
                max_ie = assign_comparison_result(max_ie,i2_ctxt,comp_gscore,m_pk); 
                // max_ie=(i+1)th row has the maximum gscore
                secret_key.Decrypt(temp_ptxt, gscore);
                cout << "gscore: " << temp_ptxt << endl; 
                secret_key.Decrypt(temp_ptxt, max_ie);
                cout << "max_ie: " << temp_ptxt << endl;

                cout<<"gscore calculate end............"<<endl;
            }
           */ 
         
        }

        // // helib::Ptxt<helib::BGV> comp(context); 
        // // comp=(qlen==maxLen2)
        // comp = equal_ptxt(qlen_ptxt,create_plaintext_from_constant(maxLen2, count, comparator), p, comparator, m_pk, secret_key); 
        // comp2=comp;
        // Ctxt h1_copy(m_pk); h1_copy=h1;
        // comp2*=-1l; comp2+=1l; 
        // H1[end]*=comp2;
        // h1_copy*=comp;
        // H1[end]+=h1_copy;            // H1[end]=comp*h1+(1-comp)*H1[end]
        Ptxt<helib::BGV>  temp_sc_ptxt(comparator.m_context);
        
        H1[end]=h1;
        secret_key.Decrypt(temp_sc_ptxt, H1[end]);
        cout<<"last H1["<<end<<"]: "<<temp_sc_ptxt<<endl;

        // F[end]*=comp2;               // F[end]=(qlen==maxLen2)*Enc(0)+(qlen!=maxLen2)*F[end]
        F[end]=zero_ctxt;
        secret_key.Decrypt(temp_sc_ptxt, F[end]);
        cout<<"last F["<<end<<"]: "<<temp_sc_ptxt<<endl;
        // OPtimization - If Enc(0) is involved in comparison,
        // we can just multiply one comparison component like above
        // and skip the zero component addition,i.e, +=comp*zero_ctxt
    
        //// maxscore calculate
        cout<<"Maxscore calculate............"<<endl;
        Ctxt comp_m(m_pk);   
        comp_m = compare_ctxt(max_sc,m,comparator,m_pk,secret_key,ord_p);
        secret_key.Decrypt(temp_sc_ptxt, m);
        cout<<"original m: "<<temp_sc_ptxt<<endl;
        secret_key.Decrypt(temp_sc_ptxt, max_sc);
        cout<<"original max_sc: "<<temp_sc_ptxt<<endl;
        max_sc = assign_comparison_result(max_sc,m,comp_m,m_pk); 
        secret_key.Decrypt(temp_sc_ptxt, max_sc);
        cout<<"max(m,max_sc): "<<temp_sc_ptxt<<endl;

        Ctxt i_1_ctxt(m_pk); i_1_ctxt=i_ctxt; i_1_ctxt+=1l;  
        Ctxt mj_1_ctxt(m_pk); mj_1_ctxt=mj; mj_1_ctxt+=1l; 
        // max_i, max_j now store index, as i=0 here actially 
        // corresponds to i=1 in matrix    
        max_i= assign_comparison_result(max_i,i_1_ctxt,comp_m,m_pk);
        max_j= assign_comparison_result(max_j,mj_1_ctxt,comp_m,m_pk);
        
        Ctxt comp_off(m_pk); comp_off=mj; comp_off-=i_ctxt;
        tmp = compare_ctxt(i_ctxt, mj, comparator, m_pk, secret_key, ord_p);
        comp_off*=tmp;
        secret_key.Decrypt(tmp_ptxt, max_i);
        cout<<"max_i: "<<tmp_ptxt<<endl;
        secret_key.Decrypt(tmp_ptxt, max_j);
        cout<<"max_j: "<<tmp_ptxt<<endl;
        secret_key.Decrypt(tmp_ptxt, comp_off);
        cout<<"|i-mj|: "<<tmp_ptxt<<endl;
/*         
        max_off= assign_comparison_result(max_off,comp_off,comp_m,m_pk);
        secret_key.Decrypt(tmp_ptxt, max_off);
        cout<<"max_off: "<<tmp_ptxt<<endl;
 */       
cout<<"Maxscore calculate end............"<<endl;

        auto stop_i = high_resolution_clock::now();
        auto duration_i = duration_cast<microseconds>(stop_i - start_i);
        // To get the value of duration use the count()
        // member function on the duration object
        cout <<"Time for DP iteration "<<i<<" is: "<< duration_i.count() << endl;

    }

    auto stop_dp = high_resolution_clock::now();
    auto duration_dp = duration_cast<microseconds>(stop_dp - start_dp);
	// To get the value of duration use the count()
	// member function on the duration object
	cout <<"Time for DP iterations: "<< duration_dp.count() << endl;
 
    // cout<<"H_x: "<<endl;
    // for(int i=0;i<maxLen1+1;i++)
    // {
    //     helib::Ptxt<helib::BGV> tmp_ptxt(comparator.m_context);
    //     for(int j=0;j<maxLen2+1;j++)
    //     {
    //         secret_key.Decrypt(tmp_ptxt, H_x[i][j]);
    //         cout<<"i: "<<i<<"; j: "<<j<<endl;
    //         cout<<tmp_ptxt<<endl<<endl;
    //     }
    //     cout<<endl;
    // }

    // cout<<"H_y: "<<endl;
    // for(int i=0;i<maxLen1+1;i++)
    // {
    //     helib::Ptxt<helib::BGV> tmp_ptxt(comparator.m_context);
    //     for(int j=0;j<maxLen2+1;j++)
    //     {
    //         secret_key.Decrypt(tmp_ptxt, H_y[i][j]);
    //         cout<<"i: "<<i<<"; j: "<<j<<endl;
    //         cout<<tmp_ptxt<<endl<<endl;
    //     }
    //     cout<<endl;
    // }

    return 0;
}
