This repository runs Banded Smith Waterman (BSW) using BGV Homomorphic Encryption. 
The current code reads bunch of nucleotide sequences of length <=4 from file number_nucl_file and performs parallel BSW with band size=3 (equivalent to Smith Waterman) for these sequence pairs. 
The sequence pairs are referred to as (reference sequence,query sequence) and length of reference sequence>=length of query sequence, as is present in numer_nucl_file. 
The current code runs BSW with weights w_match=5, w_mismatch=-3, w_gap_extend=-2, w_gap_start=0.
helibTest/main.cpp is the main code.

Prerequisities:-
Install HElib library.
Need >100 GB memory for (4,4) BSW with band size=3.

Can run (1,1) BSW by making following changes in main file helibTest/main.cpp:-
1) In line 497, change to "maxLen1=1", "maxLen2=1"
2) In line 695, change to "w=1"
3) In line 876, 877 change multiplications of m12, m11,...,m1 to 0l
4) In line 907, change to "myband+=1l"  

Run instructions:- 
1) In folder helibTest/build/, remove CMakeCache.txt and run your own "cmake .."
2) Run "make"
3) Run code using "./main" and give inputs for parameters plaintext modulus, degree of cyclotomic polynomial, number of bits in 
ciphertext modulus as 41, 885080, 8200 respectively.
 

