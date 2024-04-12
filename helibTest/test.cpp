#include <iostream>

#include <helib/helib.h>

int main(int argc, char* argv[])
{
    /*  Example of BGV scheme  */

    // Plaintext prime modulus
    unsigned long p = 81667;
    // Cyclotomic polynomial - defines phi(m)
    unsigned long m = 32109;
    // Hensel lifting (default = 1)
    unsigned long r = 1;
    // Number of bits of the modulus chain
    unsigned long bits = 500;
    // Number of columns of Key-Switching matrix (default = 2 or 3)
    unsigned long c = 2;

    std::cout << "Initialising context object..." << std::endl;
    // Initialize context
    // This object will hold information about the algebra created from the
    // previously set parameters
    helib::Context context = helib::ContextBuilder<helib::BGV>()
            .m(m)
            .p(p)
            .r(r)
            .bits(bits)
            .c(c)
            .build();

    // Print the context
    context.printout();
    std::cout << std::endl;

    // Print the security level
    std::cout << "Security: " << context.securityLevel() << std::endl;

    // Secret key management
    std::cout << "Creating secret key..." << std::endl;
    // Create a secret key associated with the context
    helib::SecKey secret_key(context);
    // Generate the secret key
    secret_key.GenSecKey();
    std::cout << "Generating key-switching matrices..." << std::endl;
    // Compute key-switching matrices that we need
    helib::addSome1DMatrices(secret_key);

    // Public key management
    // Set the secret key (upcast: SecKey is a subclass of PubKey)
    const helib::PubKey& public_key = secret_key;

    // Get the EncryptedArray of the context
    const helib::EncryptedArray& ea = context.getEA();

    // Get the number of slot (phi(m))
    long nslots = ea.size();
    std::cout << "Number of slots: " << nslots << std::endl;

    // Create a vector of long with nslots elements
    helib::Ptxt<helib::BGV> ptxt(context);
    // Set it with numbers 0..nslots - 1
    // ptxt = [0] [1] [2] ... [nslots-2] [nslots-1]
    int random;
    for (int i = 0; i < ptxt.size(); ++i) {
        random = rand() % 10;
        ptxt[i] = random;
    }

    // Print the plaintext
    std::cout << "Initial Plaintext: " << ptxt << std::endl;

    // Create a ciphertext object
    helib::Ctxt ctxt(public_key);
    // Encrypt the plaintext using the public_key
    public_key.Encrypt(ctxt, ptxt);

    /********** Operations **********/
    // Ciphertext and plaintext operations are performed
    // "entry-wise".

    // Square the ciphertext
    // [0] [1] [2] [3] [4] ... [nslots-1]
    // -> [0] [1] [4] [9] [16] ... [(nslots-1)*(nslots-1)]
    ctxt.multiplyBy(ctxt);
    // Plaintext version
    //ptxt.multiplyBy(ptxt);
    //std::cout << "Plaintext Result  for multiplication: " << ptxt << std::endl;

    helib::Ptxt<helib::BGV> test_plaintext_result(context);
    secret_key.Decrypt(test_plaintext_result, ctxt);
    std::cout << "Decrypted Result for multiplication: " << test_plaintext_result << std::endl;

    // Divide the ciphertext by itself
    // To do this we must calculate the multiplicative inverse using Fermat's
    // Little Theorem.  We calculate a^{-1} = a^{p-2} mod p, where a is non-zero
    // and p is our plaintext prime.
    // First make a copy of the ctxt using copy constructor
    helib::Ctxt ctxt_divisor(ctxt);
    // Raise the copy to the exponent p-2
    // [0] [1] [4] ... [16] -> [0] [1] [1] ... [1]
    // Note: 0 is a special case because 0^n = 0 for any power n
    ctxt_divisor.power(2);
    // a^{p-2}*a = a^{-1}*a = a / a = 1;
    helib::Ptxt<helib::BGV> test2_plaintext_result(context);
    secret_key.Decrypt(test2_plaintext_result, ctxt_divisor);
    std::cout << "Decrypted Result for power: " << test2_plaintext_result << std::endl;
    //ctxt.multiplyBy(ctxt_divisor);
    ctxt = ctxt_divisor;
    // Plaintext version
    //helib::Ptxt<helib::BGV> ptxt_divisor(ptxt);
    //ptxt_divisor.power(2);
    //std::cout << "Plaintext Result power: " << ptxt_divisor << std::endl;
    //ptxt.multiplyBy(ptxt_divisor);
    //ptxt = ptxt_divisor;
    // [0] [1] [1] ... [1] [1] -> [0] [2] [2] ... [2] [2]
    ctxt += ctxt;
    //ptxt += ptxt;

    helib::Ptxt<helib::BGV> test3_plaintext_result(context);
    secret_key.Decrypt(test3_plaintext_result, ctxt);
    std::cout << "Decrypted Result for addition: " << test3_plaintext_result << std::endl;
    //std::cout << "Plaintext Result for addition: " << ptxt << std::endl;

    // Subtract it from itself (result should be 0)
    // i.e. [0] [0] [0] [0] ... [0] [0]
//    ctxt.divideBy2();
//    helib::Ptxt<helib::BGV> test4_plaintext_result(context);
//    secret_key.Decrypt(test4_plaintext_result, ctxt);
//    std::cout << "Decrypted Result for division by 2: " << test4_plaintext_result << std::endl;
    //ctxt -= ctxt;
    // Plaintext version
    //ptxt -= ptxt;

    // Create a plaintext for decryption
    helib::Ptxt<helib::BGV> plaintext_result(context);
    // Decrypt the modified ciphertext
    secret_key.Decrypt(plaintext_result, ctxt);

    //std::cout << "Operation: 2(a*a)/(a*a) - 2(a*a)/(a*a) = 0" << std::endl;
    // Print the decrypted plaintext
    // Should be [0] [0] [0] ... [0] [0]
    //std::cout << "Decrypted Result: " << plaintext_result << std::endl;
    // Print the plaintext version result, should be the same as the ctxt version
    //std::cout << "Plaintext Result: " << ptxt << std::endl;

    // We can also add constants
    // [0] [0] [0] ... [0] [0] -> [1] [1] [1] ... [1] [1]
    ctxt.addConstant(NTL::ZZX(70l));
    helib::Ptxt<helib::BGV> test4_plaintext_result(context);
    secret_key.Decrypt(test4_plaintext_result, ctxt);
    std::cout << "Decrypted Result for addition of 70 using addConstant: " << test4_plaintext_result << std::endl;
    // Plaintext version
    //ptxt.addConstant(NTL::ZZX(1l));

    // And multiply by constants
    // [1] [1] [1] ... [1] [1]
    // -> [1*1] [1*1] [1*1] ... [1*1] [1*1] = [1] [1] [1] ... [1] [1]
    ctxt *= 10l;
    helib::Ptxt<helib::BGV> test5_plaintext_result(context);
    secret_key.Decrypt(test5_plaintext_result, ctxt);
    std::cout << "Decrypted Result for multiplication of 10: " << test5_plaintext_result << std::endl;
    // Plaintext version
    //ptxt *= 1l;

    // We can also perform ciphertext-plaintext operations
    // ctxt = [1] [1] [1] ... [1] [1], ptxt = [1] [1] [1] ... [1] [1]
    // ctxt + ptxt = [2] [2] [2] ... [2] [2]
    // Note: the output of this is also a ciphertext
    ctxt += ptxt;
    helib::Ptxt<helib::BGV> test6_plaintext_result(context);
    secret_key.Decrypt(test6_plaintext_result, ctxt);
    std::cout << "Plaintext Result: " << ptxt << std::endl;
    std::cout << "Decrypted Result for adding plaintext value: " << test6_plaintext_result << std::endl;

    // Decrypt the modified ciphertext into a new plaintext
    helib::Ptxt<helib::BGV> new_plaintext_result(context);
    secret_key.Decrypt(new_plaintext_result, ctxt);

    std::cout << "Operation: Enc{(0 + 1)*1} + (0 + 1)*1" << std::endl;
    // Print the decrypted plaintext
    // Should be [2] [2] [2] ... [2] [2]
    std::cout << "Decrypted Result: " << new_plaintext_result << std::endl;

    return 0;
}
