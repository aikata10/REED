#include "seal/seal.h"
#include <time.h>
#include <stdio.h>


using namespace seal;
using namespace std;
void run(){
    // context generation
    uint64_t logScale = 54;
    uint64_t logQ0 = logScale + 6;
    uint64_t levels = 6;

	double scale = pow(2.0, logScale);

    vector<int> rnsBits;
	rnsBits.push_back(logQ0);
	for (size_t i = 0; i < levels-1; i++) 
	{
		rnsBits.push_back(logScale);
	}
	rnsBits.push_back(logQ0);// 'special modulus'
	
	size_t logn = 13;

	size_t poly_modulus_degree = 2 * (1 << logn);
	std::cout << "poly degree: " << poly_modulus_degree << std::endl;
	EncryptionParameters parms(scheme_type::ckks);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    //parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, rnsBits));
            parms.set_coeff_modulus({576460752340123649,9007199379521537,9007199362744321,9007199353307137,9007199350161409,9007199343869953,9007199338627073,9007199282003969});
	SEALContext context(parms);
	KeyGenerator keygen(context);

    auto secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    GaloisKeys gal_keys;
	vector<int> steps;
    uint64_t factorDim = 18;
	for (int step = 1; step < factorDim; step <<= 1)
	{
		steps.push_back(step);
		steps.push_back(-step);
	}
    keygen.create_galois_keys(steps, gal_keys);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    CKKSEncoder encoder(context);
    size_t slots = encoder.slot_count();

    // homomorphic operation of 
    // g(<w. v>) 
    // where g is logit function (Taylor approximation of it), and w is trained weight, v is encrypted input
    std::vector<double> wData= {-0.00741075,1.29493,0.259481,
	0.222632,-0.0370947,0.197086,-0.207146,0.0847975,
	0.199955,-0.0321777,-0.000963437,-0.00697528,
	-0.161998,-0.188798,-0.161069,0.099619,0.130145,
	-0.10143,-0.155333,};
    std::vector<double> vData = {0,0,1,0,0,1,1,1,1,0,1,0,1,0,0,0,1,1};
    size_t len = vData.size();

    // Encoding and Encryption of input
    Plaintext vEncoded;
    encoder.encode(vData, scale, vEncoded);
    Ciphertext vEncrypted;
    encryptor.encrypt(vEncoded, vEncrypted);
    
    clock_t tic = clock();

    // Perform inner product
    // 1. element-wise multiplication of w and v
    Plaintext wEncoded;
    encoder.encode(wData, scale, wEncoded);
    evaluator.multiply_plain_inplace(vEncrypted, wEncoded);
    
    // 2. rotate-and-add: 
    // - let y be the vector (w_0*v_0, w_1*v_1, ..., w_{L-1}*v_{L-1})
    // - let L be the length of  w and v; NOTE: L doesn't have to be a power of two, but for brevity assume so.
    // - left rotate by L/2 and add will give us (y_0 + y_{L/2},  y_1 + y_{L/2+1}, ..., y_{L/2-1} + y_{L-1})
    // - left rotate by L/4 and add will give us (y_0 + y_{L/2} + y_{L/4} + y_{3L/4},  y_1 + y_{L/2+1} + y_{L/4 + 1} + y_{3L/4+1}, ...)
    // - iterate until L/(L/2) and the first element will give us \su_0^{L-1} y_0; we do not care about the others
    Ciphertext rot;
    for (size_t i = factorDim/2; i >= 1; i/=2)
    {
        evaluator.rotate_vector(vEncrypted, i, gal_keys, rot);
        evaluator.add_inplace(vEncrypted, rot);
    }
    // Rescale ciphertext
    evaluator.rescale_to_next_inplace(vEncrypted);

    // 3. Perform logit function g(x) = 1 / (1+e^x)
    // g(x) ~ 1/2 + 1/4 x + -1/48 x^3 + 1/480 x^5 + -17/80640 x^7
    // Find polynomial basis: x, x^3, x^5, x^7
    Ciphertext x2, x3, x5, x4, x7;
    // x^2
    evaluator.multiply(vEncrypted, vEncrypted, x2);
    evaluator.relinearize_inplace(x2, relin_keys);
    evaluator.rescale_to_next_inplace(x2);
    //x^3
    evaluator.mod_switch_to_next_inplace(vEncrypted);
    evaluator.multiply(x2, vEncrypted, x3);
    evaluator.relinearize_inplace(x3, relin_keys);
    evaluator.rescale_to_next_inplace(x3);
    // x^4 -> x^7 (if we find x^5*x^2 it will use more depth)
    evaluator.multiply(x2, x2, x4);
    evaluator.relinearize_inplace(x4, relin_keys);
    evaluator.rescale_to_next_inplace(x4);
    //x^5
    evaluator.mod_switch_to_next_inplace(x2);
    evaluator.multiply(x2, x3, x5);
    evaluator.relinearize_inplace(x5, relin_keys);
    evaluator.rescale_to_next_inplace(x5);
    // x^4 -> x^7 (if we find x^5*x^2 it will use more depth)
    evaluator.multiply(x4, x3, x7);
    evaluator.relinearize_inplace(x7, relin_keys);
    evaluator.rescale_to_next_inplace(x7);
    // match the level
    evaluator.mod_switch_to_inplace(vEncrypted, x7.parms_id());
    evaluator.mod_switch_to_inplace(x3, x7.parms_id());

    // Now we multiply coefficient of approximate polynomial
    std::vector<double> polyCoeff = {1./2., 1./4., 0, -1./48., 0, 1./480., 0, -17./80640.};
    // set the scale of polyCoeff
    Plaintext plainCoeff;
    // as the scaling factor or each ciphretext differs, we use independent scaling factor to multiply each polyCoeff[i]
    encoder.encode(polyCoeff[1], vEncrypted.parms_id(), scale*scale/vEncrypted.scale(), plainCoeff);
    evaluator.multiply_plain_inplace(vEncrypted, plainCoeff);

    encoder.encode(polyCoeff[3], x3.parms_id(), scale*scale/x3.scale(), plainCoeff);
    evaluator.multiply_plain_inplace(x3, plainCoeff);
    x3.scale() = vEncrypted.scale();

    encoder.encode(polyCoeff[5], x5.parms_id(), scale*scale/x5.scale(), plainCoeff);
    evaluator.multiply_plain_inplace(x5, plainCoeff);
    x5.scale() = vEncrypted.scale();

    encoder.encode(polyCoeff[7], x7.parms_id(), scale*scale/x7.scale(), plainCoeff);
    evaluator.multiply_plain_inplace(x7, plainCoeff);
    x7.scale() = vEncrypted.scale();

    evaluator.add_inplace(vEncrypted, x3);
    evaluator.add_inplace(vEncrypted, x5);
    evaluator.add_inplace(vEncrypted, x7);
    evaluator.rescale_to_next_inplace(vEncrypted);

    encoder.encode(polyCoeff[0], vEncrypted.parms_id(), vEncrypted.scale(), plainCoeff);
    evaluator.add_plain_inplace(vEncrypted, plainCoeff);
    clock_t toc = clock();

    printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

    Plaintext resultPt;
    decryptor.decrypt(vEncrypted, resultPt);
    std::vector<double> resultVec;
    encoder.decode(resultPt, resultVec);
    std::cout << "The probability of cancer is: " << resultVec[0] << std::endl;

    // validation in plaintext
    double innerProd = 0;
    for (size_t i = 0; i < factorDim; i++)
    {
        innerProd += wData[i]*vData[i];
    }
    cout << "Logit fucntion: " << 1/(1+exp(-innerProd)) << endl;
    double &x = innerProd;
    cout << "approximation: " << 1./2. + 1./4.*x + -1./48. * x*x*x + 1./480. * x*x*x*x*x -17./80640. * x*x*x*x*x*x*x << endl;

    
}

int main(){
    run();
    return 0;
}
