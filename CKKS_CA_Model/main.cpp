#include <iostream>
#include <fstream>
#include <stdint.h>
#include <vector>
#include <set>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <limits.h>
#include <chrono>
#include <thread>

using namespace std;

#include "numb.cpp"

#include "context_class.cpp"
#include "polynomial_functions_class.cpp"
#include "data_types.cpp"
#include "Scheme.cpp"

int main(){
	int i;
	unsigned random_seed=124;

	cout << endl << "-------------------------------------------------------------------- START" << endl << endl;

	Context context;
	cout << "Context created" << endl;

	Key key(context, random_seed); 
	cout << "Secret Key created" << endl;

	/*
	Evk evk(context, key, random_seed);
	cout << "Evaluation Key created" << endl;
	*/
	
	HBM hbm(context, key, random_seed);
	cout << "HBM initialized with evk Keys" << endl;

	Scheme scheme(context);
	cout << "Scheme created" << endl;

	///////////////////   Testing Homomorphic Addition and Multiplication

	double random1, random2;
	random1 = (double)(rand()%1000)/100; random2 = (double)(rand()%1000)/100;



	double m_plain1 = random1;
	double m_plain2 = random2;

	Plaintext plaintext1 = scheme.encode_single(m_plain1);
	Plaintext plaintext2 = scheme.encode_single(m_plain2);
	cout << "Ciphertexts are encoded!" << endl;

	Ciphertext ciphertext1 = scheme.encrypt(key, plaintext1);
	Ciphertext ciphertext2 = scheme.encrypt(key, plaintext2);

	Ciphertext ciphertext11 = scheme.encrypt(key, plaintext1);
	Ciphertext ciphertext21 = scheme.encrypt(key, plaintext2);
	cout << "Ciphertexts are encrypted!" << endl << endl;
	
	clock_t tStart,tEnd;

	// ------------------------------------------------------
	// write ciphertexts to file (as cXX.h)



	// ------------------------------------------------------
	// Hom. Addition
	
	cout << "------------------------------ HEadd" << endl;
	tStart = clock();
	Ciphertext ct_add_result = scheme.HEadd(ciphertext1,ciphertext2);
	tEnd = clock();

	double m_decoded_add = 0;
	Plaintext plaintext_decrypted_add = scheme.decrypt(key, ct_add_result);
	m_decoded_add = scheme.decode_single(plaintext_decrypted_add);

	cout << "m_plain1 = " << m_plain1 << endl;
	cout << "m_plain2 = " << m_plain2 << endl;

	cout << "Homomorphic add result = " << m_decoded_add << endl;
	cout << "Plain add result       = " << m_plain1+m_plain2 << endl;
	printf("Time taken in SW: %.6fs\n", (double)(tEnd - tStart)/CLOCKS_PER_SEC);
	printf("Time taken in HW: %.1fs clock cycles \n", (double)(HW_time));
	
	cout << "Error = " << m_decoded_add - (m_plain1+m_plain2) << endl;

	if( abs(m_decoded_add-(m_plain1+m_plain2))>0.001 ){
		printf("Difference is large\n");
	}
	cout << endl;

	// ------------------------------------------------------
	// Hom. Multiplication
	
	cout << "------------------------------ HEmultiply" << endl;
	tStart = clock();
	//Ciphertext ciphertext3 = scheme.HEmultiply(ciphertext1, ciphertext2, evk);
	Ciphertext ciphertext3 = scheme.HEmultiply(ciphertext1, ciphertext2, hbm);
	tEnd = clock();

	Plaintext plaintext_decrypted_mult = scheme.decrypt(key, ciphertext3);
	double m_decoded_multiply = scheme.decode_single(plaintext_decrypted_mult);

	cout << "m_plain1 = " << m_plain1 << endl;
	cout << "m_plain2 = " << m_plain2 << endl;

	cout << "Homomorphic multiplication result = " << m_decoded_multiply << endl;
	cout << "Plain multiplication result = " << m_plain1*m_plain2 << endl;
	printf("Time taken in SW: %.6fs\n", (double)(tEnd - tStart)/CLOCKS_PER_SEC);
	printf("Time taken in HW: %.1fs clock cycles \n", (double)(HW_time));

	cout << "Error = " << m_decoded_multiply - m_plain1*m_plain2 << endl;
	
	if( abs(m_decoded_multiply-m_plain1*m_plain2)>0.001 ){
		printf("Difference is large\n");
	}

	cout << endl << "-------------------------------------------------------------------- END" << endl << endl;

	return 0;
}






