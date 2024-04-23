
class Scheme:polynomial_functions{

	public:
	
	Context& context;

	Scheme(Context& context) : context(context) {} 
 
	Plaintext encode_single(double m_plain);

	Ciphertext encrypt(Key& key, Plaintext& plaintext);	

	Plaintext decrypt(Key& key, Ciphertext& ciphertext);

	double decode_single(Plaintext& plaintext);

	Ciphertext HEadd(Ciphertext& ciphertext1, Ciphertext& ciphertext2);	

	/*
    	Ciphertext HEmultiply(Ciphertext& ciphertext1, Ciphertext& ciphertext2, Evk& evk);
	Ciphertext relinearize(non_linearized_Ciphertext& ciphertext, Evk& evk);
	void chiplet_relinearize(int i, int chiplet_ctr, non_linearized_Ciphertext& nlc, Evk& hbm,uint64_t *ct_a_shares,uint64_t *ct_b_shares);
    */


	Ciphertext HEmultiply(Ciphertext& ciphertext1, Ciphertext& ciphertext2, HBM& hbm);
	Ciphertext relinearize(non_linearized_Ciphertext& ciphertext, HBM& hbm);
        void chiplet_relinearize(int i, int chiplet_ctr, non_linearized_Ciphertext& nlc, uint64_t *PRNG_a, uint64_t *HBM_b,uint64_t *ct_a_shares,uint64_t *ct_b_shares);

	Ciphertext modDown_inplace(uint64_t *shares_qp1, uint64_t *shares_qp2, non_linearized_Ciphertext& nlc);

	void rescaling_inplace(Ciphertext& ct, int level_in);
};


Plaintext Scheme::encode_single(double m_plain) {

	Plaintext plaintext;
	
	double m_mul_penc = m_plain * context.penc;

	for (int i = 0; i < L; i++) {
		uint64_t* op = plaintext.m + i * POLDEG;
		op[0] = m_mul_penc >= 0 ? m_mul_penc : context.qVec[i] + m_mul_penc;
	}	
	
	return plaintext;
}

double Scheme::decode_single(Plaintext& plaintext) {

	double m_plain = (double) plaintext.m[0]/context.penc;
	
	return m_plain;
}

Ciphertext Scheme::encrypt(Key& key, Plaintext& plaintext){

	Ciphertext c(context);	

	uint64_t *e0 = new uint64_t[POLDEG]();	
	uint64_t *e1 = new uint64_t[POLDEG]();	
	uint64_t *v = new uint64_t[POLDEG]();	
	uint64_t temp3, temp4;
	uint64_t *op1, *op2, *result;

	// Step1: sample small error polynomial e0
	for(int i=0; i<POLDEG; i++)
		e0[i] = rand()%8; 

	// Step2: sample small error polynomial e1
	for(int i=0; i<POLDEG; i++)
		e1[i] = rand()%8; 

	// Step3: sample sbinary polynomial v
	for(int i=0; i<POLDEG; i++){
		temp3 = rand()%2;
		temp4 = rand()%2;
		v[i] = temp3*temp4; 
	}


	bitreverse_sequential(plaintext.m);
	uint64_t* temp1 = new uint64_t[POLDEG]();
	uint64_t* temp2 = new uint64_t[POLDEG]();

	


	//c0_q = (v*b_pk + m_plain + e0)%firreducible;
	for(int i=0; i<L; i++){

		//Step1: Compute v*b_pk where b_pk is in NTT
		copy(v, v+POLDEG, temp1);
		dit_fwd_ntt_qi(context, temp1, i);

		op1 = key.b_pk_q_shares + i*POLDEG;
		op2 = temp1;
		result = c.c0_q_shares + i*POLDEG;
		diadic_multiplication_qi(context, result, op1, temp1, i);

		// Step2: temp2 <-- NTT(m_plain+e0)
		op1 = plaintext.m;
		op2 = e0;
		result = temp2;
		polynomial_addition_qi(context, result, op1, op2, i);
		dit_fwd_ntt_qi(context, result, i);

		// Step3: NTT(v*b_pk + m_plain + e0)
		op1 = c.c0_q_shares + i*POLDEG;
		op2 = temp2;
		result = c.c0_q_shares + i*POLDEG;
		polynomial_addition_qi(context, result, op1, op2, i);
	}


	//c1_q = (v*a_pk + e1)%firreducible;
	for(int i=0; i<L; i++){

		//Step1: Compute v*a_pk where b_pk is in NTT
		copy(v, v+POLDEG, temp1);
		dit_fwd_ntt_qi(context, temp1, i);

		op1 = key.a_pk_q_shares + i*POLDEG;
		op2 = temp1;
		result = c.c1_q_shares + i*POLDEG;
		diadic_multiplication_qi(context, result, op1, op2, i);

		//Step2: temp2 <-- NTT(e1)
		copy(e1, e1+POLDEG, temp2);
		dit_fwd_ntt_qi(context, temp2, i);

		// Step3: NTT(v*a_pk + e1)
		op1 = c.c1_q_shares + i*POLDEG;
		op2 = temp2;
		result = c.c1_q_shares + i*POLDEG;
		polynomial_addition_qi(context, result, op1, op2, i);
	}

	delete[] e0;
	delete[] e1;
	delete[] v;
	delete[] temp1;
	delete[] temp2;
	return c;
}


Plaintext Scheme::decrypt(Key& key, Ciphertext& ciphertext) {
	
	Plaintext plaintext;
	uint64_t *op1, *op2, *result;


	//m_prime = (c0_dec_q_shares[0] + s*c1_dec_q_shares[0])%firreducible

	op1 = key.s;
	op2 = ciphertext.c1_q_shares + 0;
	result = plaintext.m;
	//polynomial_multiplication_qi(context, result, op1, op2, 0);
	diadic_multiplication_qi(context, result, op1, op2, 0);
	
	op1 = ciphertext.c0_q_shares + 0;
	op2 = plaintext.m;
	result = plaintext.m;
	polynomial_addition_qi(context, result, op1, op2, 0);

	dif_inv_ntt_qi(context, result, 0);
	bitreverse_sequential(plaintext.m);

	return plaintext;
}

Ciphertext Scheme::HEadd(Ciphertext& ciphertext1, Ciphertext& ciphertext2) {

	Ciphertext ciphertext(context);	// ciphertext = HE(m1+m2)
	
	HW_time=N2;//Initial Read cost
	
	for(int i=0; i<Seq_Proc; i++){
		HW_time=HW_time+N2;
		
		#pragma omp parallel for 
		for(int chiplet_ctr=0; chiplet_ctr<Num_Chiplets;chiplet_ctr++){	
			uint64_t *op1, *op2, *op3, *op4, *result;	
			
			int t= i+Seq_Proc*chiplet_ctr;
			if(t<L){
			// c0+c0', c1+c1'	
			op1 = ciphertext1.c0_q_shares + t*POLDEG;
			op2 = ciphertext2.c0_q_shares + t*POLDEG;		
			result = ciphertext.c0_q_shares + t*POLDEG;
			polynomial_addition_qi(context, result, op1, op2, t);

			op3 = ciphertext1.c1_q_shares + t*POLDEG;
			op4 = ciphertext2.c1_q_shares + t*POLDEG;		
			result = ciphertext.c1_q_shares + t*POLDEG;
			polynomial_addition_qi(context, result, op3, op4, t);
			}
		}
	}



	cout << "DONE" << endl;
	cout << endl;

	return ciphertext;
}


		

Ciphertext Scheme::HEmultiply(Ciphertext& ciphertext1, Ciphertext& ciphertext2, HBM& hbm) {

	non_linearized_Ciphertext nlciphertext(context);	// ciphertext = HE(m1+m2)

	


	cout << "--------------- HEmultiply - Dyadic" << endl;

	// d0_q_shares[j] = (c0_q_shares[j]*c0_prime_q_shares[j])%firreducible
	HW_time=HW_time+2*N2;//Initial Read cost
	for(int i=0; i<Seq_Proc; i++){
		HW_time=HW_time+2*N2; // Computation Cost for four (two parallel) multiplications and one accumulation (MAC)
		
		#pragma omp parallel for 
		for(int chiplet_ctr=0; chiplet_ctr<Num_Chiplets;chiplet_ctr++){
			
			uint64_t *op1, *op2, *result;
			uint64_t *temp_poly = new uint64_t[POLDEG]();
	
			int t= i+Seq_Proc*chiplet_ctr;
			//
			if(t<L){
				// c0*c0', c0*c1'
				op1 = ciphertext1.c0_q_shares + t*POLDEG;
				op2 = ciphertext2.c0_q_shares + t*POLDEG;
				result = nlciphertext.d0_q_shares + t*POLDEG;
				diadic_multiplication_qi(context, result, op1, op2, t);
				
				op2 = ciphertext2.c1_q_shares + t*POLDEG;
				result = nlciphertext.d1_q_shares + t*POLDEG;
				diadic_multiplication_qi(context, result, op1, op2, t);
				
				// c1*c1', (c0*c1')+c1*c0'
				op1 = ciphertext1.c1_q_shares + t*POLDEG;
				op2 = ciphertext2.c0_q_shares + t*POLDEG;
				diadic_multiplication_qi(context, temp_poly, op1, op2, t);
				
				op2 = ciphertext2.c1_q_shares + t*POLDEG;
				result = nlciphertext.d2_q_shares + t*POLDEG;
				diadic_multiplication_qi(context, result, op1, op2, t);

				op1 = nlciphertext.d1_q_shares + t*POLDEG;
				op2 = temp_poly;
				result = nlciphertext.d1_q_shares + t*POLDEG;
				polynomial_addition_qi(context, result, op1, op2, t);			
			}
			delete[] temp_poly;
		}
	}

	

	return relinearize(nlciphertext, hbm);
}



int flag=0;

void Scheme::chiplet_relinearize(int i, int chiplet_ctr, non_linearized_Ciphertext& nlc, uint64_t *PRNG_a, uint64_t *HBM_b,uint64_t *ct_a_shares,uint64_t *ct_b_shares) {

	uint64_t *temp_poly = new uint64_t[POLDEG]();
	uint64_t *temp_poly1 = new uint64_t[POLDEG]();
	uint64_t *temp_poly2 = new uint64_t[POLDEG]();
	uint64_t *op1, *op2, *result1, *result2;
	int t= i+Seq_Proc*chiplet_ctr;
			
			
	if(t<L){	
		// Compute Inv-NTT of d2. The shares will be used for key-switching.
		dif_inv_ntt_qi(context, nlc.d2_q_shares + t*POLDEG, t);	
	}
			
	//SYNC to ensure INTT has finished across all threads
	#pragma omp flush (flag)
	flag=flag+1;
	bool all_flag_one=1;
	while(flag<Num_Chiplets){
		;
	}


	//Communication-Computation parallel loop taking results from all chiplets
	
	for(int comm=0; comm<Num_Chiplets;comm++){
		
		
		for(int j=0; j<Seq_Proc; j++){	
			
			
		int t_c= Seq_Proc*comm+i;
		int t_j=j+Seq_Proc*chiplet_ctr;

		if(t_c<L){
			//Cost of NTT+MAC working in parallel as well as pipeline
			if(chiplet_ctr==0) HW_time=HW_time+N2; 
			
			//DEBUG:cout << " --- " << comm << " " << t_j << " " << t_c << endl;
			std::copy(nlc.d2_q_shares + t_c*POLDEG, nlc.d2_q_shares + t_c*POLDEG + POLDEG,
						temp_poly);

			dit_fwd_ntt_qi(context, temp_poly, t_j);
			
			
			//op1 = hbm.b_evk_q_shares + t_c*POLDEG*(L+1) + t_j*POLDEG;
			//op2 = hbm.a_evk_q_shares + t_c*POLDEG*(L+1) + t_j*POLDEG;	

			op1 = HBM_b + t_c*POLDEG*(Seq_Proc) + j*POLDEG;
			op2 = PRNG_a + t_c*POLDEG*(Seq_Proc) + j*POLDEG;
			
			diadic_multiplication_qi(context, temp_poly1, temp_poly, op1, t_j);
			diadic_multiplication_qi(context, temp_poly2, temp_poly, op2, t_j);


			op1 = ct_b_shares + j*POLDEG;
			op2 = ct_a_shares + j*POLDEG;
			result1 = ct_b_shares + j*POLDEG;
			result2 = ct_a_shares + j*POLDEG;

			polynomial_addition_qi(context, result1, op1, temp_poly1, t_j);
			polynomial_addition_qi(context, result2, op2, temp_poly2, t_j);
		
			}
		}
	}
	delete [] temp_poly;
	delete [] temp_poly1;
	delete [] temp_poly2;
	#pragma omp flush (flag)
	flag=0;
}


Ciphertext Scheme::relinearize(non_linearized_Ciphertext& nlc, HBM& hbm) {

	cout << "--------------- HEmultiply - relinearize" << endl;
        
	uint64_t *ct_a_shares = new uint64_t[(L + 1)*POLDEG]();
	uint64_t *ct_b_shares = new uint64_t[(L + 1)*POLDEG]();
	

	// Initialize ct_a_shares and ct_b_shares to 0s
	for(int i=0; i<(L + 1)*POLDEG; i++){
		ct_a_shares[i] = 0;
		ct_b_shares[i] = 0;
	}

	cout << "--------------- HEmultiply - relinearize-ct_b calculation" << endl;
	cout << "--------------- HEmultiply - relinearize-ct_a calculation" << endl;

	// Compute ct_b_qi = sum( d2_share0*evkb0 + d2_share1*evkb1 + ...) mod q0, q1, ..., p 
	// Compute ct_a_qi = sum( d2_share0*evka0 + d2_share1*evka1 + ...) mod q0, q1, ..., p 

	for(int i=0; i<Seq_Proc; i++){
		HW_time=HW_time+3*N2;	//INTT cost + half NTT cost
		#pragma omp parallel for  
		for(int chiplet_ctr=0; chiplet_ctr<Num_Chiplets;chiplet_ctr++){

			// THIS LOOP invoked parallel Chiplet execution with respective HBM-data
			    uint64_t *PRNG_a = hbm.a_evk_q_shares[chiplet_ctr];
				uint64_t *HBM_b = hbm.b_evk_q_shares[chiplet_ctr];
				chiplet_relinearize(i,chiplet_ctr,nlc,PRNG_a,HBM_b,ct_a_shares+Seq_Proc*chiplet_ctr*POLDEG,ct_b_shares+Seq_Proc*chiplet_ctr*POLDEG);
		}
	
	}



	cout << "--------------- HEmultiply - relinearize-modDown_inplace (a,b)" << endl;
	Ciphertext ciphertext=modDown_inplace(ct_b_shares,ct_a_shares,nlc);

    
	
	cout << "--------------- HEmultiply - relinearize-rescaling_inplace" << endl;

	rescaling_inplace(ciphertext, L-1);

	cout << "--------------- HEmultiply - relinearize-final result" << endl;



	cout << "DONE" << endl;
	cout << endl;


	return ciphertext;
}



Ciphertext Scheme::modDown_inplace(uint64_t *shares_qp1, uint64_t *shares_qp2, non_linearized_Ciphertext& nlc) {

	
	Ciphertext ciphertext(context);
	// Inv-NTT of the last share
	HW_time=HW_time+3*N2;	//INTT cost + half NTT cost
	dif_inv_ntt_pi(context, shares_qp1 + L*POLDEG, 0);
	
	for(int i=0; i<Seq_Proc; i++){
		HW_time=HW_time+N2;	//NTT cost
		
		#pragma omp parallel for  
		for(int chiplet_ctr=0; chiplet_ctr<Num_Chiplets;chiplet_ctr++){
		
				
			int t= i+Seq_Proc*chiplet_ctr;
			
			if(t<L){
				uint64_t *op1, *op2, *result;
				uint64_t *temp_poly = new uint64_t[POLDEG]();

				std::copy(shares_qp1+L*POLDEG, shares_qp1+L*POLDEG+POLDEG,
						  temp_poly);
				dit_fwd_ntt_qi(context, temp_poly, t);
				
				op1 = shares_qp1+t*POLDEG;
				op2 = temp_poly;
				result = op1;
				
				polynomial_subtraction_qi(context, result, op1, op2, t);
				for(int j=0; j<POLDEG; j++)
					mulMod(result[j], context.P_inv_mod_qi[t], result[j], context.qVec[t]);
				polynomial_addition_qi(context, ciphertext.c0_q_shares + t*POLDEG, result, nlc.d0_q_shares + t*POLDEG, t);
				
				delete [] temp_poly;
		
			}
		}
	}
	
	// Inv-NTT of the last share
	HW_time=HW_time+3*N2;	//INTT cost + half NTT cost
	dif_inv_ntt_pi(context, shares_qp2 + L*POLDEG, 0);
	
	for(int i=0; i<Seq_Proc; i++){
		HW_time=HW_time+N2;	//NTT cost
		
		#pragma omp parallel for  
		for(int chiplet_ctr=0; chiplet_ctr<Num_Chiplets;chiplet_ctr++){
		
				
			int t= i+Seq_Proc*chiplet_ctr;
			
			if(t<L){
				uint64_t *op1, *op2, *result;
				uint64_t *temp_poly = new uint64_t[POLDEG]();

				std::copy(shares_qp2+L*POLDEG, shares_qp2+L*POLDEG+POLDEG,
						  temp_poly);
				dit_fwd_ntt_qi(context, temp_poly, t);
				
				op1 = shares_qp2+t*POLDEG;
				op2 = temp_poly;
				result = op1;
				
				polynomial_subtraction_qi(context, result, op1, op2, t);
				for(int j=0; j<POLDEG; j++)
					mulMod(result[j], context.P_inv_mod_qi[t], result[j], context.qVec[t]);
				polynomial_addition_qi(context, ciphertext.c1_q_shares + t*POLDEG, result, nlc.d1_q_shares + t*POLDEG, t);
				
				delete [] temp_poly;
		
			}
		}
	}

	return ciphertext;
		
}


void Scheme::rescaling_inplace(Ciphertext& ct, int level_in) {
	
	

	HW_time=HW_time+3*N2;	//INTT cost + half NTT cost
	// Inv-NTT of the last share of c0
	dif_inv_ntt_qi(context, ct.c0_q_shares+level_in*POLDEG, level_in);
	
	for(int i=0; i<Seq_Proc; i++){
		HW_time=HW_time+N2;	//NTT cost + Mul+Sub cost
		
		#pragma omp parallel for  
		for(int chiplet_ctr=0; chiplet_ctr<Num_Chiplets;chiplet_ctr++){
		
				
			int t= i+Seq_Proc*chiplet_ctr;
			
			if(t<level_in){
				uint64_t *op1, *op2, *result;
				uint64_t *temp_poly = new uint64_t[POLDEG]();

				std::copy(ct.c0_q_shares+level_in*POLDEG, ct.c0_q_shares+level_in*POLDEG+POLDEG,
						  temp_poly);
				dit_fwd_ntt_qi(context, temp_poly, t);
				
				op1 = ct.c0_q_shares+t*POLDEG;
				op2 = temp_poly;
				result = op1;
				polynomial_subtraction_qi(context, result, op1, op2, t);	
				
				for(int j=0; j<POLDEG; j++)
					mulMod(result[j], context.ql_inv_mod_qj[level_in][t], result[j], context.qVec[t]);
					
				delete [] temp_poly;
				}
			}
		}
		
	
	HW_time=HW_time+3*N2;	//INTT cost + half NTT cost
	// Inv-NTT of the last share of c1
	dif_inv_ntt_qi(context, ct.c1_q_shares+level_in*POLDEG, level_in);
	
	for(int i=0; i<Seq_Proc; i++){
		HW_time=HW_time+N2;	//NTT cost + Mul+Sub cost
		
		#pragma omp parallel for  
		for(int chiplet_ctr=0; chiplet_ctr<Num_Chiplets;chiplet_ctr++){
		
				
			int t= i+Seq_Proc*chiplet_ctr;
			
			if(t<level_in){
				uint64_t *op1, *op2, *result;
				uint64_t *temp_poly = new uint64_t[POLDEG]();

				std::copy(ct.c1_q_shares+level_in*POLDEG, ct.c1_q_shares+level_in*POLDEG+POLDEG,
						  temp_poly);
				dit_fwd_ntt_qi(context, temp_poly, t);
				
				op1 = ct.c1_q_shares+t*POLDEG;
				op2 = temp_poly;
				result = op1;
				polynomial_subtraction_qi(context, result, op1, op2, t);	
				
				for(int j=0; j<POLDEG; j++)
					mulMod(result[j], context.ql_inv_mod_qj[level_in][t], result[j], context.qVec[t]);
				delete [] temp_poly;
			}
		}
		
	}

	
}	
	


/*
Ciphertext Scheme::HEmultiply(Ciphertext& ciphertext1, Ciphertext& ciphertext2, Evk& evk) {

	non_linearized_Ciphertext nlciphertext(context);	// ciphertext = HE(m1+m2)

	


	cout << "--------------- HEmultiply - Dyadic" << endl;

	// d0_q_shares[j] = (c0_q_shares[j]*c0_prime_q_shares[j])%firreducible
	HW_time=HW_time+2*N2;//Initial Read cost
	for(int i=0; i<Seq_Proc; i++){
		HW_time=HW_time+2*N2; // Computation Cost for four (two parallel) multiplications and one accumulation (MAC)
		
		#pragma omp parallel for 
		for(int chiplet_ctr=0; chiplet_ctr<Num_Chiplets;chiplet_ctr++){
			
			uint64_t *op1, *op2, *result;
			uint64_t *temp_poly = new uint64_t[POLDEG]();
	
			int t= i+Seq_Proc*chiplet_ctr;
			//
			if(t<L){
				// c0*c0', c0*c1'
				op1 = ciphertext1.c0_q_shares + t*POLDEG;
				op2 = ciphertext2.c0_q_shares + t*POLDEG;
				result = nlciphertext.d0_q_shares + t*POLDEG;
				diadic_multiplication_qi(context, result, op1, op2, t);
				
				op2 = ciphertext2.c1_q_shares + t*POLDEG;
				result = nlciphertext.d1_q_shares + t*POLDEG;
				diadic_multiplication_qi(context, result, op1, op2, t);
				
				// c1*c1', (c0*c1')+c1*c0'
				op1 = ciphertext1.c1_q_shares + t*POLDEG;
				op2 = ciphertext2.c0_q_shares + t*POLDEG;
				diadic_multiplication_qi(context, temp_poly, op1, op2, t);
				
				op2 = ciphertext2.c1_q_shares + t*POLDEG;
				result = nlciphertext.d2_q_shares + t*POLDEG;
				diadic_multiplication_qi(context, result, op1, op2, t);

				op1 = nlciphertext.d1_q_shares + t*POLDEG;
				op2 = temp_poly;
				result = nlciphertext.d1_q_shares + t*POLDEG;
				polynomial_addition_qi(context, result, op1, op2, t);			
			}
			delete[] temp_poly;
		}
	}

	

	return relinearize(nlciphertext, evk);
}

Ciphertext Scheme::relinearize(non_linearized_Ciphertext& nlc, Evk& evk) {

	cout << "--------------- HEmultiply - relinearize" << endl;
        
	uint64_t *ct_a_shares = new uint64_t[(L + 1)*POLDEG]();
	uint64_t *ct_b_shares = new uint64_t[(L + 1)*POLDEG]();
	

	// Initialize ct_a_shares and ct_b_shares to 0s
	for(int i=0; i<(L + 1)*POLDEG; i++){
		ct_a_shares[i] = 0;
		ct_b_shares[i] = 0;
	}

	

	cout << "--------------- HEmultiply - relinearize-ct_b calculation" << endl;
	cout << "--------------- HEmultiply - relinearize-ct_a calculation" << endl;

	// Compute ct_b_qi = sum( d2_share0*evkb0 + d2_share1*evkb1 + ...) mod q0, q1, ..., p 
	// Compute ct_a_qi = sum( d2_share0*evka0 + d2_share1*evka1 + ...) mod q0, q1, ..., p 

	for(int i=0; i<Seq_Proc; i++){
		HW_time=HW_time+3*N2;	//INTT cost + half NTT cost
		#pragma omp parallel for  
		for(int chiplet_ctr=0; chiplet_ctr<Num_Chiplets;chiplet_ctr++){
			    chiplet_relinearize(i,chiplet_ctr,nlc,evk,ct_a_shares+Seq_Proc*chiplet_ctr*POLDEG,ct_b_shares+Seq_Proc*chiplet_ctr*POLDEG);
		}
		cout << "--------------- " << endl;
	}



	cout << "--------------- HEmultiply - relinearize-modDown_inplace (a,b)" << endl;
	Ciphertext ciphertext=modDown_inplace(ct_b_shares,ct_a_shares,nlc);

    
	
	cout << "--------------- HEmultiply - relinearize-rescaling_inplace" << endl;

	rescaling_inplace(ciphertext, L-1);

	cout << "--------------- HEmultiply - relinearize-final result" << endl;



	cout << "DONE" << endl;
	cout << endl;


	return ciphertext;
}





void Scheme::chiplet_relinearize(int i, int chiplet_ctr, non_linearized_Ciphertext& nlc, Evk &evk,uint64_t *ct_a_shares,uint64_t *ct_b_shares) {

	
			int t= i+Seq_Proc*chiplet_ctr;
			
			
			if(t<L){	
				// Compute Inv-NTT of d2. The shares will be used for key-switching.
				dif_inv_ntt_qi(context, nlc.d2_q_shares + t*POLDEG, t);	
			}
			
			//SYNC to ensure INTT has finished across all threads
			#pragma omp flush (flag)
			flag=flag+1;
			bool all_flag_one=1;
			while(flag<Num_Chiplets){
				;
			}

				
				uint64_t *temp_poly = new uint64_t[POLDEG]();
				uint64_t *temp_poly1 = new uint64_t[POLDEG]();
				uint64_t *temp_poly2 = new uint64_t[POLDEG]();
				uint64_t *op1, *op2, *result1, *result2;
				
				
				
				
				//Communication-Computation parallel loop taking results from all chiplets
				
				for(int comm=0; comm<Num_Chiplets;comm++){
					
					
				   for(int j=0; j<Seq_Proc; j++){	
						
						
					int t_c= Seq_Proc*comm+i;
					int t_j=j+Seq_Proc*chiplet_ctr;

					if(t_c<L){
						//Cost of NTT+MAC working in parallel as well as pipeline
						if(chiplet_ctr==0) HW_time=HW_time+N2; 
						
						//DEBUG:cout << " --- " << comm << " " << t_j << " " << t_c << endl;
						std::copy(nlc.d2_q_shares + t_c*POLDEG, nlc.d2_q_shares + t_c*POLDEG + POLDEG,
							  	  temp_poly);

						dit_fwd_ntt_qi(context, temp_poly, t_j);
						
						
						op1 = evk.b_evk_q_shares + t_c*POLDEG*(L+1) + t_j*POLDEG;
						op2 = evk.a_evk_q_shares + t_c*POLDEG*(L+1) + t_j*POLDEG;	
						diadic_multiplication_qi(context, temp_poly1, temp_poly, op1, t_j);
						diadic_multiplication_qi(context, temp_poly2, temp_poly, op2, t_j);


						op1 = ct_b_shares + j*POLDEG;
						op2 = ct_a_shares + j*POLDEG;
						result1 = ct_b_shares + j*POLDEG;
						result2 = ct_a_shares + j*POLDEG;

						

						polynomial_addition_qi(context, result1, op1, temp_poly1, t_j);
						polynomial_addition_qi(context, result2, op2, temp_poly2, t_j);
					
						}
					}
				}
				delete [] temp_poly;
				delete [] temp_poly1;
				delete [] temp_poly2;
	    		#pragma omp flush (flag)
	    		flag=0;
}
*/
