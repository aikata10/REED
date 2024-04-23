// key = [s, b_pk, a_pk]
class Key:polynomial_functions{

	public:
	uint64_t *s;				// Pointer to Secret key; Note: s is not in RNS as it is composed of small coefficients. 
	uint64_t *b_pk_q_shares;	// Pointer to RNS polynomial vector b_pk
	uint64_t *a_pk_q_shares;	// Pointer to RNS polynomial vector a_pk

	Key(Context& context, unsigned seed);	// constructor for creating pk-sk pair
	~Key();
};


Key::Key(Context& context, unsigned seed) {
	s = new uint64_t[POLDEG]();
	b_pk_q_shares = new uint64_t[POLDEG*L]();
	a_pk_q_shares = new uint64_t[POLDEG*L]();

	FILE *fp;

	fp = fopen("sage_generated_keys/sk", "r");
	for(int i=0; i<POLDEG; i++)
		fscanf(fp, "%lu", &s[i]);		
	fclose(fp);

	bitreverse_sequential(s);
	dit_fwd_ntt_qi(context, s, 0);


	fp = fopen("sage_generated_keys/a_pk", "r");
	for(int j=0; j<L; j++){
		for(int i=0; i<POLDEG; i++)
			fscanf(fp, "%lu", &a_pk_q_shares[j*POLDEG+i]);

		bitreverse_sequential(a_pk_q_shares+j*POLDEG);
		dit_fwd_ntt_qi(context, a_pk_q_shares+j*POLDEG, j);
	}

	fclose(fp);
	
	fp = fopen("sage_generated_keys/b_pk", "r");
	for(int j=0; j<L; j++){
		for(int i=0; i<POLDEG; i++)
			fscanf(fp, "%lu", &b_pk_q_shares[j*POLDEG+i]);		

		bitreverse_sequential(b_pk_q_shares+j*POLDEG);
		dit_fwd_ntt_qi(context, b_pk_q_shares+j*POLDEG, j);
	}


	fclose(fp);

}

Key::~Key(){
	delete[] s;
	delete[] b_pk_q_shares;
	delete[] a_pk_q_shares; 
}



class HBM:polynomial_functions{
	
	public:

	uint64_t *b_evk_q_shares[Num_Chiplets];	 
	uint64_t *a_evk_q_shares[Num_Chiplets];

	HBM(Context& context, Key& key, unsigned seed);			// constructor for creating pk-sk pair
	~HBM();
};
HBM::HBM(Context& context, Key& key, unsigned seed) {
	// HERE WE CAN SEE HOW DATA IS DISTRIBUTED ACROSS HBM
	for(int i=0;i<Num_Chiplets;i++){
		b_evk_q_shares[i] = new uint64_t[L*(Seq_Proc)*POLDEG]();
		a_evk_q_shares[i] = new uint64_t[L*(Seq_Proc)*POLDEG]();
	}
	FILE *fp;
	uint64_t *op;
	
	fp = fopen("sage_generated_keys/a_evk", "r");
	for(int i=0; i<L; i++){
		for(int chiplet_ctr=0; chiplet_ctr<Num_Chiplets;chiplet_ctr++){
			for(int j=0; j<Seq_Proc; j++){
				op = a_evk_q_shares[chiplet_ctr] + i*POLDEG*(Seq_Proc) + j*POLDEG;
				int t= j+chiplet_ctr*Seq_Proc;
				if(t<=L){
				for(int k=0; k<POLDEG; k++)
					fscanf(fp, "%lu", &op[k]);
				bitreverse_sequential(op);
				if(t==L)
					dit_fwd_ntt_pi(context, op, 0);
				else
					dit_fwd_ntt_qi(context, op, t);
			}
			}
		}
	}
	fclose(fp);

	fp = fopen("sage_generated_keys/b_evk", "r");
	for(int i=0; i<L; i++){
		for(int chiplet_ctr=0; chiplet_ctr<Num_Chiplets;chiplet_ctr++){
			for(int j=0; j<Seq_Proc; j++){
				op = b_evk_q_shares[chiplet_ctr] + i*POLDEG*(Seq_Proc) + j*POLDEG;
			for(int k=0; k<POLDEG; k++)
				fscanf(fp, "%lu", &op[k]);
			int t= j+chiplet_ctr*Seq_Proc;
			bitreverse_sequential(op);
			if(t==L)
				dit_fwd_ntt_pi(context, op, 0);
			else
				dit_fwd_ntt_qi(context, op, t);
			}
		}
	}
	fclose(fp);

	cout << "Keys loaded from file" << endl;

}

HBM::~HBM(){
	for(int i=0;i<Num_Chiplets;i++){
			delete[] b_evk_q_shares[i];
			delete[] a_evk_q_shares[i];
	}
}

class Plaintext{

	public:
	uint64_t *m;

	Plaintext();
	~Plaintext();
};

Plaintext::Plaintext(){
	m = new uint64_t[POLDEG*L]();
}

Plaintext::~Plaintext(){
	delete[] m;
}

// ciphertext = [c0_q_shares, c1_q_shares]
class Ciphertext:polynomial_functions{

	public:
	uint64_t level;
	uint64_t *c0_q_shares;	// Pointer to RNS polynomial vector c0_pk
	uint64_t *c1_q_shares;	// Pointer to RNS polynomial vector c1_pk

	Ciphertext(Context& context);			// constructor for ciphertext
	~Ciphertext();
};

Ciphertext::Ciphertext(Context& context){
	level = L-1;
	c0_q_shares = new uint64_t[POLDEG*L]();
	c1_q_shares = new uint64_t[POLDEG*L]();
}

Ciphertext::~Ciphertext(){
	delete[] c0_q_shares;
	delete[] c1_q_shares; 
}

// non_linearized_ciphertext = [d0_q_shares, d1_q_shares, d2_q_shares]
class non_linearized_Ciphertext:polynomial_functions{

	public:

	uint64_t *d0_q_shares;	
	uint64_t *d1_q_shares;	
	uint64_t *d2_q_shares;	

	non_linearized_Ciphertext(Context& context);			// constructor for ciphertext
	~non_linearized_Ciphertext();
};

non_linearized_Ciphertext::non_linearized_Ciphertext(Context& context){
	d0_q_shares = new uint64_t[POLDEG*L]();
	d1_q_shares = new uint64_t[POLDEG*L]();
	d2_q_shares = new uint64_t[POLDEG*L]();
}

non_linearized_Ciphertext::~non_linearized_Ciphertext(){
	delete[] d0_q_shares;
	delete[] d1_q_shares; 
	delete[] d2_q_shares; 
}



/*
// evk = [b_evk_q_shares, b_evk_p_shares, a_evk_q_shares, a_evk_p_shares]
class Evk:polynomial_functions{
	
	public:

	uint64_t *b_evk_q_shares;	 
	uint64_t *a_evk_q_shares;

	Evk(Context& context, Key& key, unsigned seed);			// constructor for creating pk-sk pair
	~Evk();
};

Evk::Evk(Context& context, Key& key, unsigned seed) {
	b_evk_q_shares = new uint64_t[L*(L+1)*POLDEG]();
	a_evk_q_shares = new uint64_t[L*(L+1)*POLDEG]();

	FILE *fp;
	uint64_t *op;

	fp = fopen("sage_generated_keys/a_evk", "r");
	for(int i=0; i<L; i++)
		for(int j=0; j<L+1; j++){
			op = a_evk_q_shares + i*POLDEG*(L+1) + j*POLDEG;
			for(int k=0; k<POLDEG; k++)
				fscanf(fp, "%lu", &op[k]);

			bitreverse_sequential(op);
			if(j==L)
				dit_fwd_ntt_pi(context, op, 0);
			else
				dit_fwd_ntt_qi(context, op, j);
		}
	fclose(fp);

	fp = fopen("sage_generated_keys/b_evk", "r");
	for(int i=0; i<L; i++)
		for(int j=0; j<L+1; j++){
			op = b_evk_q_shares + i*POLDEG*(L+1) + j*POLDEG;
			for(int k=0; k<POLDEG; k++)
				fscanf(fp, "%lu", &op[k]);

			bitreverse_sequential(op);
			if(j==L)
				dit_fwd_ntt_pi(context, op, 0);
			else
				dit_fwd_ntt_qi(context, op, j);
		}
	fclose(fp);

	cout << "Keys loaded from file" << endl;

}




Evk::~Evk(){
	delete[] b_evk_q_shares;
	delete[] a_evk_q_shares;
}
*/

