#define POLDEG 65536
#define LOG_POLDEG 16

#define N1 64
#define N2 POLDEG/N1
#define Num_Chiplets 4
#define Seq_Proc int((L+1)/Num_Chiplets)


#define M  2*POLDEG
#define L 7
#define K 1


long HW_time;


class Context{

	public:
	
	friend class polynomial_functions;

	uint64_t P_inv_mod_qi[L];
	uint64_t P_mod_qi[L];

	// New parameters
	uint64_t qVec[L] = {576460752340123649, 9007199282003969, 9007199284101121, 9007199318704129, 9007199338627073, 9007199343869953, 9007199350161409};
	uint64_t pVec[1] = {9007199379521537};
	uint64_t moduli[L+1] = {576460752340123649, 9007199282003969, 9007199284101121, 9007199318704129, 9007199338627073, 9007199343869953, 9007199350161409, 9007199379521537};
	
	uint64_t B_const[L][L+1] = {{1,0,0,0,0,0,0,6128844380391920},
	                            {0,1,0,0,0,0,0,3258278962458579},
								{0,0,1,0,0,0,0,1698135864756356},
								{0,0,0,1,0,0,0,7692416588738419},
								{0,0,0,0,1,0,0,7850798679818722},
								{0,0,0,0,0,1,0,2613029004197884},
								{0,0,0,0,0,0,1,1945665158767199}};
	
	uint64_t penc = 9007199254740992;	// Power of two for encode and decode.

	uint64_t ql_inv_mod_qj[L][L];

	uint64_t qVec_omega_initial[L][LOG_POLDEG], qVec_inv_omega_initial[L][LOG_POLDEG], qVec_primrt_table[L][LOG_POLDEG], qVec_inv_primrt_table[L][LOG_POLDEG];
	uint64_t pVec_omega_initial[K][LOG_POLDEG], pVec_inv_omega_initial[K][LOG_POLDEG], pVec_primrt_table[K][LOG_POLDEG], pVec_inv_primrt_table[K][LOG_POLDEG];

	uint64_t qVec_primrt_2Nth[L], qVec_inv_primrt_2Nth[L];
	uint64_t pVec_primrt_2Nth[K], pVec_inv_primrt_2Nth[K];

	Context();

};


Context::Context(){

	uint64_t qRoots[L],qRootsInv[L];
	uint64_t pRoots[K],pRootsInv[K];
	uint64_t temp;

	for (long i = 0; i < L; ++i){
		qRoots[i] = findMthRootOfUnity(2*M, qVec[i]);
		qRootsInv[i] = invMod(qRoots[i], qVec[i]);


	}

	for (long i = 0; i < L; ++i){
		mulMod(temp, qRoots[i], qRoots[i], qVec[i]);
		qVec_omega_initial[i][LOG_POLDEG-1] = temp/*qRoots[i]*/;	// M=2N th primitive root of unity

		qVec_inv_omega_initial[i][LOG_POLDEG-1] = invMod(qVec_omega_initial[i][LOG_POLDEG-1], qVec[i]);

		qVec_inv_primrt_2Nth[i] = invMod(qRoots[i], qVec[i]); 
		mulMod(qVec_inv_primrt_2Nth[i], qVec_inv_primrt_2Nth[i], qVec_inv_primrt_2Nth[i], qVec[i]);

		mulMod(temp, qRoots[i], qRoots[i], qVec[i]);
		mulMod(temp, temp, temp, qVec[i]);
		qVec_primrt_table[i][LOG_POLDEG-1] = temp;

		qVec_inv_primrt_table[i][LOG_POLDEG-1] = invMod(qVec_primrt_table[i][LOG_POLDEG-1], qVec[i]);
	}

	for (long i = 0; i < L; ++i){
		for(long j=LOG_POLDEG-2; j>=0; j--){
			mulMod(temp, qVec_primrt_table[i][j+1], qVec_primrt_table[i][j+1], qVec[i]);
			qVec_primrt_table[i][j] = temp;

			mulMod(temp, qVec_inv_primrt_table[i][j+1], qVec_inv_primrt_table[i][j+1], qVec[i]);
			qVec_inv_primrt_table[i][j] = temp;

			mulMod(temp, qVec_omega_initial[i][j+1], qVec_omega_initial[i][j+1], qVec[i]);
			qVec_omega_initial[i][j] = temp;

			mulMod(temp, qVec_inv_omega_initial[i][j+1], qVec_inv_omega_initial[i][j+1], qVec[i]);
			qVec_inv_omega_initial[i][j] = temp;
		}
	}



	for (long i=L-1; i>0; i--){
		for(long j=0; j<i; j++){
			ql_inv_mod_qj[i][j] = invMod(qVec[i], qVec[j]);
		}
	}

	for (long i=0; i<L; i++){
			P_inv_mod_qi[i] = invMod(pVec[0], qVec[i]);
			P_mod_qi[i] = (i==0) ? pVec[0] : pVec[0] - qVec[i]; 
	}


	for (long i = 0; i < K; ++i) {
		pRoots[i] = findMthRootOfUnity(2*M, pVec[i]);
		pRootsInv[i] = invMod(pRoots[i], pVec[i]);
		

	}

	for (long i = 0; i < K; ++i){
		mulMod(temp, pRoots[i], pRoots[i], pVec[i]);
		pVec_omega_initial[i][LOG_POLDEG-1] = temp/*pRoots[i]*/;	// M=2N th primitive root of unity

		pVec_inv_omega_initial[i][LOG_POLDEG-1] = invMod(pVec_omega_initial[i][LOG_POLDEG-1], pVec[i]);

		pVec_inv_primrt_2Nth[i] = invMod(pRoots[i], pVec[i]); 
		mulMod(pVec_inv_primrt_2Nth[i], pVec_inv_primrt_2Nth[i], pVec_inv_primrt_2Nth[i], pVec[i]);

		mulMod(temp, pRoots[i], pRoots[i], pVec[i]);
		mulMod(temp, temp, temp, pVec[i]);
		pVec_primrt_table[i][LOG_POLDEG-1] = temp;

		pVec_inv_primrt_table[i][LOG_POLDEG-1] = invMod(pVec_primrt_table[i][LOG_POLDEG-1], pVec[i]);
	}

	for (long i = 0; i < K; ++i){
		for(long j=LOG_POLDEG-2; j>=0; j--){
			mulMod(temp, pVec_primrt_table[i][j+1], pVec_primrt_table[i][j+1], pVec[i]);
			pVec_primrt_table[i][j] = temp;

			mulMod(temp, pVec_inv_primrt_table[i][j+1], pVec_inv_primrt_table[i][j+1], pVec[i]);
			pVec_inv_primrt_table[i][j] = temp;

			mulMod(temp, pVec_omega_initial[i][j+1], pVec_omega_initial[i][j+1], pVec[i]);
			pVec_omega_initial[i][j] = temp;

			mulMod(temp, pVec_inv_omega_initial[i][j+1], pVec_inv_omega_initial[i][j+1], pVec[i]);
			pVec_inv_omega_initial[i][j] = temp;
		}
	}

	
}

