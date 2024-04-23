
class polynomial_functions{

	public:
	
	void bitreverse_sequential(uint64_t *a);	

	void dit_fwd_ntt_qi(Context&  context, uint64_t* a, long index);

	void dit_fwd_ntt_pi(Context&  context, uint64_t* a, long index);

	void dif_inv_ntt_qi(Context&  context, uint64_t *a, long index);

	void dif_inv_ntt_pi(Context&  context, uint64_t *a, long index);

	void polynomial_addition_qi(Context&  context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index);

	void polynomial_addition_pi(Context&  context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index);

	void polynomial_subtraction_qi(Context&  context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index);

	void polynomial_subtraction_pi(Context&  context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index);

	void polynomial_multiplication_qi(Context&  context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index);

	void polynomial_multiplication_pi(Context&  context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index);

	void diadic_multiplication_qi(Context&  context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index);

	void diadic_multiplication_pi(Context&  context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index);

	uint64_t tiny_hash(uint64_t *op1);

	friend class Context;
	
};

void polynomial_functions::bitreverse_sequential(uint64_t* a)
{
	int i;
	int bit1, bit2, bit3, bit4, bit5, bit6, bit7, bit8, bit9, bit10, bit11, bit12, bit13, bit14, bit15, bit16;
	int	swp_index;
	uint64_t temp;

	for(i=0; i<65536; i++)
	{
		bit1 = i%2;
		bit2 = (i>>1)%2;
		bit3 = (i>>2)%2;
		bit4 = (i>>3)%2;
		bit5 = (i>>4)%2;
		bit6 = (i>>5)%2;
		bit7 = (i>>6)%2;
		bit8 = (i>>7)%2;
		bit9 = (i>>8)%2;
		bit10 = (i>>9)%2;
		bit11 = (i>>10)%2;
		bit12 = (i>>11)%2;
		bit13 = (i>>12)%2;
		bit14 = (i>>13)%2;
		bit15 = (i>>14)%2;
		bit16 = (i>>15)%2;

		swp_index = bit1*32768 + bit2*16384 + bit3*8192 + bit4*4096 + bit5*2048 + bit6*1024 + bit7*512 + bit8*256 + bit9*128 + bit10*64 + bit11*32 + bit12*16 + bit13*8 + bit14*4 + bit15*2 + bit16;      
	   
		if(swp_index>i)
		{
			temp = a[i];
			a[i] = a[swp_index];
			a[swp_index] = temp;
		}      
	}  
}


void polynomial_functions::dit_fwd_ntt_qi(Context&  context, uint64_t* a, long index)
{
	int i, j, k, m;
	uint64_t u1, t1, u2, t2;
	uint64_t primrt, omega, temp;

	uint64_t modulus;
	if (index<L)
		modulus = context.qVec[index];
	else
		modulus = context.pVec[index-L];

	if(index<L){
	for(m=2; m<=POLDEG; m=m*2)
	{
		if(m==2) primrt = context.qVec_primrt_table[index][0];
		if(m==4) primrt = context.qVec_primrt_table[index][1];
		if(m==8) primrt = context.qVec_primrt_table[index][2];
		if(m==16) primrt = context.qVec_primrt_table[index][3];
		if(m==32) primrt = context.qVec_primrt_table[index][4];
		if(m==64) primrt = context.qVec_primrt_table[index][5];
		if(m==128) primrt = context.qVec_primrt_table[index][6];
		if(m==256) primrt = context.qVec_primrt_table[index][7];
		if(m==512) primrt = context.qVec_primrt_table[index][8];
		if(m==1024) primrt = context.qVec_primrt_table[index][9];
		if(m==2048) primrt = context.qVec_primrt_table[index][10];
		if(m==4096) primrt = context.qVec_primrt_table[index][11];
		if(m==8192) primrt = context.qVec_primrt_table[index][12];
		if(m==16384) primrt = context.qVec_primrt_table[index][13];
		if(m==32768) primrt = context.qVec_primrt_table[index][14];
		if(m==65536) primrt = context.qVec_primrt_table[index][15];

		if(m==2) omega = context.qVec_omega_initial[index][0];
		if(m==4) omega = context.qVec_omega_initial[index][1];
		if(m==8) omega = context.qVec_omega_initial[index][2];
		if(m==16) omega = context.qVec_omega_initial[index][3];
		if(m==32) omega = context.qVec_omega_initial[index][4];
		if(m==64) omega = context.qVec_omega_initial[index][5];
		if(m==128) omega = context.qVec_omega_initial[index][6];
		if(m==256) omega = context.qVec_omega_initial[index][7];
		if(m==512) omega = context.qVec_omega_initial[index][8];
		if(m==1024) omega = context.qVec_omega_initial[index][9];
		if(m==2048) omega = context.qVec_omega_initial[index][10];
		if(m==4096) omega = context.qVec_omega_initial[index][11];
		if(m==8192) omega = context.qVec_omega_initial[index][12];
		if(m==16384) omega = context.qVec_omega_initial[index][13];
		if(m==32768) omega = context.qVec_omega_initial[index][14];
		if(m==65536) omega = context.qVec_omega_initial[index][15];

		for(j=0; j<m/2; j++)
		{
			for(k=0; k<POLDEG; k=k+m)
			{
				t1 = a[k+j+m/2];
				u1 = a[k+j];

				mulMod(temp, omega, t1, modulus);
				t1 = temp;

				addMod(temp,u1,t1, modulus);
				a[k+j] = temp;

				subMod(temp,u1,t1, modulus);
				a[k+j+m/2] = temp;
			}      
			mulMod(omega, omega, primrt, modulus);
		}
	}
	}
	else{
	index=index-L;
	for(m=2; m<=POLDEG; m=m*2)
	{
		if(m==2) primrt = context.pVec_primrt_table[index][0];
		if(m==4) primrt = context.pVec_primrt_table[index][1];
		if(m==8) primrt = context.pVec_primrt_table[index][2];
		if(m==16) primrt = context.pVec_primrt_table[index][3];
		if(m==32) primrt = context.pVec_primrt_table[index][4];
		if(m==64) primrt = context.pVec_primrt_table[index][5];
		if(m==128) primrt = context.pVec_primrt_table[index][6];
		if(m==256) primrt = context.pVec_primrt_table[index][7];
		if(m==512) primrt = context.pVec_primrt_table[index][8];
		if(m==1024) primrt = context.pVec_primrt_table[index][9];
		if(m==2048) primrt = context.pVec_primrt_table[index][10];
		if(m==4096) primrt = context.pVec_primrt_table[index][11];
		if(m==8192) primrt = context.pVec_primrt_table[index][12];
		if(m==16384) primrt = context.pVec_primrt_table[index][13];
		if(m==32768) primrt = context.pVec_primrt_table[index][14];
		if(m==65536) primrt = context.pVec_primrt_table[index][15];

		if(m==2) omega = context.pVec_omega_initial[index][0];
		if(m==4) omega = context.pVec_omega_initial[index][1];
		if(m==8) omega = context.pVec_omega_initial[index][2];
		if(m==16) omega = context.pVec_omega_initial[index][3];
		if(m==32) omega = context.pVec_omega_initial[index][4];
		if(m==64) omega = context.pVec_omega_initial[index][5];
		if(m==128) omega = context.pVec_omega_initial[index][6];
		if(m==256) omega = context.pVec_omega_initial[index][7];
		if(m==512) omega = context.pVec_omega_initial[index][8];
		if(m==1024) omega = context.pVec_omega_initial[index][9];
		if(m==2048) omega = context.pVec_omega_initial[index][10];
		if(m==4096) omega = context.pVec_omega_initial[index][11];
		if(m==8192) omega = context.pVec_omega_initial[index][12];
		if(m==16384) omega = context.pVec_omega_initial[index][13];
		if(m==32768) omega = context.pVec_omega_initial[index][14];
		if(m==65536) omega = context.pVec_omega_initial[index][15];

		for(j=0; j<m/2; j++)
		{
			for(k=0; k<POLDEG; k=k+m)
			{
				t1 = a[k+j+m/2];
				u1 = a[k+j];

				mulMod(temp, omega, t1, modulus);
				t1 = temp;

				addMod(temp,u1,t1, modulus);
				a[k+j] = temp;

				subMod(temp,u1,t1, modulus);
				a[k+j+m/2] = temp;
			}      
			mulMod(omega, omega, primrt, modulus);
		}
	}
	
	}
}

void polynomial_functions::dit_fwd_ntt_pi(Context&  context, uint64_t* a, long index)
{
	int i, j, k, m;
	uint64_t u1, t1, u2, t2;
	uint64_t primrt, omega, temp;

	uint64_t modulus = context.pVec[index];

	for(m=2; m<=POLDEG; m=m*2)
	{
		if(m==2) primrt = context.pVec_primrt_table[index][0];
		if(m==4) primrt = context.pVec_primrt_table[index][1];
		if(m==8) primrt = context.pVec_primrt_table[index][2];
		if(m==16) primrt = context.pVec_primrt_table[index][3];
		if(m==32) primrt = context.pVec_primrt_table[index][4];
		if(m==64) primrt = context.pVec_primrt_table[index][5];
		if(m==128) primrt = context.pVec_primrt_table[index][6];
		if(m==256) primrt = context.pVec_primrt_table[index][7];
		if(m==512) primrt = context.pVec_primrt_table[index][8];
		if(m==1024) primrt = context.pVec_primrt_table[index][9];
		if(m==2048) primrt = context.pVec_primrt_table[index][10];
		if(m==4096) primrt = context.pVec_primrt_table[index][11];
		if(m==8192) primrt = context.pVec_primrt_table[index][12];
		if(m==16384) primrt = context.pVec_primrt_table[index][13];
		if(m==32768) primrt = context.pVec_primrt_table[index][14];
		if(m==65536) primrt = context.pVec_primrt_table[index][15];

		if(m==2) omega = context.pVec_omega_initial[index][0];
		if(m==4) omega = context.pVec_omega_initial[index][1];
		if(m==8) omega = context.pVec_omega_initial[index][2];
		if(m==16) omega = context.pVec_omega_initial[index][3];
		if(m==32) omega = context.pVec_omega_initial[index][4];
		if(m==64) omega = context.pVec_omega_initial[index][5];
		if(m==128) omega = context.pVec_omega_initial[index][6];
		if(m==256) omega = context.pVec_omega_initial[index][7];
		if(m==512) omega = context.pVec_omega_initial[index][8];
		if(m==1024) omega = context.pVec_omega_initial[index][9];
		if(m==2048) omega = context.pVec_omega_initial[index][10];
		if(m==4096) omega = context.pVec_omega_initial[index][11];
		if(m==8192) omega = context.pVec_omega_initial[index][12];
		if(m==16384) omega = context.pVec_omega_initial[index][13];
		if(m==32768) omega = context.pVec_omega_initial[index][14];
		if(m==65536) omega = context.pVec_omega_initial[index][15];

		for(j=0; j<m/2; j++)
		{
			for(k=0; k<POLDEG; k=k+m)
			{
				t1 = a[k+j+m/2];
				u1 = a[k+j];

				mulMod(temp, omega, t1, context.pVec[index]);
				t1 = temp;

				addMod(temp,u1,t1, context.pVec[index]);
				a[k+j] = temp;

				subMod(temp,u1,t1, context.pVec[index]);
				a[k+j+m/2] = temp;
			}      
			mulMod(omega, omega, primrt, context.pVec[index]);
		}
	}
}


void polynomial_functions::dif_inv_ntt_qi(Context&  context, uint64_t *a, long index)
{
	int i, j, k, m;
	uint64_t u1, t1, u2, t2;
	uint64_t primrt, omega, temp;

	uint64_t modulus = context.qVec[index];

	uint64_t inv_2 = invMod(2, modulus);

	for(m=POLDEG; m>=2; m=m/2)
	{
		if(m==2) primrt = context.qVec_inv_primrt_table[index][0];
		if(m==4) primrt = context.qVec_inv_primrt_table[index][1];
		if(m==8) primrt = context.qVec_inv_primrt_table[index][2];
		if(m==16) primrt = context.qVec_inv_primrt_table[index][3];
		if(m==32) primrt = context.qVec_inv_primrt_table[index][4];
		if(m==64) primrt = context.qVec_inv_primrt_table[index][5];
		if(m==128) primrt = context.qVec_inv_primrt_table[index][6];
		if(m==256) primrt = context.qVec_inv_primrt_table[index][7];
		if(m==512) primrt = context.qVec_inv_primrt_table[index][8];
		if(m==1024) primrt = context.qVec_inv_primrt_table[index][9];
		if(m==2048) primrt = context.qVec_inv_primrt_table[index][10];
		if(m==4096) primrt = context.qVec_inv_primrt_table[index][11];
		if(m==8192) primrt = context.qVec_inv_primrt_table[index][12];
		if(m==16384) primrt = context.qVec_inv_primrt_table[index][13];
		if(m==32768) primrt = context.qVec_inv_primrt_table[index][14];
		if(m==65536) primrt = context.qVec_inv_primrt_table[index][15];

		//omega = 1;
		/*
		if(m==2)     mulMod(omega, context.qVec_inv_primrt_table[index][1] , inv_2, modulus);
		if(m==4)     mulMod(omega, context.qVec_inv_primrt_table[index][2] , inv_2, modulus);
		if(m==8)     mulMod(omega, context.qVec_inv_primrt_table[index][3] , inv_2, modulus);
		if(m==16)    mulMod(omega, context.qVec_inv_primrt_table[index][4] , inv_2, modulus);
		if(m==32)    mulMod(omega, context.qVec_inv_primrt_table[index][5] , inv_2, modulus);
		if(m==64)    mulMod(omega, context.qVec_inv_primrt_table[index][6] , inv_2, modulus);
		if(m==128)   mulMod(omega, context.qVec_inv_primrt_table[index][7] , inv_2, modulus);
		if(m==256)   mulMod(omega, context.qVec_inv_primrt_table[index][8] , inv_2, modulus);
		if(m==512)   mulMod(omega, context.qVec_inv_primrt_table[index][9] , inv_2, modulus);
		if(m==1024)  mulMod(omega, context.qVec_inv_primrt_table[index][10], inv_2, modulus);
		if(m==2048)  mulMod(omega, context.qVec_inv_primrt_table[index][11], inv_2, modulus);
		if(m==4096)  mulMod(omega, context.qVec_inv_primrt_table[index][12], inv_2, modulus);
		if(m==8192)  mulMod(omega, context.qVec_inv_primrt_table[index][13], inv_2, modulus);
		if(m==16384) mulMod(omega, context.qVec_inv_primrt_15th[index], inv_2, modulus);
		*/
		if(m==2)     mulMod(omega, context.qVec_inv_omega_initial[index][0] , inv_2, modulus);
		if(m==4)     mulMod(omega, context.qVec_inv_omega_initial[index][1] , inv_2, modulus);
		if(m==8)     mulMod(omega, context.qVec_inv_omega_initial[index][2] , inv_2, modulus);
		if(m==16)    mulMod(omega, context.qVec_inv_omega_initial[index][3] , inv_2, modulus);
		if(m==32)    mulMod(omega, context.qVec_inv_omega_initial[index][4] , inv_2, modulus);
		if(m==64)    mulMod(omega, context.qVec_inv_omega_initial[index][5] , inv_2, modulus);
		if(m==128)   mulMod(omega, context.qVec_inv_omega_initial[index][6] , inv_2, modulus);
		if(m==256)   mulMod(omega, context.qVec_inv_omega_initial[index][7] , inv_2, modulus);
		if(m==512)   mulMod(omega, context.qVec_inv_omega_initial[index][8] , inv_2, modulus);
		if(m==1024)  mulMod(omega, context.qVec_inv_omega_initial[index][9] , inv_2, modulus);
		if(m==2048)  mulMod(omega, context.qVec_inv_omega_initial[index][10], inv_2, modulus);
		if(m==4096)  mulMod(omega, context.qVec_inv_omega_initial[index][11], inv_2, modulus);
		if(m==8192)  mulMod(omega, context.qVec_inv_omega_initial[index][12], inv_2, modulus);
		if(m==16384) mulMod(omega, context.qVec_inv_omega_initial[index][13], inv_2, modulus);
		if(m==32768) mulMod(omega, context.qVec_inv_omega_initial[index][14], inv_2, modulus);
		if(m==65536) mulMod(omega, context.qVec_inv_omega_initial[index][15], inv_2, modulus);

		for(j=0; j<m/2; j++)
		{
			for(k=0; k<POLDEG; k=k+m)
			{
				t1 = a[k+j+m/2];
				u1 = a[k+j];

				addMod(temp,u1,t1,context.qVec[index]);
				//a[k+j]=temp;
				a[k+j]= (temp >> 1) + (temp & 0x1)*((context.qVec[index]+1)>>1);

				subMod(temp,u1,t1,context.qVec[index]);
				t1=temp;

				mulMod(temp, omega, t1, context.qVec[index]);
				a[k+j+m/2]=temp;
			}      
			mulMod(omega, omega, primrt, context.qVec[index]);
		}
	}

	
	// Start: Generate ROM for post-scaling in bit-reverse order
	/*
	omega = 1;
	m = 16384;
	uint64_t omega_arr[16384];
	uint64_t inv_m = invMod(m, context.qVec[index]);
	mulMod(temp, omega, inv_m, context.qVec[index]);	// omega = inv_primrt_15th * (n^-1)

	omega_arr[0] = temp;
	for(j=1; j<m; j++)
	{
		mulMod(temp, omega_arr[j-1], context.qVec_inv_primrt_15th[index], context.qVec[index]);
		omega_arr[j] = temp;
	}
	bitreverse_sequential(omega_arr);
	// End: Generate ROM for post-scaling in bit-reverse order


	for(j=0; j<m; j++)
	{
		mulMod(temp, a[j], omega_arr[j], context.qVec[index]);
		a[j] = temp;
	}
	*/
}

void polynomial_functions::dif_inv_ntt_pi(Context&  context, uint64_t *a, long index)
{
	int i, j, k, m;
	uint64_t u1, t1, u2, t2;
	uint64_t primrt, omega, temp;

	uint64_t modulus = context.pVec[index];

	uint64_t inv_2 = invMod(2, modulus);

	for(m=POLDEG; m>1; m=m/2)
	{
		if(m==2) primrt = context.pVec_inv_primrt_table[index][0];
		if(m==4) primrt = context.pVec_inv_primrt_table[index][1];
		if(m==8) primrt = context.pVec_inv_primrt_table[index][2];
		if(m==16) primrt = context.pVec_inv_primrt_table[index][3];
		if(m==32) primrt = context.pVec_inv_primrt_table[index][4];
		if(m==64) primrt = context.pVec_inv_primrt_table[index][5];
		if(m==128) primrt = context.pVec_inv_primrt_table[index][6];
		if(m==256) primrt = context.pVec_inv_primrt_table[index][7];
		if(m==512) primrt = context.pVec_inv_primrt_table[index][8];
		if(m==1024) primrt = context.pVec_inv_primrt_table[index][9];
		if(m==2048) primrt = context.pVec_inv_primrt_table[index][10];
		if(m==4096) primrt = context.pVec_inv_primrt_table[index][11];
		if(m==8192) primrt = context.pVec_inv_primrt_table[index][12];
		if(m==16384) primrt = context.pVec_inv_primrt_table[index][13];
		if(m==32768) primrt = context.pVec_inv_primrt_table[index][14];
		if(m==65536) primrt = context.pVec_inv_primrt_table[index][15];

	    //omega = 1;
		/*
		if(m==2)     mulMod(omega, context.pVec_inv_primrt_table[index][1] , inv_2, modulus);
		if(m==4)     mulMod(omega, context.pVec_inv_primrt_table[index][2] , inv_2, modulus);
		if(m==8)     mulMod(omega, context.pVec_inv_primrt_table[index][3] , inv_2, modulus);
		if(m==16)    mulMod(omega, context.pVec_inv_primrt_table[index][4] , inv_2, modulus);
		if(m==32)    mulMod(omega, context.pVec_inv_primrt_table[index][5] , inv_2, modulus);
		if(m==64)    mulMod(omega, context.pVec_inv_primrt_table[index][6] , inv_2, modulus);
		if(m==128)   mulMod(omega, context.pVec_inv_primrt_table[index][7] , inv_2, modulus);
		if(m==256)   mulMod(omega, context.pVec_inv_primrt_table[index][8] , inv_2, modulus);
		if(m==512)   mulMod(omega, context.pVec_inv_primrt_table[index][9] , inv_2, modulus);
		if(m==1024)  mulMod(omega, context.pVec_inv_primrt_table[index][10], inv_2, modulus);
		if(m==2048)  mulMod(omega, context.pVec_inv_primrt_table[index][11], inv_2, modulus);
		if(m==4096)  mulMod(omega, context.pVec_inv_primrt_table[index][12], inv_2, modulus);
		if(m==8192)  mulMod(omega, context.pVec_inv_primrt_table[index][13], inv_2, modulus);
		if(m==16384) mulMod(omega, context.pVec_inv_primrt_15th[index], inv_2, modulus);
		*/
		if(m==2)     mulMod(omega, context.pVec_inv_omega_initial[index][0] , inv_2, modulus);
		if(m==4)     mulMod(omega, context.pVec_inv_omega_initial[index][1] , inv_2, modulus);
		if(m==8)     mulMod(omega, context.pVec_inv_omega_initial[index][2] , inv_2, modulus);
		if(m==16)    mulMod(omega, context.pVec_inv_omega_initial[index][3] , inv_2, modulus);
		if(m==32)    mulMod(omega, context.pVec_inv_omega_initial[index][4] , inv_2, modulus);
		if(m==64)    mulMod(omega, context.pVec_inv_omega_initial[index][5] , inv_2, modulus);
		if(m==128)   mulMod(omega, context.pVec_inv_omega_initial[index][6] , inv_2, modulus);
		if(m==256)   mulMod(omega, context.pVec_inv_omega_initial[index][7] , inv_2, modulus);
		if(m==512)   mulMod(omega, context.pVec_inv_omega_initial[index][8] , inv_2, modulus);
		if(m==1024)  mulMod(omega, context.pVec_inv_omega_initial[index][9] , inv_2, modulus);
		if(m==2048)  mulMod(omega, context.pVec_inv_omega_initial[index][10], inv_2, modulus);
		if(m==4096)  mulMod(omega, context.pVec_inv_omega_initial[index][11], inv_2, modulus);
		if(m==8192)  mulMod(omega, context.pVec_inv_omega_initial[index][12], inv_2, modulus);
		if(m==16384) mulMod(omega, context.pVec_inv_omega_initial[index][13], inv_2, modulus);
		if(m==32768) mulMod(omega, context.pVec_inv_omega_initial[index][14], inv_2, modulus);
		if(m==65536) mulMod(omega, context.pVec_inv_omega_initial[index][15], inv_2, modulus);

		for(j=0; j<m/2; j++)
		{
			for(k=0; k<POLDEG; k=k+m)
			{
				t1 = a[k+j+m/2];
				u1 = a[k+j];

				addMod(temp,u1,t1,context.pVec[index]);
				//a[k+j]=temp;
				a[k+j]= (temp >> 1) + (temp & 0x1)*((context.pVec[index]+1)>>1);

				subMod(temp,u1,t1,context.pVec[index]);
				t1=temp;

				mulMod(temp, omega, t1, context.pVec[index]);
				a[k+j+m/2]=temp;

			}      
			mulMod(omega, omega, primrt, context.pVec[index]);
		}
	}


	// Start: Generate ROM for post-scaling in bit-reverse order
	/*
	omega = 1;
	m = 16384;
	uint64_t omega_arr[16384];
	uint64_t inv_m = invMod(m, context.pVec[index]);
	mulMod(temp, omega, inv_m, context.pVec[index]);	// omega = inv_primrt_15th * (n^-1)

	omega_arr[0] = temp;
	for(j=1; j<m; j++)
	{
		mulMod(temp, omega_arr[j-1], context.pVec_inv_primrt_15th[index], context.pVec[index]);
		omega_arr[j] = temp;
	}
	bitreverse_sequential(omega_arr);
	// End: Generate ROM for post-scaling in bit-reverse order

	for(j=0; j<m; j++)
	{
		mulMod(temp, a[j], omega_arr[j], context.pVec[index]);
		a[j] = temp;
	}
	*/
}


void polynomial_functions::polynomial_addition_qi(Context& context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index)
{
	int j;

	if (index<L){
	for(j=0; j<POLDEG; j++)
		addMod(result[j], op1[j], op2[j], context.qVec[index]); 
		}
	else{
	for(j=0; j<POLDEG; j++)
		addMod(result[j], op1[j], op2[j], context.pVec[index-L]);
	}	
}

void polynomial_functions::polynomial_subtraction_qi(Context& context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index)
{
	int j;
	
	if (index<L){
	for(j=0; j<POLDEG; j++)
		subMod(result[j], op1[j], op2[j], context.qVec[index]); 
		}
	else{
	for(j=0; j<POLDEG; j++)
		subMod(result[j], op1[j], op2[j], context.pVec[index-L]);
	}	
}

void polynomial_functions::polynomial_addition_pi(Context& context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index)
{
	int j;

	for(j=0; j<POLDEG; j++)
		addMod(result[j], op1[j], op2[j], context.pVec[index]); 	
}

void polynomial_functions::polynomial_subtraction_pi(Context& context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index)
{
	int j;

	for(j=0; j<POLDEG; j++)
		subMod(result[j], op1[j], op2[j], context.pVec[index]); 	
}

void polynomial_functions::polynomial_multiplication_qi(Context& context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index)
{
	int j;
	uint64_t temp;

	uint64_t* temp1 = new uint64_t[POLDEG]();
	copy(op1, op1 + POLDEG, temp1);

	uint64_t* temp2 = new uint64_t[POLDEG]();
	copy(op2, op2 + POLDEG, temp2);

	// It is important to keep input in bit-reverse for correct multiplication, although DiT-DiF doesn't require bitreverse inbetween.
	bitreverse_sequential(temp1);	
	bitreverse_sequential(temp2);

	dit_fwd_ntt_qi(context, temp1, index);
	dit_fwd_ntt_qi(context, temp2, index);

	// Coefficient wise multiply
	for(j=0; j<POLDEG; j++) {
		mulMod(result[j], temp1[j], temp2[j], context.qVec[index]); 
	}

	dif_inv_ntt_qi(context, result, index);

	bitreverse_sequential(result);

	delete[] temp1;
	delete[] temp2;
}

void polynomial_functions::diadic_multiplication_qi(Context& context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index)
{
	// Coefficient wise multiply
	if (index<L){
	for(int j=0; j<POLDEG; j++)
		mulMod(result[j], op1[j], op2[j], context.qVec[index]); 
		}
	else{
	for(int j=0; j<POLDEG; j++)
		mulMod(result[j], op1[j], op2[j], context.pVec[index-L]);
	}
}

void polynomial_functions::diadic_multiplication_pi(Context& context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index)
{
	// Coefficient wise multiply
	for(long j=0; j<POLDEG; j++) {
		mulMod(result[j], op1[j], op2[j], context.pVec[index]); 
	}
}


void polynomial_functions::polynomial_multiplication_pi(Context& context, uint64_t *result, uint64_t *op1, uint64_t *op2, long index)
{
	int j;

	uint64_t* temp1 = new uint64_t[POLDEG]();
	copy(op1, op1 +  POLDEG, temp1);

	uint64_t* temp2 = new uint64_t[POLDEG]();
	copy(op2, op2 + POLDEG, temp2);

	// It is important to keep input in bit-reverse for correct multiplication, although DiT-DiF doesn't require bitreverse inbetween.
	bitreverse_sequential(temp1);	
	bitreverse_sequential(temp2);

	dit_fwd_ntt_pi(context, temp1, index);
	dit_fwd_ntt_pi(context, temp2, index);

	// Coefficient wise multiply
	for(j=0; j<POLDEG; j++)
		mulMod(result[j], temp1[j], temp2[j], context.pVec[index]); 	

	dif_inv_ntt_pi(context, result, index);

	bitreverse_sequential(result);

	delete[] temp1;
	delete[] temp2;
}


