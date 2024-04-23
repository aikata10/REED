
from parameters import *
from random import randint

# 
# -- Proposed Hybrid NTT flow:
#
# MINTT: UNROLLNTT(N-->N) --> HP(N-->N)      --> PIPENTT(N-->R) --> PP(R-->R)        --> Write as N
# MFNTT: PP(N-->N)        --> PIPENTT(N-->R) --> HP(R-->R)      --> UNROLLNTT(R-->R) --> Write as N
#
# -- Definitions:
#
# FNTT     : Forward Number Theoretic Transform
# INTT     : Inverse Number Theoretic Transform
# MFNTT    : Merged FNTT (FNTT with negative-wrapped convolution)
# MINTT    : Merged INTT (INTT with negative-wrapped convolution)
# UNROLLNTT: FNTT with fully-unrolled design approach
# PIPENTT  : FNTT with pipelined (SDF) design approach
# HP       : Hadamard Product
# PP       : Pre-/Post-processing
# N        : Normal-ordered polynomial
# R        : Bit-reverse-ordered polynomial
# DIT      : Decimation-in-time
# DIF      : Decimation-in-frequency
#

# -------------------------------- 0. MFNTT reference
A        = [i for i in range(n)]                 # Input A (in normal order)
A_fntt_r = MFNTT_DIT_NR(A,psi_table,q)           # A_fntt_r (in bitreversed order)
A_fntt_n = arr_reverse(A_fntt_r)                 # A_fntt_n (in normal order)
A_intt_n_matrix = arr_array2matrix(A_fntt_n,k,m) # A_intt_n_matrix (in matrix format) 

# Print reference data
if USE_TOY_EXAMPLE: print("--->>> Reference computation:")
if USE_TOY_EXAMPLE: print("A:")
if USE_TOY_EXAMPLE: print(A)
if USE_TOY_EXAMPLE: print()
if USE_TOY_EXAMPLE: print("MFNTT(A):")
if USE_TOY_EXAMPLE: print(A_fntt_n)
if USE_TOY_EXAMPLE: print()
if USE_TOY_EXAMPLE: print("Input MFNTT(A) in matrix form:")
if USE_TOY_EXAMPLE: print_matrix_inmem_row(A_intt_n_matrix,print_dim,print_mat,True)

# -------------------------------- 1. Apply our Hybrid MINTT
print("-"*119 + " Hybrid MINTT")
print()

# 1.1. Re-order data for using NTT circuit for INTT (A_fntt_n (in normal order) has m columns of k-sized vector)
#      This can be done reading data in proper order (without extra hw/cycle cost)
A_intt_n_rearranged_matrix = arr_array2matrix(arr_rearrange(A_fntt_n),k,m)
print("--->>> Input in matrix form (re-ordered - this can be done reading data from memory in proper order in HW):")
print_matrix_inmem_row(A_intt_n_rearranged_matrix,print_dim,print_mat)

# 1.2. UNROLLNTT(N-->N)
#      Apply NTT to m columns of k-sized vector
print("--->>> Applying {}-pt FNTTs to {} columns using unrolled NTT approach (sequentially):".format(k,m))
for i in range(m):
    A_intt_n_rearranged_matrix[i] = arr_reverse(FNTT_DIF_NR(A_intt_n_rearranged_matrix[i],w_k_table,q))
    print_matrix_inmem_row(A_intt_n_rearranged_matrix,print_dim,print_mat,False,True,i,"FNTT")

# 1.3. HP(N-->N)
#      Hadamard Product
print("--->>> Performing Hadamard Product (sequentially):")
for i in range(m):
    for j in range(k):
        A_intt_n_rearranged_matrix[i][j] = (A_intt_n_rearranged_matrix[i][j]*W_H_intt[i][j])%q
    print_matrix_inmem_row(A_intt_n_rearranged_matrix,print_dim,print_mat,False,True,i," HP ")

# 1.4. !!! (Here we use transpose just to mimic data movement - this happens naturally in HW) !!!
A_intt_n_rearranged_matrix = [list(i) for i in zip(*A_intt_n_rearranged_matrix)]

# 1.5. PIPENTT(N-->R)
#      Apply NTT to k columns of m-sized vector
print("--->>> Applying {}-pt FNTTs to {} rows using pipelined NTT approach (in parallel):".format(m,k))
for i in range(k):
    A_intt_n_rearranged_matrix[i] = FNTT_DIF_NR(A_intt_n_rearranged_matrix[i],w_m_table,q)
print_matrix_inmem_col(A_intt_n_rearranged_matrix,print_dim,print_mat,False,True,-1,"FNTT")

# 1.6. PP(R-->R)
#      Post-processing
print("--->>> Performing Post-processing (sequentially):")
for i in range(m):
    for j in range(k):
        A_intt_n_rearranged_matrix[j][i] = (A_intt_n_rearranged_matrix[j][i]*PP_intt[j][i])%q
    print_matrix_inmem_col(A_intt_n_rearranged_matrix,print_dim,print_mat,False,True,i," PP ")

# 1.7. Write as N
#      This can be done writing data in proper order (without extra hw/cycle cost)
print("--->>> Output in matrix form (re-ordered - this can be done writing data back to memory in proper order in HW):")
for i in range(k):
    A_intt_n_rearranged_matrix[i] = arr_reverse(A_intt_n_rearranged_matrix[i])
print_matrix_inmem_col(A_intt_n_rearranged_matrix,print_dim,print_mat,True)

# 1.8. Convert from matrix form to polynomial form (in normal order)
A_intt_n_rearranged = arr_matrix2array(A_intt_n_rearranged_matrix)

# Sanity check
if sum([abs(p0-p1) for p0,p1 in zip(A,A_intt_n_rearranged)]) == 0:
    print("MINTT TEST: OK")
else:
    print("MINTT TEST: FAIL")

print()
print("-"*119 + " Hybrid MFNTT")
print()

# -------------------------------- 2. Apply our Hybrid MFNTT
# 2.1. A_fntt_n_matrix has k columns of m-sized vector
print("--->>> Input in matrix form:")
A_fntt_n_matrix = A_intt_n_rearranged_matrix
print_matrix_inmem_col(A_fntt_n_matrix,print_dim,print_mat, True)

# 2.2. PP(N-->N)
#      Pre-processing   
print("--->>> Performing Pre-processing (sequentially):")     
for i in range(m):
    for j in range(k):
        A_fntt_n_matrix[j][i] = (A_fntt_n_matrix[j][i]*PP_fntt[j][i])%q
    print_matrix_inmem_col(A_fntt_n_matrix,print_dim,print_mat,False,True,i," PP ")

# 2.3. PIPENTT(N-->R) 
#      Apply NTT to k columns of m-sized vector
print("--->>> Applying {}-pt FNTTs to {} rows using pipelined NTT approach (in parallel):".format(m,k))
for i in range(k):
    A_fntt_n_matrix[i] = FNTT_DIF_NR(A_fntt_n_matrix[i],w_m_table,q)
print_matrix_inmem_col(A_fntt_n_matrix,print_dim,print_mat,False,True,-1,"FNTT")

# 2.4. HP(R-->R)
#      Hadamard Product 
print("--->>> Performing Hadamard Product (sequentially):")     
for i in range(m):
    for j in range(k):
        A_fntt_n_matrix[j][i] = (A_fntt_n_matrix[j][i]*W_H_fntt[j][i])%q
    print_matrix_inmem_col(A_intt_n_rearranged_matrix,print_dim,print_mat,False,True,i," HP ")

# 2.5. !!! (Here we use transpose just to mimic data movement - this happens naturally in HW) !!!
A_fntt_n_matrix = [list(i) for i in zip(*A_fntt_n_matrix)]

# 2.6. UNROLLNTT(R-->R) 
#      Apply NTT to m columns of k-sized vector
print("--->>> Applying {}-pt FNTTs to {} columns using unrolled NTT approach (sequentially):".format(k,m))
for i in range(m):
    A_fntt_n_matrix[i] = arr_reverse(FNTT_DIF_NR(A_fntt_n_matrix[i],w_k_table,q))
    print_matrix_inmem_row(A_fntt_n_matrix,print_dim,print_mat,False,True,i,"FNTT")

# 2.7. Write as N
#      This can be done writing data in proper order (without extra hw/cycle cost)
print("--->>> Output in matrix form (re-ordered - this can be done writing data to memory in proper order in HW):")
A_fntt_n_matrix = [list(i) for i in zip(*A_fntt_n_matrix)]
for i in range(k):
    A_fntt_n_matrix[i] = arr_reverse(A_fntt_n_matrix[i])
A_fntt_n_matrix = [list(i) for i in zip(*A_fntt_n_matrix)] # mimicing data write
print_matrix_inmem_row(A_fntt_n_matrix,print_dim,print_mat,True)

# 2.8. Convert from matrix form to polynomial form (in normal order)
A_fntt_n_2 = arr_matrix2array(A_fntt_n_matrix)

# Sanity check
if sum([abs(p0-p1) for p0,p1 in zip(A_fntt_n,A_fntt_n_2)]) == 0:
    print("MFNTT TEST: OK")
else:
    print("MFNTT TEST: FAIL")

print()
