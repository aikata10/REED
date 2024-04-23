
from helper import *

USE_TOY_EXAMPLE = True

# -- Parameters
# Scheme parameters
if USE_TOY_EXAMPLE:
    n       = 64          # polynomial size (n)
    n_dim   = [16, 4]
    q       = 257         # prime
    psi     = 26          # 2n-th root of unity
    psi_inv = 89          # modular inverse of 2n-th root of unity (in mod q)
    w       = 162         # n-th root of unity
    w_inv   = 211         # modular inverse of n-th root of unity (in mod q)
    n_inv   = pow(n,-1,q) # modular inverse of n (in mod q)

    m       = n_dim[0]
    k       = n_dim[1]
else:
    n       = 65536
    n_dim   = [4096, 16]
    q       = 9007199256051713 # prime
    psi     = 178665121121     # 2n-th root of unity
    psi_inv = 3740941444629126 # modular inverse of 2n-th root of unity (in mod q)
    w       = 8579309779131170 # n-th root of unity
    w_inv   = 6964258013800667 # modular inverse of n-th root of unity (in mod q)
    n_inv   = pow(n,-1,q)      # modular inverse of n (in mod q)

    m       = n_dim[0]
    k       = n_dim[1]

assert(n == (n_dim[0]*n_dim[1]))
assert((q % (2*n)) == 1)
assert((psi**(2*n) % q) == 1)
assert((w**n % q) == 1)
assert((psi*psi_inv % q) == 1)
assert((w*w_inv % q) == 1)

print("n      : {}".format(n))
print("n_dim  : {}".format(n_dim))
print("q      : {}".format(q))
print("psi    : {}".format(psi))
print("psi_inv: {}".format(psi_inv))
print("w      : {}".format(w))
print("w_inv  : {}".format(w_inv))
print("n_inv  : {}".format(n_inv))
print()
if USE_TOY_EXAMPLE: print("-"*119)
if USE_TOY_EXAMPLE: print()

# Print
print_dim = False
print_mat = USE_TOY_EXAMPLE

# Twiddle factors for large FNTT/INTT
psi_table     = [1] * n
psi_inv_table = [1] * n
for i in range(1,n):
    psi_table[i]     = (psi_table[i-1]*psi) % q
    psi_inv_table[i] = (psi_inv_table[i-1]*psi_inv) % q

w_table     = [1] * (n//2)
w_inv_table = [1] * (n//2)
for i in range(1,n//2):
    w_table[i]     = (w_table[i-1]*w) % q
    w_inv_table[i] = (w_inv_table[i-1]*w_inv) % q

# Twiddle factors for small FNTT/INTT
w_m = (w**k) % q
w_k = (w**m) % q
w_m_inv = (w_inv**k) % q
w_k_inv = (w_inv**m) % q

w_m_table     = [1] * (m//2)
w_m_inv_table = [1] * (m//2)
for i in range(1,m//2):
    w_m_table[i]     = (w_m_table[i-1]*w_m) % q
    w_m_inv_table[i] = (w_m_inv_table[i-1]*w_m_inv) % q

w_k_table     = [1] * (k//2)
w_k_inv_table = [1] * (k//2)
for i in range(1,k//2):
    w_k_table[i]     = (w_k_table[i-1]*w_k) % q
    w_k_inv_table[i] = (w_k_inv_table[i-1]*w_k_inv) % q

# Hadamard matrix 
# (for MINTT, m columns of k-sized vectors)
# (for FINTT, k columns of m-sized vectors)
W_H_intt = [[1 for i in range(k)] for j in range(m)]
for i in range(1,m):
    W_m = (w ** i) % q
    for j in range(1,k):
        W_H_intt[i][j] = (W_H_intt[i][j-1] * W_m) % q

W_H_fntt = [list(i) for i in zip(*W_H_intt)]
for i in range(k):
    W_H_fntt[i] = arr_reverse(W_H_fntt[i])

# Pre/post-processing array
# (for MINTT, k columns of m-sized vectors)
# (for MFNTT, k columns of m-sized vectors)
PP_intt = [[1 for i in range(m)] for j in range(k)]
for i in range(0,k):
    psi_inv_k = (psi_inv ** k) % q
    PP_intt[i][0] = (psi_inv ** i)*n_inv % q
    for j in range(1,m):
        PP_intt[i][j] = (PP_intt[i][j-1] * psi_inv_k) % q
    PP_intt[i] = arr_reverse(PP_intt[i])

PP_fntt = [[1 for i in range(m)] for j in range(k)]
for i in range(0,k):
    psi_k = (psi ** k) % q
    PP_fntt[i][0] = (psi ** i) % q
    for j in range(1,m):
        PP_fntt[i][j] = (PP_fntt[i][j-1] * psi_k) % q
