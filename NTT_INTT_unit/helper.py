
from colorama import init as colorama_init
from colorama import Fore
from colorama import Style
from math import log,ceil

VERBOSE = 0

# ----------------------------------------------------------------------------
# Bit-Reverse integer
def int_reverse(a,n):
    b = ('{:0'+str(n)+'b}').format(a)
    return int(b[::-1],2)

# Bit-Reverse array
def arr_reverse(a):
    n = len(a)
    b = [0]*n
    log_n = int(log(n,2))
    for i in range(len(a)):
        i_rev = int_reverse(i,log_n)
        b[i_rev] = a[i]
    return b

# Re-arrange (order) array as [0,1,...,n-2,n-1] --> [0,n-1,...,2,1]
def arr_rearrange(a):
    b = [a[0]] + list(reversed(a[1:]))
    return b

# Bit-reverse and Re-arrange (order) array
def arr_reverse_rearrange(a,n):
    b = arr_reverse(a,n)
    c = arr_rearrange(b)
    return c

# Div-by-n array
def arr_divbyn(a,n,q):
    n_inv = modinv(n,q)
    b = [(i*n_inv)%q for i in a]
    return b

# array to matrix (n -->> m * k, m:row, k:column)
def arr_array2matrix(a,m,k):
    # Converts an array of size n=m*k into matrix [m x k] in row-major order
    # Output will be an array of matrix columns
    """
    Example: N = m * k = 2 * 4
    [0 1 2 3 4 5 6 7] = [0 1 2 3]
                        [4 5 6 7]
    In:  [0 1 2 3 4 5 6 7]
    Out: [[0,4], [1,5], [2,6], [3,7]]
    """

    b = [[0 for i in range(m)] for j in range(k)]

    for i in range(m):
        for j in range(k):
            b[j][i] = a[i*k+j]

    return b

# matrix to array (m * k -->> n, m:row, k:column)
def arr_matrix2array(a):
    # Converts a matrix of size [m x k] in row-major order into an array of size n=m*k
    # Input is an array of matrix columns
    """
    Example: N = m * k = 2 * 4
    [0 1 2 3] = [0 1 2 3 4 5 6 7]
    [4 5 6 7]
    In:  [[0,4], [1,5], [2,6], [3,7]]
    Out: [0 1 2 3 4 5 6 7]
    """
    k = len(a)    # number of columns
    m = len(a[0]) # number of rows

    b = [0]*(m*k)

    for i in range(m):
        for j in range(k):
            b[i*k+j] = a[j][i]
    
    return b

def print_matrix(a, print_dim, print_mat):
    # prints [m x k] matrix
    # Input is an array of matrix columns
    k = len(a)    # number of columns
    m = len(a[0]) # number of rows

    S = ""
    for i in range(m):
        for j in range(k):
            S = S + ("{}".format(a[j][i])).ljust(10)
        S = S + "\n"

    if print_dim: print("Matrix dimension (row x column): {} x {}".format(m,k))
    if print_mat: print(S)

def print_matrix_inmem_row(a, print_dim, print_mat, is_in_memory=False, print_step=False, step_num=0, step_label=""):
    # prints [m x k] matrix
    # Input is an array of matrix columns
    # Columns are distributed over memories
    k = len(a)    # number of columns
    m = len(a[0]) # number of rows

    S = ""
    if print_step:
        for j in range(step_num):
            S = S + " ".ljust(7)
        S = S + step_label + "\n"
        for j in range(step_num):
            S = S + " ".ljust(7)
        S = S + " || \n"
        for j in range(step_num):
            S = S + " ".ljust(7)
        S = S + " \/ \n"

    if is_in_memory:
        S = S + "addr:".ljust(7)
        for j in range(k):
            S = S + ("{}".format(j)).ljust(7) 
        S = S + "\n"

    for i in range(m):
        if is_in_memory:
            for j in range(k+1):
                S = S + "-"*7
        else:
            for j in range(k):
                S = S + "-"*7
        S = S + "\n"

        if is_in_memory: 
            S = S + ("MEM{}  ".format(i)).ljust(6) + "|"
        else:
            S = S + "|"

        for j in range(k):
            t_str = ("{}".format(a[j][i])).ljust(6)
            if print_step and (j == step_num):
                S = S + (f"{Fore.GREEN}"+t_str+f"{Style.RESET_ALL}").ljust(6) + "|"
            else:
                S = S + ("{}".format(a[j][i])).ljust(6) + "|"
        S = S + "\n"
    
    if is_in_memory:
        for j in range(k+1):
            S = S + "-"*7
    else:
        for j in range(k):
            S = S + "-"*7
    S = S + "\n"

    if print_dim: print("Matrix dimension (row x column): {} x {}".format(m,k))
    if print_mat: print(S)

def print_matrix_inmem_col(a, print_dim, print_mat, is_in_memory=False, print_step=False, step_num=0, step_label=""):
    # prints [m x k] matrix
    # Input is an array of matrix columns
    # Each column is in one memory
    m = len(a)    # number of columns
    k = len(a[0]) # number of rows

    S = ""
    if print_step and (step_num != -1):
        for j in range(step_num):
            S = S + " ".ljust(7)
        S = S + step_label + "\n"
        for j in range(step_num):
            S = S + " ".ljust(7)
        S = S + " || \n"
        for j in range(step_num):
            S = S + " ".ljust(7)
        S = S + " \/ \n"

    if is_in_memory:
        S = S + "addr:".ljust(7)
        for j in range(k):
            S = S + ("{}".format(j)).ljust(7) 
        S = S + "\n"

    for i in range(m):
        if is_in_memory:
            for j in range(k+1):
                S = S + "-"*7
        else:
            for j in range(k):
                S = S + "-"*7
        S = S + "\n"

        if print_step and (step_num == -1):
            if is_in_memory: 
                S = S + ("MEM{}  ".format(i)).ljust(6) + "|"
            else:
                S = S + "|"

            for j in range(k):
                t_str = ("{}".format(a[i][j])).ljust(6)
                S = S + (f"{Fore.GREEN}"+t_str+f"{Style.RESET_ALL}").ljust(6) + "|"
            S = S + "<= " + step_label + "\n"
        else:
            if is_in_memory: 
                S = S + ("MEM{}  ".format(i)).ljust(6) + "|"
            else:
                S = S + "|"

            for j in range(k):
                t_str = ("{}".format(a[i][j])).ljust(6)
                if print_step and (j == step_num):
                    S = S + (f"{Fore.GREEN}"+t_str+f"{Style.RESET_ALL}").ljust(6) + "|"
                else:
                    S = S + ("{}".format(a[i][j])).ljust(6) + "|"
            S = S + "\n"
    if is_in_memory:
        for j in range(k+1):
            S = S + "-"*7
    else:
        for j in range(k):
            S = S + "-"*7
    S = S + "\n"

    if print_dim: print("Matrix dimension (row x column): {} x {}".format(m,k))
    if print_mat: print(S)

# ----------------------------------------------------------------------------
# Merged NTT with pre-processing (optimized) (iterative)
# This is not NTT, this is pre-processing + NTT
# (see: https://eprint.iacr.org/2016/504.pdf)
# A: input polynomial (standard order)
# Psi: 2n-th root of unity
# q: modulus
# B: output polynomial (bit-reversed order)
def MFNTT_DIT_NR(A,Psi_table,q):
    N = len(A)
    B = [_ for _ in A]

    if VERBOSE: print("--- MFNTT_DIT_NR")

    l = int(log(N,2))

    t = N
    m = 1
    while(m<N):
        t = int(t/2)
        for i in range(m):
            j1 = 2*i*t
            j2 = j1 + t - 1
            Psi_pow = int_reverse(m+i,l)
            #S = (Psi**Psi_pow) % q
            S = Psi_table[Psi_pow]
            for j in range(j1,j2+1):
                U = B[j]
                V = (B[j+t]*S) % q

                B[j]   = (U+V) % q
                B[j+t] = (U-V) % q

                if VERBOSE: print("A={}, B={}, W={}".format(j,j+t,Psi_pow))
        m = 2*m

        if VERBOSE: print("---")

    return B

# Iterative Radix-2 Decimation-in-Frequency (DIF) (GS) NTT - NR
# A: input polynomial (standard order)
# W: twiddle factor
# q: modulus
# B: output polynomial (bit-reversed order)
def FNTT_DIF_NR(A,W_table,q):
    N = len(A)
    B = [_ for _ in A]

    if VERBOSE: print("--- FNTT_DIF_NR")

    m = 1
    v = N

    while v>1:
        s = int(v/2)
        for k in range(m):
            jf = k * v
            jl = jf + s - 1
            jt = 0
            for j in range(jf,jl+1):
                #TW = pow(W,jt,q)
                TW = W_table[jt]

                temp = B[j]
                temp2= B[j+s]

                B[j  ] = (temp + B[j+s]) % q
                B[j+s] = (temp - B[j+s])*TW % q

                if VERBOSE: print("A={}, B={}, W={}".format(j,j+s,jt))
                
                jt = jt + m
        m = 2*m
        v = s

        if VERBOSE: print("---")

    return B

