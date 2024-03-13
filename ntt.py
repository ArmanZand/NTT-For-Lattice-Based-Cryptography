#Implementation of the NTT algorithms described by the paper:
#Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography
#Patrick Longa and Michael Naehrig
#Microsoft Research, USA
#{plonga,mnaehrig}@microsoft.com

import math, os


def secure_next_int( min: int, max: int) -> int:
    '''
    Find and returns a random integer between the given range with cryptographically secure RNG.
    
    Parameters
    ----------
    min : int
        Minimum possible random value to be returned.
    max : int
        Maximum possible random value to be returned.    
    
    Returns
    -------
    result : int
        The cryptographically secure random integer between the given range.
    
    Examples
    --------
    >>> secure_next_int(2,23)
    9
    '''
    if min == max: return min
    if min > max: raise Exception("Random number range max cannot be less than min.")
    i_range = max - min
    if i_range == 0: return min
    mask = i_range
    numSize = (i_range.bit_length() +7) // 8
    iter = int(math.log(numSize*8, 2))
    for i in range(iter):
        mask |= mask >> pow(2, i)
    result = 0
    while True:
        result = int.from_bytes(os.urandom(numSize), byteorder="little") & mask
        if (result <= i_range): break
    result += min
    return result

def get_prim_root(n, q, r=None):
    '''
    Finds a possible primitive nth root of unity given n and mod q.
    This is done by checking if a random integer x, x^{n/2} mod q is
    congruent to q - 1. The function will try 50 random integers to
    find the primitive nth root of unity.
    
    Parameters
    ----------
    n : int
        The length of the polynomials which should be some power of 2.
    q : int
        The modulus value.
    r : int
        Factor of n to q - 1.
        
    Returns
    -------
    candidate : int
        A primitive nth root of unity discovered.
        
    Examples
    --------
    >>> get_prim_roots(8,17)
    6
    '''
    n = 2*n
    if r==None:
        r = (q -1) // n    
    tries = 50
    for _ in range(tries):
        uri = secure_next_int(2,q-1)
        candidate = pow(uri,r, q)
        if pow(candidate, n//2, q) == q-1:
            return candidate
    raise Exception("Failed to find primitive nth root of unity.")

def speedup_FNTT(a, q, psi_table_rev):
    #NTT based on the Cooley-Tukey butterfly
    '''
    Implementation of the Foward NTT based on the Cooley-Tukey
    butterfly. A precomputed table of the powers of psi in
    bit-reversed order. Compared to other algorithms, the
    input array and output array does not need to be converted
    in bit-reversed order.
    
    Parameters
    ----------
    a : List[int]
        The array to forward transform.
    q : int
        Modulus value.
    psi_table_rev : List[int]
        List of the powers of psi in bit-reversed order.
        
    Returns
    -------
    A : List[int]
        The forward number theoretic transform of the input a.
    '''
    n = len(a)
    t = n
    A = [x for x in a]
    m = 1
    while m < n:
        t >>= 1
        for i in range(m):
            j_1 = 2 * i * t
            j_2 = j_1 + t - 1
            S = psi_table_rev[m + i]
            for j in range(j_1, j_2+1):
                U = A[j]
                V = A[j + t] * S
                A[j] = (U + V) % q
                A[j + t] = (U - V) % q
        m<<=1
    return A
def speedup_INTT(a, q, psi_inv_table_rev):
    '''
    Implementation of the Inverse NTT based on the Gentleman-Sande 
    butterfly. Similarly to the Forward NTT method `speedup_FNTT()`, 
    a precomputed table of the powers of psi in bit-reversed order 
    is used. The inputs and outputs also don't need their bit orders
    changed.
    
    Parameters
    ----------
    a : List[int]
        An array of integers that has been foward transformed.
    q : int
        Modulus value.
    psi_inv_table_rev : List[int]
        An array of the powers of the  modular multiplicative inverse 
        of psi.

    Returns
    -------
    A : List[int]
        The inverse NTT of the given input array a.
    '''
    n = len(a)
    A = [x for x in a]
    t = 1
    m = n
    while m > 1:
        j_1 = 0
        h = m >> 1
        for i in range(h):
            j_2 = j_1 + t - 1
            S = psi_inv_table_rev[h + i]
            for j in range(j_1, j_2+1):
                U = A[j]
                V = A[j + t]
                A[j] = (U + V) % q
                A[j + t] = ((U - V) * S) % q
            j_1 += t << 1
        t<<=1
        m>>=1
    n_inv = pow(n,-1, q)
    for j in range(n):
        A[j] = (A[j] * n_inv) % q
    return A
    
def bit_reverse_order(a):
    '''
    Reorders the given array in reverse-bit order.
    
    Parameters
    ----------
    a : List
        Any arbitrary array.
    
    Returns
    -------
    result : List
        Given array a where the elements are ordered in bit-reverse.
        
    Examples
    --------
    >>> a = [1, 2, 3, 4, 5, 6, 8]
    >>> bit_reverse_order(a)
    [1, 5, 3, 7, 2, 6, 4, 8]
    '''
    num_bits = len(bin(len(a) - 1)) - 2
    result = [0] * len(a)
    for i in range(len(a)):
        rev_index = int(bin(i)[2:].zfill(num_bits)[::-1], 2)
        result[rev_index] = a[i]
    return result
