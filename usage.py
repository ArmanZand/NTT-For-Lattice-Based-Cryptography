from ntt import *
# Toy example and usage
# Basic parameters
n = 8
q = 17
psi = get_prim_root(n,q)
psi_inv = pow(psi,-1,q)

# The two polynomials to multiply with NTT.
A = [1,2,3,4,5,6,7,8]
B = [1,2,3,4,5,6,7,8]

# Precompute powers of psi to speedup main NTT process.
psi_table = [1] * n
psi_inv_table = [1] * n
for i in range(1, n):
    psi_table[i] = ((psi_table[i-1] * psi) % q)
    psi_inv_table[i] = ((psi_inv_table[i-1] * psi_inv) % q)
    
# Change the lists into bit-reverse order.
psi_table = bit_reverse_order(psi_table)
psi_inv_table = bit_reverse_order(psi_inv_table)

# Forward Number Theoretic Transform of arrays A and B.
A = speedup_FNTT(A, q, psi_table)
B = speedup_FNTT(B, q, psi_table)
# Multiply the elements in mod q.
AB = [(x*y) % q for x,y in zip(A,B)]
# Reverse Number Theoretic Transform the product.
C = speedup_INTT(AB, q, psi_inv_table)
print(C)