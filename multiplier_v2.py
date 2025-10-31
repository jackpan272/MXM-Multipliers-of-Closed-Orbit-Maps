import sympy as sp

z, m = sp.symbols('z m')
n = 1
a = sp.symbols(f'a1:{n+1}', complex=True)
#a_conj = [sp.conjugate(ai) for ai in a]
a_bar = sp.symbols(f'a_bar1:{n+1}', complex=True)
#NUMERICAL VALUES FOR TESTING


def F(w):
    expr = w
    for i in range(n):
        expr *= (w - a[i]) / (1 - w * a_bar[i])
    return expr

F2 = F(F(z))
F2 = F2.ratsimp()
F2 = sp.cancel(F2)
# Use together() to combine over common denominator, then extract numerator
num_A = sp.numer(F2)-z*sp.denom(F2)
#print(num_A)
poly_A = sp.Poly(num_A, z, m)
print("Degree of A:", sp.degree(poly_A))

# F'
#F_prime = sp.simplify(sp.diff(Fz, z))

# F(F(z))
#F2 = sp.simplify(F(Fz))

# (F²)'(z)

F2_prime = sp.diff(F2, z)# get the polynomial
F2_prime = F2_prime.ratsimp()
F2_prime = sp.cancel(F2_prime)

num_F2prime = sp.numer(F2_prime)
den_F2prime = sp.denom(F2_prime)
#num_B = sp.numer(sp.together(F2_prime ))
poly_B = sp.Poly(num_F2prime-m*den_F2prime, z, m)
print("Degree of B:",sp.degree(poly_B))
'''
print("Computing F'(z)...")
F_prime_z = sp.diff(Fz, z)

print("Computing F'(F(z))...")
F_prime_Fz = F_prime_z.subs(z, Fz)

print("Applying chain rule: F'(F(z)) * F'(z)...")
F2_prime = F_prime_Fz * F_prime_z

print("Extracting numerator...")
num_B = sp.numer(sp.together(F2_prime - m))
'''

#resultant = poly_A.resultant(poly_B, z)
#Enter numerical values for parameters OR calculate resultant from scratch
resultant = sp.resultant(poly_A, poly_B, z)
resultant = sp.Poly(resultant, m)
print("Degree of resultant:",sp.degree(resultant))