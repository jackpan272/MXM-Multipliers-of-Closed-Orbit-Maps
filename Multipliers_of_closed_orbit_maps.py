import sympy as sp


n = 3




a = sp.symbols('a1:%d' % n) 
a_conj=sp.symbols('a_conj1:%d' % n)
z = sp.symbols('z')

#define F as a function not an expression to allow for composition 

'''F = z
for i in range(len(a)):
    F *= (z - a[i]) / (1 - z * a_conj[i])
    '''

def FoF(z):
    F =z
    for i in range(len(a)):
        F *= (z - a[i]) / (1 - z * a_conj[i])
    FoF=F.subs(z, F)

    return FoF


F_comp = F.subs(z, F)


F_comp_prime = sp.diff(F_comp, z)


m = int(sp.simplify(n**2 + 2*n - (n**2 + n)/2))


x = sp.symbols('x0:%d' % m)


def symmetric_poly(k, vars):
    if k == 0:
        return 1
    elif k > len(vars):
        return 0
    else:
        from itertools import combinations
        return sum(sp.prod(comb) for comb in combinations(vars, k))
    

c=sp.symbols('c')
print(symmetric_poly(4, z)) 


sym_polys = [symmetric_poly(k, x) for k in range(1, m+1)]


z_points = sp.symbols('z0:%d' % m)
subs_values = [F_comp_prime.subs(z, zp) for zp in z_points]

subs_dict = {x[i]: subs_values[i] for i in range(m)}


sym_polys_evaluated = [p.subs(subs_dict) for p in sym_polys]


print("F(z) =", F)
print("(F∘F)'(z) =", F_comp_prime)
print("Symmetric polynomials (evaluated):")
for i, poly in enumerate(sym_polys_evaluated, 1):
    print(f"e_{i} =", poly)

    #write down the jacobian of the sym_polys_evaluated wrt a and a_conj
    #Need period two points to calculate the Jacobian, don't have the period two periods, why?
    #might have to write the period two points in terms of a and a_conj
    #see if the period two points are independent of a and a_conj