import sympy as sp
from itertools import combinations


n = 3

a = sp.symbols('a1:%d' % n) 

a_conj=sp.symbols('a_conj1:%d' % n)

z = sp.symbols('z')


def F(w):
    """F(w) = w * ∏(w - a_i)/(1 - w * a_conj_i)"""
    result = w
    for i in range(len(a)):
        result *= (w - a[i]) / (1 - w * a_conj[i])
    return result

# Define the multiplier (F ∘ F)'(z)
def F_composition_F_prime(z_i):
    """Compute (F∘F)'(z) = derivative of F(F(z)) with respect to z"""
    F_comp = F(F(z_i))
    return sp.diff(F_comp, z_i)


m = int(n**2 + 2*n - (n**2 + n)/2)

# Define elementary symmetric polynomials function
def elementary_symmetric_poly(k, vars):
    """Compute k-th elementary symmetric polynomial in vars"""
    if k == 0:
        return 0
    elif k > len(vars):
        return 0
    else:
        return sum(sp.prod(comb) for comb in combinations(vars, k))

# Define multiplier variables (representing multipliers at period-2 points)
multipliers = sp.symbols('x0:%d'%m)

#Need to find the algebraic relations between the multipliers and the parameters a_i, a_conj_i

# Pass the multipliers in the elementary symmetric polynomials
sym_polys = [elementary_symmetric_poly(k, multipliers) for k in range(1, m)]

# Take the jacobian with respect to the parameters a_1 through a_n
all_parameters = list(a) + list(a_conj)


# Jacobian matrix
jacobian_matrix = []
for poly in sym_polys[:m]:  # Just first 3 polynomials to avoid clutter
    jacobian_row = []
    for param in all_parameters:
        # This would be the partial derivative ∂poly/∂param
        # But since poly contains symbolic multipliers, this is zero
        # In the actual computation, you'd substitute the multiplier expressions first
        partial_deriv = sp.diff(poly, param)
        jacobian_row.append(partial_deriv)
    jacobian_matrix.append(jacobian_row)