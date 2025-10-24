import sympy as sp
from itertools import combinations

def sym(poly, m0):
    m = sp.Symbol('m')
    if not isinstance(poly, sp.Poly):
        poly = sp.Poly(poly, m)
    N = poly.degree()
    coeffs = poly.all_coeffs()
    sig = coeffs[:]
    res = [0] * (N + 1)

    # Start from res[N] = sig[N]
    res[N] = sig[N]/m0

    # Compute res from N-1 down to 0
    for i in range(N-1, 0, -1):
        res[i] = sp.simplify(((sig[i] - res[i+1]) / m0))
    #res[0] = sig[0]-m0
    return res

num = 7
m = sp.symbols('m1:8')  # m1, m2, ..., m7

# Function to compute k-th elementary symmetric polynomial
def elementary_symmetric(vars, k):
    return sum(sp.prod(c) for c in combinations(vars, k))

# Generate all elementary symmetric polynomials
E = [elementary_symmetric(m, i) for i in range(1, len(m) + 1)]
#print(E)
   
def test_sym():
    m = sp.symbols('m')
    #m1,m2, m3 = sp.symbols('m1 m2 m3')
    #poly = (m1+m2+m3)*m**2 + (m1*m2 + m2*m3 + m3*m1)*m + (m1*m2*m3)
    #m0 = m1
    # Define variables
    m1, m2, m3, m4 = sp.symbols('m1 m2 m3 m4')

    # Define elementary symmetric polynomials
    #e1 = m1 + m2 + m3 + m4
    #e2 = m1*m2 + m1*m3 + m1*m4 + m2*m3 + m2*m4 + m3*m4
    #e3 = m1*m2*m3 + m1*m2*m4 + m1*m3*m4 + m2*m3*m4
    #e4 = m1*m2*m3*m4
    poly = 0
    for i in range(num):
        poly += m**(num-1-i) * E[i]
    m0 = m1
    result = sym(poly, m0)
    print("Symmetric polynomial coefficients:", result)

test_sym()