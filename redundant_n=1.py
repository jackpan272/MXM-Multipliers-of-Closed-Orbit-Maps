import sympy as sp

def all_period2_candidates(a, a_conj):
    z = sp.symbols('z')
    F = z * (z - a) / (1 - z * a_conj)
    numer = sp.fraction(sp.simplify(F.subs(z, F) - z))[0]
    all_pts = sp.solve(numer, z)
    # Optionally remove duplicate/strange symbolic solutions, or zeroth point
    return [pt for pt in all_pts if pt != 0]

a, a_conj = sp.symbols('a a_conj')
final_pts = all_period2_candidates(a, a_conj)
print("All roots of F(F(z)) = z (excl. 0):", final_pts)

def F(w):
    return w * (w - a) / (1 - w * a_conj)

z = sp.symbols('z')
F2 = F(F(z))
F2_prime = sp.simplify(sp.diff(F2, z))

multipliers = [sp.simplify(F2_prime.subs(z, pt)) for pt in final_pts]

if len(multipliers) >= 3:
    m1, m2, m3 = multipliers[:3]  # first three points
    e1 = sp.simplify(m1 + m2 + m3)
    e2 = sp.simplify(m1*m2 + m2*m3 + m1*m3)
    e3 = sp.simplify(m1*m2*m3)

    jacobian = sp.Matrix([
        [sp.diff(e1, a), sp.diff(e1, a_conj)],
        [sp.diff(e2, a), sp.diff(e2, a_conj)],
        [sp.diff(e3, a), sp.diff(e3, a_conj)]
    ])

    print("Multipliers at points:", [m1, m2, m3])
    print("e1 =", e1, ", e2 =", e2, ", e3 =", e3)
    print("Jacobian matrix:\n", jacobian)
    rank = jacobian.rank()
    rows, cols = jacobian.shape
    is_full_rank = (rank == min(rows, cols))
    print(f"Jacobian rank: {rank}, full rank: {is_full_rank}")
else:
    print(f"ERROR: Only {len(multipliers)} points found, need at least 3.")
