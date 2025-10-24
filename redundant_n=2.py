import sympy as sp

def all_period2_candidates_n2(a1, a2, a1_conj, a2_conj):
    z = sp.symbols('z')
    # Blaschke product F(z) for n=2
    F = z * (z - a1) / (1 - z * a1_conj) * (z - a2) / (1 - z * a2_conj)
    numer = sp.fraction(sp.simplify(F.subs(z, F) - z))[0]
    all_pts = sp.solve(numer, z)
    # Remove 0, if present
    return [pt for pt in all_pts if pt != 0]

a1, a2, a1_conj, a2_conj = sp.symbols('a1 a2 a1_conj a2_conj')
final_pts = all_period2_candidates_n2(a1, a2, a1_conj, a2_conj)
print("All roots of F(F(z)) = z (excl. 0):", final_pts)

def F(w):
    return w * (w - a1) / (1 - w * a1_conj) * (w - a2) / (1 - w * a2_conj)

z = sp.symbols('z')
F2 = F(F(z))
F2_prime = sp.simplify(sp.diff(F2, z))

multipliers = [sp.simplify(F2_prime.subs(z, pt)) for pt in final_pts]

if len(multipliers) >= 3:
    m1, m2, m3 = multipliers[:3]
    e1 = sp.simplify(m1 + m2 + m3)
    e2 = sp.simplify(m1*m2 + m2*m3 + m1*m3)
    e3 = sp.simplify(m1*m2*m3)

    diff_params = [a1, a2, a1_conj, a2_conj]
    jacobian = sp.Matrix([
        [sp.diff(e1, param) for param in diff_params],
        [sp.diff(e2, param) for param in diff_params],
        [sp.diff(e3, param) for param in diff_params],
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
