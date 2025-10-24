import sympy as sp

def chenxi_subtracted_period2_points_n2(a1, a2, a1_conj, a2_conj):
    z = sp.symbols('z')
    # Define F(z) for n=2
    F = z
    F *= (z - a1) / (1 - z * a1_conj)
    F *= (z - a2) / (1 - z * a2_conj)
    # F(F(z)) - z = 0
    numer = sp.fraction(sp.simplify(F.subs(z, F) - z))[0]
    #figure out complex solutions
    all_pts = sp.solve(numer, z)
    # Fixed points
    fixed_pts = sp.solve(sp.simplify(F - z), z)
    fixed_pts = [pt for pt in fixed_pts if pt != 0]
    # Non-fixed period-2 points
    period2_pts = [pt for pt in all_pts if not any(sp.simplify(pt - fp) == 0 for fp in fixed_pts) and pt != 0]
    # Find representatives to exclude
    representatives = set()
    used = set()
    for pt in period2_pts:
        if pt in used:
            continue
        F_at_pt = sp.simplify(F.subs(z, pt))
        if F_at_pt in period2_pts and F_at_pt != pt:
            pair = tuple(sorted([pt, F_at_pt], key=lambda x: sp.expand_complex(x)))
            representatives.add(pair[0])
            used.add(pt)
            used.add(F_at_pt)
    final_period2_pts = [pt for pt in all_pts if pt not in representatives and pt!=0]
    return final_period2_pts

a2, a2_conj,a3,a3_conj = sp.symbols('a2 a2_conj a3 a3_conj')
print(len(chenxi_subtracted_period2_points_n2(a2,a2_conj,a3,a3_conj)))


# Define symbols for n=2
a1, a2, a1_conj, a2_conj = sp.symbols('a1 a2 a1_conj a2_conj')
final_pts = chenxi_subtracted_period2_points_n2(a1, a2, a1_conj, a2_conj)
print("Final period-2 points after Chenxi subtraction (n=2):", final_pts)

# Define F(z) for n=2
z = sp.symbols('z')
def F(w):
    return w * (w - a1) / (1 - w * a1_conj) * (w - a2) / (1 - w * a2_conj)

F2 = F(F(z))
F2_prime = sp.simplify(sp.diff(F2, z))

# Compute multipliers at the period-2 points
multipliers = [sp.simplify(F2_prime.subs(z, pt)) for pt in final_pts]

# Elementary symmetric polynomials
# For m period-2 points, e1 = sum, e2 = sum of products, etc.
e1 = sp.simplify(sum(multipliers))
e2 = sp.simplify(sp.expand(sum(multipliers[i]*multipliers[j] for i in range(len(multipliers)) for j in range(i+1, len(multipliers)))))

# Jacobian with respect to all four parameters
diff_params = [a1, a2, a1_conj, a2_conj]
jacobian = sp.Matrix([
    [sp.diff(e1, param) for param in diff_params],
    [sp.diff(e2, param) for param in diff_params]
])
determinant = sp.simplify(jacobian.det())

print("Multipliers at period-2 points:", multipliers)
print("Elementary symmetric polynomials: e1 =", e1, ", e2 =", e2)
print("Jacobian matrix:\n", jacobian)
print("Jacobian determinant:\n", determinant)
