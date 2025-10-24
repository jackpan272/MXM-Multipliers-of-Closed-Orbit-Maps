import sympy as sp
    #find period 2 pts when n=1
def chenxi_subtracted_period2_points_n1(a, a_conj):
    z = sp.symbols('z')
    # Define F(z)
    F = z * (z - a) / (1 - z * a_conj)
    # F(F(z)) - z = 0
    numer = sp.simplify(F.subs(z, F) - z)
    all_pts = sp.solve(numer, z)
    # Fixed points
    fixed_pts = sp.solve(sp.simplify(F - z), z)
    # Non-fixed period-2 points
    period2_pts = [pt for pt in all_pts if not any(sp.simplify(pt - fp) == 0 for fp in fixed_pts)]
    # Find representatives to exclude
    representatives = set()
    used = set()
    for pt in period2_pts:
        if pt in used:
            continue
        F_at_pt = sp.simplify(F.subs(z, pt))
        if F_at_pt in period2_pts and F_at_pt != pt:
            # Choose a representative deterministically
            pair = tuple(sorted([pt, F_at_pt], key=lambda x: sp.expand_complex(x)))
            representatives.add(pair[0])  # always omit the first from the pair
            used.add(pt)
            used.add(F_at_pt)
    # Return: non-fixed period-2 points minus representatives
    final_period2_pts = [pt for pt in all_pts if pt not in representatives]
    return final_period2_pts

a, a_conj = sp.symbols('a a_conj')
final_pts = chenxi_subtracted_period2_points_n1(a, a_conj)
#print("Final period-2 points after Chenxi subtraction (n=1):", final_pts)
print("Final period-2 points after Chenxi subtraction (n=1) length:", len(final_pts))

def F(w):
    return w * (w - a) / (1 - w * a_conj)

z = sp.symbols('z')
F2 = F(F(z))
F2_prime = sp.simplify(sp.diff(F2, z))
print(F2_prime)



# Compute multipliers at the two period-2 points
m1 = sp.simplify(F2_prime.subs(z, final_pts[0]))
m2 = sp.simplify(F2_prime.subs(z, final_pts[1]))

# Elementary symmetric polynomials
e1 = sp.simplify(m1 + m2)
e2 = sp.simplify(m1 * m2)

# Jacobian with respect to a and a_conj
jacobian = sp.Matrix([
    [sp.diff(e1, a), sp.diff(e1, a_conj)],
    [sp.diff(e2, a), sp.diff(e2, a_conj)]
])
determinant = sp.simplify(jacobian.det())

actual_final_pts = [final_pts[0]]
for i in range(1, len(final_pts)):
    for j in range(len(actual_final_pts)):
        if final_pts[i] == actual_final_pts[j] or F(final_pts[i]) == actual_final_pts[j]:
            break
        else:
            actual_final_pts.append(final_pts[i])
    #if final_pts[i] not in actual_final_pts and F(final_pts[i]) not in actual_final_pts:
     #   actual_final_pts.append(final_pts[i])
print("Actual final period-2 points after Chenxi subtraction (n=1) length:", len(actual_final_pts))

print((F(final_pts[3])-final_pts[2]).simplify())

for i in range(len(final_pts)):
    print(sp.simplify(F(F(final_pts[i]))-final_pts[i]))

for i in range(len(final_pts)):
    print(sp.simplify(F(final_pts[i])-final_pts[i]))


#print("Multipliers at period-2 points:", m1, m2)
#print("Elementary symmetric polynomials: e1 =", e1, ", e2 =", e2)
#print("Jacobian row number:\n", len(jacobian[:, 0]))
#print("Jacobian column number:\n", len(jacobian[0, :]))
#print("Jacobian determinant:\n", determinant)