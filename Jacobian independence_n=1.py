'''
this code doesn't exclude 0 as a fixed point 
import sympy as sp
    #find period 2 pts when n=1
def chenxi_subtracted_period2_points_n1(a, a_conj):
    z = sp.symbols('z')
    # Define F(z)
    F = z * (z - a) / (1 - z * a_conj)
    # F(F(z)) - z = 0
    numer = sp.fraction(sp.simplify(F.subs(z, F) - z))[0]
    all_pts = sp.solve(numer, z)
    # Fixed points
    fixed_pts = sp.solve(sp.simplify(F - z), z)
    # Non-fixed period-2 points
    nonfixedperiod2_pts = [pt for pt in all_pts if not any(sp.simplify(pt - fp) == 0 for fp in fixed_pts)]
    # Find representatives to exclude
    representatives = set()
    used = set()
    for pt in nonfixedperiod2_pts:
        if pt in used:
            continue
        F_at_pt = sp.simplify(F.subs(z, pt))
        if F_at_pt in nonfixedperiod2_pts and F_at_pt != pt:
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
print("Final period-2 points after Chenxi subtraction (n=1):", final_pts)
# Print the number of final points using a for-loop and list them with indices
count = 0
for i, pt in enumerate(final_pts, start=1):
    print(f"Point {i}:", pt)
    count += 1
print("Number of final points:", count)

def F(w):
    return w * (w - a) / (1 - w * a_conj)

z = sp.symbols('z')
F2 = F(F(z))
F2_prime = sp.simplify(sp.diff(F2, z))

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

print("Multipliers at period-2 points:", m1, m2)
print("Elementary symmetric polynomials: e1 =", e1, ", e2 =", e2)
print("Jacobian matrix:\n", jacobian)
print("Jacobian determinant:\n", determinant)

def period2_pts_verification(a, a_conj):
    import sympy as sp
    z = sp.symbols('z')
    F = z * (z - a) / (1 - z * a_conj)
    F2 = F.subs(z, F)
    numer = sp.fraction(sp.simplify(F2 - z))[0]
    all_pts = sp.solve(numer, z)
    period2_pts = []
    for pt in all_pts:
        test1 = sp.simplify(F2.subs(z, pt) - pt)
        if test1 == 0:
            period2_pts.append(pt)
    return period2_pts

# Example usage:
a, a_conj = sp.symbols('a a_conj')
pts = period2_pts_verification(a, a_conj)
print("period-2 points:", pts)
'''

import sympy as sp

def chenxi_subtracted_period2_points_n1(a, a_conj):
    z = sp.symbols('z')
    F = z * (z - a) / (1 - z * a_conj)
    numer = sp.fraction(sp.simplify(F.subs(z, F) - z))[0]
    all_pts = sp.solve(numer, z)
    # Fixed points
    fixed_pts = sp.solve(sp.simplify(F - z), z)
    # Remove 0 from fixed points explicitly
    fixed_pts = [pt for pt in fixed_pts if pt != 0]
    # Non-fixed period-2 points (removing 0 just to be safe)
    nonfixedperiod2_pts = [pt for pt in all_pts if not any(sp.simplify(pt - fp) == 0 for fp in fixed_pts) and pt != 0]
    # Find representatives to exclude
    representatives = set()
    used = set()
    for pt in nonfixedperiod2_pts:
        if pt in used:
            continue
        F_at_pt = sp.simplify(F.subs(z, pt))
        if F_at_pt in nonfixedperiod2_pts and F_at_pt != pt:
            pair = tuple(sorted([pt, F_at_pt], key=lambda x: sp.expand_complex(x)))
            representatives.add(pair[0])
            used.add(pt)
            used.add(F_at_pt)
    final_period2_pts = [pt for pt in all_pts if pt not in representatives and pt != 0]
    return final_period2_pts

a2, a2_conj = sp.symbols('a a_conj')
print(len(chenxi_subtracted_period2_points_n1(a2,a2_conj)))

a, a_conj = sp.symbols('a a_conj')
final_pts = chenxi_subtracted_period2_points_n1(a, a_conj)
print("Final period-2 points after Chenxi subtraction (n=1):", final_pts)



#try not reducing the points and check if it's full rank



def F(w):
    return w * (w - a) / (1 - w * a_conj)

z = sp.symbols('z')
F2 = F(F(z))
F2_prime = sp.simplify(sp.diff(F2, z))

m1 = sp.simplify(F2_prime.subs(z, final_pts[0]))
m2 = sp.simplify(F2_prime.subs(z, final_pts[1]))

e1 = sp.simplify(m1 + m2)
e2 = sp.simplify(m1 * m2)

jacobian = sp.Matrix([
    [sp.diff(e1, a), sp.diff(e1, a_conj)],
    [sp.diff(e2, a), sp.diff(e2, a_conj)]
])
determinant = sp.simplify(jacobian.det())

print("Multipliers at period-2 points:", m1, m2)
print("Elementary symmetric polynomials: e1 =", e1, ", e2 =", e2)
print("Jacobian matrix:\n", jacobian)
print("Jacobian determinant:\n", determinant)

def period2_pts_verification(a, a_conj):
    z = sp.symbols('z')
    F = z * (z - a) / (1 - z * a_conj)
    F2 = F.subs(z, F)
    numer = sp.fraction(sp.simplify(F2 - z))[0]
    all_pts = sp.solve(numer, z)
    period2_pts = []
    for pt in all_pts:
        test1 = sp.simplify(F2.subs(z, pt) - pt)
        if test1 == 0 and pt != 0:  # Remove 0 here too
            period2_pts.append(pt)
    return period2_pts

pts = period2_pts_verification(a, a_conj)
print("period-2 points:", pts)
