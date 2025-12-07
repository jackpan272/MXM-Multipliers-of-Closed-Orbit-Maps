"""
ANALYTIC SOLUTION: n=1 case
Compute period-2 points and their multipliers symbolically using SymPy.

For n=1, the map is: f(z) = z * (z - a) / (1 - z*ā)
The period-2 equation f(f(z)) = z is a degree-4 polynomial (quartic),
which is solvable in closed form.
"""

import sympy as sp
from sympy import symbols, expand, simplify, solve, diff, conjugate
import numpy as np

# ============================================================================
# SYMBOLIC COMPUTATION (n=1 case)
# ============================================================================

print("="*80)
print("ANALYTIC n=1 CASE: Period-2 Points and Multipliers via SymPy")
print("="*80)

# Define symbolic variables
z = symbols('z', complex=True)
a = symbols('a', complex=True)

# Define the map f(z) = z * (z - a) / (1 - z*ā)
# For numerical testing, we'll pick a specific value of a after building the formula
numerator = z * (z - a)
denominator = 1 - z * conjugate(a)
f = numerator / denominator

print(f"\nMap: f(z) = z(z - a) / (1 - zā)")

# Compose f(f(z))
f_of_z = f.subs(z, z)  # This is already f(z)
f_of_f = f.subs(z, f)  # f(f(z))

print(f"\nComputing f(f(z))... (this may take a moment)")
f_of_f_simplified = simplify(f_of_f)

# Period-2 equation: f(f(z)) = z
period2_equation = f_of_f - z
period2_equation_simplified = simplify(period2_equation)

print(f"\nPeriod-2 equation: f(f(z)) - z = 0")
print(f"(Simplified form shown below)")

# ============================================================================
# NUMERICAL EXAMPLE: a = 0.7 * exp(i*π/4)
# ============================================================================

print(f"\n{'='*80}")
print(f"NUMERICAL EXAMPLE: a = 0.7 * exp(iπ/4)")
print(f"{'='*80}")

# Specific value: a = 0.7 * e^(i*π/4)
a_val = 0.2 * sp.exp(1j * sp.pi / 8)
a_val_complex = complex(a_val)

print(f"\na = {a_val_complex:.10f}")
print(f"|a| = {abs(a_val_complex):.6f}")

# Substitute a into the period-2 equation
period2_eq_numerical = period2_equation_simplified.subs(a, a_val)

print(f"\nSolving period-2 equation for this value of a...")
print(f"(Finding roots of a quartic polynomial)")

# Solve for z
try:
    period2_roots = solve(period2_eq_numerical, z)
    
    # Filter out z=0 and other trivial solutions, keep period-2 points on |z|=1
    period2_points = []
    for root in period2_roots:
        root_val = complex(root.evalf())
        # Check if on unit circle (|z| ≈ 1) and not trivial
        if abs(abs(root_val) - 1) < 0.01 and abs(root_val) > 1e-6:
            period2_points.append(root_val)
    
    print(f"\nFound {len(period2_points)} period-2 points on the unit circle:")
    for i, z_i in enumerate(period2_points, 1):
        print(f"  z_{i} = {z_i:.10f}  (|z_{i}| = {abs(z_i):.10f})")
    
except Exception as e:
    print(f"Error solving: {e}")
    period2_points = []

# ============================================================================
# MULTIPLIER COMPUTATION
# ============================================================================

if period2_points:
    print(f"\n{'='*80}")
    print(f"MULTIPLIER COMPUTATION: m_j = (f∘f)'(z_j)")
    print(f"{'='*80}")
    
    # Compute the derivative symbolically
    f_of_f_symbolic = f.subs(z, f.subs(z, z))
    f_of_f_derivative = diff(f_of_f_symbolic, z)
    
    print(f"\nComputing d/dz[f(f(z))]... (this may take a moment)")
    f_of_f_derivative_simplified = simplify(f_of_f_derivative)
    
    # Evaluate at each period-2 point
    print(f"\nMultipliers at each period-2 point:")
    print(f"-" * 80)
    
    multipliers = []
    for i, z_i in enumerate(period2_points[:5], 1):  # Show first 5
        try:
            m_i = f_of_f_derivative_simplified.subs([(z, z_i), (a, a_val)])
            m_i_val = complex(m_i.evalf())
            multipliers.append(m_i_val)
            
            print(f"\nz_{i} = {z_i:.10f}")
            print(f"  m_{i} = {m_i_val:.10f}")
            print(f"  |m_{i}| = {abs(m_i_val):.10f}")
            print(f"  Real part: {m_i_val.real:.10f}")
            print(f"  Imag part: {m_i_val.imag:.2e}")
            
        except Exception as e:
            print(f"\nz_{i} = {z_i:.10f}")
            print(f"  Error computing multiplier: {e}")

# ============================================================================
# COMPARISON: Numerical vs Symbolic
# ============================================================================

print(f"\n{'='*80}")
print(f"COMPARISON: Numerical FD vs Symbolic Derivative")
print(f"{'='*80}")

# Numerical finite-difference method for comparison
def map_f_numeric(z_in, a_in):
    """Numerically evaluate f(z) = z(z-a)/(1-zā)"""
    return z_in * (z_in - a_in) / (1 - z_in * np.conj(a_in))

def f_of_f_numeric(z_in, a_in):
    """Numerically evaluate f(f(z))"""
    return map_f_numeric(map_f_numeric(z_in, a_in), a_in)

def deriv_numeric_fd(z_in, a_in, h=1e-8):
    """Numerical derivative via central difference"""
    dz = complex(h, h)
    return (f_of_f_numeric(z_in + dz, a_in) - f_of_f_numeric(z_in - dz, a_in)) / (2 * dz)

if period2_points:
    print(f"\nFor a = {a_val_complex:.6f}:\n")
    for i, z_i in enumerate(period2_points[:3], 1):
        m_numeric = deriv_numeric_fd(z_i, a_val_complex)
        
        # Get symbolic result if available
        if i <= len(multipliers):
            m_symbolic = multipliers[i-1]
            error = abs(m_numeric - m_symbolic)
            
            print(f"z_{i} = {z_i:.8f}")
            print(f"  Symbolic:  {m_symbolic:.10f}")
            print(f"  Numerical: {m_numeric:.10f}")
            print(f"  Error:     {error:.2e}\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print(f"{'='*80}")
print(f"SUMMARY: Analytic Solution for n=1")
print(f"{'='*80}\n")

if period2_points and multipliers:
    print(f"✓ For n=1, the period-2 equation is a QUARTIC (degree 4)")
    print(f"✓ SymPy can solve it in closed form")
    print(f"✓ Multipliers at period-2 points can be computed exactly")
    print(f"✓ For n=2+, the polynomial degree grows exponentially → infeasible")
else:
    print(f"✗ Could not compute complete analytic solution")
    print(f"  (SymPy may have numerical issues with this parameter choice)")

print(f"\n{'='*80}\n")
