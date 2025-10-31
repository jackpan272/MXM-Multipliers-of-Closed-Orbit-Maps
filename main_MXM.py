# pipeline.py
import sympy as sp
from multiplier_v2 import resultant, a, a_bar  # Import the outputs
from symmetric_poly import sym, E, m # Import the symmetric polynomial function
from Jacobian_Method import compute_jacobian, compute_rank


result_from_symmetric = sym(resultant, m[0])

# Use result in Jacobian_Method
# Convert result to list of expressions for Jacobian computation
functions = result_from_symmetric  # This is already a list
variables = []
variables.extend(a)          # Add a1, a2, ..., an
variables.extend(a_bar)   # Add conjugates

# Compute Jacobian and rank
jacobian = compute_jacobian(functions, variables)
rank = compute_rank(jacobian)

print("Jacobian matrix:")
sp.init_printing(use_unicode=True)

# Or for LaTeX-style output
sp.init_printing(use_latex='mathjax')
sp.pprint(jacobian)
print(f"\nRank: {rank}")