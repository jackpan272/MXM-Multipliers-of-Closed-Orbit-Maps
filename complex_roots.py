import sympy as sp

# Define the variable
z = sp.symbols('z')

polynomial = z**4 + 1


solutions = sp.solve(polynomial, z)

print("Solutions (symbolic):")
for sol in solutions:
    print(sol)
