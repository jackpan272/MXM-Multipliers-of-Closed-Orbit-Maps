import numpy as np
import sympy as sp

def map_f(z, a):
    """Evaluate f(z) = z * Π_k (z - a_k) / (1 - z * conj(a_k))"""
    result = z
    for ak in a:
        result *= (z - ak) / (1 - z * np.conj(ak))
    return result

def F(z, a):
    """Second iterate minus identity: F(z) = f(f(z)) - z."""
    return map_f(map_f(z, a), a) - z

def newton_complex(g, z0, a, tol=1e-13, max_iter=500):
    """Newton iteration for complex-valued functions using finite-difference derivative."""
    z = z0
    for _ in range(max_iter):
        # numerical derivative g'(z)
        h = 1e-6
        g_z = g(z, a)
        if abs(g_z) < tol:
            return z

        g_prime = (g(z + h, a) - g_z)/h
        if g_prime == 0:
            return None

        z_new = z - g_z/g_prime
        if abs(z_new - z) < tol:
            return z_new
        z = z_new
    return None

def find_period2_points(a, n_samples):
    """Search unit disk for period-2 points."""
    candidates = []

    # uniform random samples in the disk
    for _ in range(n_samples ):  # oversample to account for rejections
        r = np.sqrt(np.random.rand())
        t = 2*np.pi*np.random.rand()
        z0 = r * np.exp(1j * t)

        z_root = newton_complex(F, z0, a)
        if z_root is None:
            continue

        # enforce |z|<1 (domain) and remove diverged values
        if abs(z_root) >= 1:
            continue

        # ensure it's not a fixed point: f(z) != z
        if abs(map_f(z_root, a) - z_root) < 1e-6:
            continue

        candidates.append(z_root)

    # cluster solutions that are numerically identical
    final_points = []
    for z in candidates:
        if not any(abs(z - w) < 1e-4 for w in final_points):
            final_points.append(z)

    return final_points

def perturb_a(a, eps=0.05):
    """
    Perturb each complex value in a by a random complex amount of magnitude < eps,
    then renormalize if necessary to keep |a_k| < 1.
    
    Parameters:
        a   - numpy array of complex numbers with |a_k| < 1
        eps - maximum magnitude of perturbation
        
    Returns:
        numpy array of perturbed complex numbers with |a_k| < 1
    """
    perturbed = []
    for ak in a:
        # random direction + small random magnitude < eps
        r = eps * np.sqrt(np.random.rand())
        theta = 2*np.pi*np.random.rand()
        delta = r * np.exp(1j * theta)

        new_val = ak + delta
        
        # if it leaves the unit disk, project it back inside
        if abs(new_val) >= 1:
            new_val = new_val / abs(new_val) * (1 - 1e-6)

        perturbed.append(new_val)
    
    return np.array(perturbed)

# ------------------------- Example Usage -------------------------
if __name__ == "__main__":
    # example parameter values:
    #a = np.array([0.2 + 0.1j, -0.1 + 0.3j])
    a = np.array([
    0.3 + 0.2j,
    -0.4 + 0.1j,
    0.1 - 0.5j,
    -0.2 - 0.3j
    ])
    n = 3
    #number_of_multipliers = n**2 + 2*n  - (n**2 + n)/2
    pts = find_period2_points(a, n_samples=2000)
    print("Period-2 points:")
    for p in pts:
        print(p)
    print(pts)
    print()
    print()
    norms = abs(np.array(pts))
    print("Norms of found period-2 points:", norms)
    compose_f_twice_vals = [F(z,a) for z in pts]
    print("f(f(z)) values of found period-2 points:", compose_f_twice_vals)

    z, m = sp.symbols('z m')
    #a = sp.symbols(f'a1:{n+1}', complex=True)
    #a_conj = [sp.conjugate(ai) for ai in a]
    #a_bar = sp.symbols(f'a_bar1:{n+1}', complex=True)

    def FF(w, a_par=a):
        expr = w
        for i in range(n):
            expr *= (w - a[i]) / (1 - w * a[i].conjugate())
        return expr

    F2 = FF(FF(z))
    #F2 = F2.ratsimp()
    F2 = sp.cancel(F2)

    F2_prime = sp.diff(F2, z)# get the polynomial
    old_m = [F2_prime.subs(z, zi).evalf() for zi in pts]
    print(old_m)
    print()
    print()

    a_perturbed = perturb_a(a, eps=0.05)
    print("Perturbed a values:", a_perturbed)
    pert_pts = []
    for z0 in pts:  # old period-2 points
        z_new = newton_complex(F, z0, a_perturbed)
        if z_new is None or abs(z_new) >= 1:
            # fallback: keep old value if Newton fails (optional)
            z_new = z0
        pert_pts.append(z_new)

    F2_par = FF(FF(z, a_par=a_perturbed))
    F2_par = sp.cancel(F2_par)
    F2_prime_par = sp.diff(F2_par, z)# get the polynomial
    #pert_pts = find_period2_points(a_perturbed, n_samples=number_of_multipliers)
    new_m = [F2_prime_par.subs(z, zi).evalf() for zi in pert_pts]
    print()
    print()
    print("Comparison of multipliers before and after perturbation:")
    print(np.array(new_m)-np.array(old_m))
    print("length of new ", len(new_m), " length of old ", len(old_m))

