"""
NUMERICAL PERTURBATION ANALYSIS - COMPLETE FILE
===============================================

This is the full implementation. Copy this entire file as-is.
No modifications needed. Just save as: numerical_perturbation.py
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict, Callable
import warnings

warnings.filterwarnings('ignore')


class NumericalMultiplierAnalysis:
    """
    Analyze multipliers of period-2 points under parameter perturbations.
    """
    
    def __init__(self, n: int = 2, verbose: bool = False):
        """Initialize analysis."""
        self.n = n
        self.verbose = verbose
        self.num_period2_points = 2**(n+1) - 2
    
    def map_f(self, z: complex, a: np.ndarray) -> complex:
        """Evaluate f(z) = z * ∏_k (z - a_k) / (1 - z * conj(a_k))"""
        result = z
        for k in range(len(a)):
            result *= (z - a[k]) / (1 - z * np.conj(a[k]))
        return result
    
    def map_f_derivative(self, z: complex, a: np.ndarray, 
                        step: float = 1e-8) -> complex:
        """Numerically compute f'(z) using finite differences."""
        dz = complex(step, step)
        f_plus = self.map_f(z + dz, a)
        f_minus = self.map_f(z, a)
        return (f_plus - f_minus) / (dz)
    
    def compose_f_twice(self, z: complex, a: np.ndarray) -> complex:
        """Compute f(f(z))"""
        return self.map_f(self.map_f(z, a), a)
    
    def compose_f_twice_derivative(self, z: complex, a: np.ndarray,
                                step: float = 1e-8) -> complex:
  
        fz = self.map_f(z, a)                    # Compute f(z)
        f_prime_z = self.map_f_derivative(z, a)  # Compute f'(z)
        f_prime_fz = self.map_f_derivative(fz, a)  # Compute f'(f(z))
    
    
        return f_prime_fz * f_prime_z


    
    def find_period2_points(self, a: np.ndarray, 
                           num_guesses: int = 500) -> List[complex]:
        """Find period-2 points numerically."""
        period2_points = []
        tolerance = 1e-12
        
        for _ in range(num_guesses):
            theta = 2 * np.pi * np.random.random()
            z0 = complex(np.cos(theta), np.sin(theta))
            
            try:
                def residual_real(x):
                    z = complex(x[0], x[1])
                    res = self.compose_f_twice(z, a) - z
                    return [res.real, res.imag]
                
                sol = fsolve(residual_real, [z0.real, z0.imag], full_output=True)
                z_found = complex(sol[0][0], sol[0][1])
                info = sol[1]
                
                residual = np.linalg.norm(info['fvec'])
                if residual > tolerance:
                    continue
                
                fz = self.map_f(z_found, a)
                if abs(fz - z_found) < tolerance:
                    continue
                
                is_new = True
                for existing in period2_points:
                    if abs(z_found - existing) < tolerance:
                        is_new = False
                        break
                
                if is_new:
                    z_found = z_found / abs(z_found)
                    period2_points.append(z_found)
            
            except Exception:
                continue
 
        return period2_points
    
    def compute_multipliers(self, period2_points: List[complex],
                           a: np.ndarray) -> np.ndarray:
        """Compute multipliers m_j = (f∘f)'(z_j)"""
        multipliers = np.array([
            self.compose_f_twice_derivative(z, a) 
            for z in period2_points
        ])
        return multipliers


 
    

    

    
    def perturb_parameter(self, a: np.ndarray, 
                         perturbation_magnitude: float) -> np.ndarray:
        """Perturb parameters by small random amount."""
        a_perturbed = a.copy()
        for k in range(len(a)):
            delta = perturbation_magnitude * abs(a[k]) * np.exp(
                2j * np.pi * np.random.random()
            )
            a_perturbed[k] = a[k] + delta
            
            if abs(a_perturbed[k]) >= 0.99:
                a_perturbed[k] = 0.99 * a_perturbed[k] / abs(a_perturbed[k])
        
        return a_perturbed
    
    def measure_multiplier_sensitivity(self, a: np.ndarray,
                                      perturbation_magnitude: float = 0.01,
                                      num_trials: int = 100) -> Dict:
        """Measure how multipliers change under perturbations."""
        
        if self.verbose:
            print(f"\n{'='*70}")
            print(f"MULTIPLIER SENSITIVITY ANALYSIS")
            print(f"{'='*70}")
            print(f"Finding period-2 points for base parameters...")
        
        period2_points = self.find_period2_points(a, num_guesses=100)
        
        if not period2_points:
            print("ERROR: Could not find period-2 points!")
            return None
        
        if self.verbose:
            print(f"Found {len(period2_points)} period-2 points")
        
        base_multipliers = self.compute_multipliers(period2_points, a)
        
        if self.verbose:
            print(f"\nBase multipliers:")
            for i, m in enumerate(base_multipliers):
                print(f"  m_{i+1} = {m:.6f}")
        
        all_changes = []
        perturbed_multiplier_sets = []
        changed_count = 0
        
        for trial in range(num_trials):
            a_perturbed = self.perturb_parameter(a, perturbation_magnitude)
            
            try:
                perturbed_multipliers = self.compute_multipliers(
                    period2_points, a_perturbed
                )
                perturbed_multiplier_sets.append(perturbed_multipliers)
                
                changes = []
                for m_base, m_pert in zip(base_multipliers, perturbed_multipliers):
                    if abs(m_base) > 1e-10:
                        relative_change = 100 * abs(m_pert - m_base) / abs(m_base)
                        changes.append(relative_change)
                    else:
                        changes.append(100 * abs(m_pert))
                
                all_changes.extend(changes)
                
                if any(c > 0.1 for c in changes):
                    changed_count += 1
            
            except Exception as e:
                continue
        
        all_changes = np.array(all_changes)
        
        results = {
            'multiplier_changes': all_changes.tolist(),
            'num_changed': changed_count,
            'num_trials': num_trials,
            'num_period2_points': len(period2_points),
            'mean_change': float(np.mean(all_changes)) if len(all_changes) > 0 else 0.0,
            'std_change': float(np.std(all_changes)) if len(all_changes) > 0 else 0.0,
            'max_change': float(np.max(all_changes)) if len(all_changes) > 0 else 0.0,
            'min_change': float(np.min(all_changes)) if len(all_changes) > 0 else 0.0,
            'median_change': float(np.median(all_changes)) if len(all_changes) > 0 else 0.0,
            'base_multipliers': base_multipliers,
            'perturbed_multipliers': perturbed_multiplier_sets,
            'perturbation_magnitude': perturbation_magnitude,
        }
        
        if self.verbose:
            print(f"\n{'─'*70}")
            print(f"RESULTS (Perturbation: {perturbation_magnitude*100:.2f}%)")
            print(f"{'─'*70}")
            print(f"Trials completed: {num_trials}")
            print(f"Multiplier changes (relative %):")
            print(f"  Mean:   {results['mean_change']:.4f}%")
            print(f"  Median: {results['median_change']:.4f}%")
            print(f"  Std:    {results['std_change']:.4f}%")
            print(f"  Min:    {results['min_change']:.4f}%")
            print(f"  Max:    {results['max_change']:.4f}%")
            print(f"  Trials with >0.1% change: {changed_count}/{num_trials}")
            print(f"{'─'*70}\n")
        
        return results
    
    def sweep_perturbation_magnitudes(self, a: np.ndarray,
                                     perturbation_range: np.ndarray,
                                     num_trials_per_magnitude: int = 50) -> Dict:
        """Sweep over different perturbation magnitudes."""
        sweep_results = {}
        
        print(f"\n{'='*70}")
        print(f"PERTURBATION MAGNITUDE SWEEP")
        print(f"{'='*70}")
        
        for mag in perturbation_range:
            print(f"\nTesting perturbation magnitude: {mag*100:.3f}%")
            
            results = self.measure_multiplier_sensitivity(
                a, 
                perturbation_magnitude=mag,
                num_trials=num_trials_per_magnitude
            )
            
            sweep_results[mag] = results
        
        return sweep_results
    
    def statistical_test(self, a: np.ndarray, 
                        perturbation_magnitude: float = 0.01,
                        num_trials: int = 100,
                        significance_level: float = 0.05) -> Dict:
        """Statistical hypothesis test."""
        from scipy import stats
        
        results = self.measure_multiplier_sensitivity(
            a, 
            perturbation_magnitude=perturbation_magnitude,
            num_trials=num_trials
        )
        
        changes = np.array(results['multiplier_changes'])
        
        t_statistic, p_value = stats.ttest_1samp(changes, 0.0)
        
        reject_h0 = p_value < significance_level
        
        test_results = {
            'null_hypothesis': 'Multipliers do NOT change under perturbation',
            'alternative_hypothesis': 'Multipliers DO change under perturbation',
            't_statistic': float(t_statistic),
            'p_value': float(p_value),
            'significance_level': significance_level,
            'reject_null_hypothesis': bool(reject_h0),
            'conclusion': 'Multipliers ARE sensitive' if reject_h0 else 'Cannot confirm sensitivity',
            'mean_change_percent': results['mean_change'],
            'num_observations': len(changes),
        }
        
        print(f"\n{'='*70}")
        print(f"STATISTICAL HYPOTHESIS TEST")
        print(f"{'='*70}")
        print(f"H0: {test_results['null_hypothesis']}")
        print(f"H1: {test_results['alternative_hypothesis']}")
        print(f"\nResults:")
        print(f"  t-statistic: {t_statistic:.4f}")
        print(f"  p-value:     {p_value:.6f}")
        print(f"  α (sig level): {significance_level}")
        print(f"  Mean change: {results['mean_change']:.4f}%")
        print(f"\nConclusion: {test_results['conclusion']}")
        print(f"Reject H0: {reject_h0}")
        print(f"{'='*70}\n")
        
        return test_results

   