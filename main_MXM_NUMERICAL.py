
import numpy as np
import sympy as sp
from typing import Dict, List, Tuple
import pickle

# ============================================================================
# PHASE 1: NUMERICAL ANALYSIS 
# ============================================================================

from numerical_perturbation import NumericalMultiplierAnalysis
from corrected_statistical_test import sufficient_change_test, compare_thresholds


def phase_1_numerical_sensitivity(n: int = 2, num_trials: int = 100) -> Dict:
    """
    Parameters
    ----------
    n : int
        Number of parameters (default n=2)
    num_trials : int
        Number of perturbation trials
    
    Returns
    -------
    Dict with multiplier sensitivities
    """
    
    print("\n" + "="*80)
    print("PHASE 1: NUMERICAL PERTURBATION ANALYSIS")
    print("="*80)
    
    # Initialize analyzer
    analyzer = NumericalMultiplierAnalysis(n=n, verbose=True)
    
    # Generate random parameters in unit disk (safe defaults)
    np.random.seed(42)
    a = 0.7 * np.exp(2j * np.pi * np.random.random(n))
    
    print(f"\nGenerated {n} random parameters:")
    for k, ak in enumerate(a):
        print(f"  a_{k+1} = {ak:.6f}  (|a_{k+1}| = {abs(ak):.6f})")
    
    # ─────────────────────────────────────────────────────────────────────
    # Step 1a: Measure sensitivity at 1% perturbation
    # ─────────────────────────────────────────────────────────────────────
    
    results_1pct = analyzer.measure_multiplier_sensitivity(
        a=a,
        perturbation_magnitude=0.01,  # 1% perturbation
        num_trials=num_trials
    )
    
    print(f"\n{'─'*80}")
    print("SENSITIVITY AT 1% PERTURBATION:")
    print(f"{'─'*80}")
    print(f"Mean change:   {results_1pct['mean_change']:.6f}%")
    print(f"Std dev:       {results_1pct['std_change']:.6f}%")
    print(f"Min/Max:       {results_1pct['min_change']:.6f}% / {results_1pct['max_change']:.6f}%")
    
    # ─────────────────────────────────────────────────────────────────────
    # Step 1b: Statistical test with CORRECTED hypothesis
    # ─────────────────────────────────────────────────────────────────────
    
    print(f"\n{'─'*80}")
    print("HYPOTHESIS TEST (Corrected):")
    print(f"{'─'*80}")
    
    test_results = sufficient_change_test(
        np.array(results_1pct['multiplier_changes']),
        sufficiency_threshold=0.1,  # YOUR threshold
        confidence_level=0.95,
        verbose=True
    )
    
    # ─────────────────────────────────────────────────────────────────────
    # Step 1c: Sweep different perturbation magnitudes
    # ─────────────────────────────────────────────────────────────────────
    
    print(f"\n{'─'*80}")
    print("PERTURBATION MAGNITUDE SWEEP:")
    print(f"{'─'*80}")
    
    sweep_results = analyzer.sweep_perturbation_magnitudes(
        a=a,
        perturbation_range=np.array([0.001, 0.005, 0.01, 0.05]),
        num_trials_per_magnitude=50
    )
    
    # ─────────────────────────────────────────────────────────────────────
    # PHASE 1 CONCLUSION
    # ─────────────────────────────────────────────────────────────────────
    
    phase1_result = {
        'parameters': a,
        'period2_points': analyzer.find_period2_points(a, num_guesses=100),
        'multipliers': results_1pct['base_multipliers'],
        'sensitivity_1pct': results_1pct,
        'statistical_test': test_results,
        'sweep_results': sweep_results,
        'conclusion': (
            "✓ SUFFICIENT" if test_results['reject_null_hypothesis']
            else "✗ INSUFFICIENT"
        )
    }
    
    print(f"\n{'─'*80}")
    print("PHASE 1 CONCLUSION:")
    print(f"{'─'*80}")
    print(f"Multipliers ARE sensitive: {phase1_result['conclusion']}")
    print(f"Statistical significance: p-value = {test_results['p_value']:.6f}")
    
    return phase1_result


# ============================================================================
# PHASE 2: JACOBIAN RANK ANALYSIS (YOUR EXISTING CODE, MODIFIED)
# ============================================================================

def phase_2_jacobian_rank(
    phase1_result: Dict,
    n: int = 2,
    num_sample_points: int = 10
) -> Dict:
    """
    Phase 2: Theoretical validation via Jacobian rank analysis.
    
    Uses finite-difference Jacobian for any n.
    
    Parameters
    ----------
    phase1_result : Dict
        Output from Phase 1
    n : int
        Number of parameters
    num_sample_points : int
        How many random parameter sets to test
    
    Returns
    -------
    Dict with Jacobian analysis results
    """
    
    print("\n" + "="*80)
    print("PHASE 2: JACOBIAN RANK ANALYSIS")
    print("="*80)
    
    print("\nUsing finite-difference Jacobian...")
    phase2_result = _phase2_finite_diff(phase1_result, n)
    
    return phase2_result


def _phase2_finite_diff(phase1_result: Dict, n: int, epsilon: float = 1e-6) -> Dict:
    """
    Jacobian analysis using finite differences (fast, for any n).
    
    Numerically estimates ∂m_j/∂a_k using finite differences.
    """
    
    print(f"\nEstimating Jacobian via finite differences (ε={epsilon})...")
    
    parameters = phase1_result['parameters']
    period2_points = phase1_result['period2_points']
    
    if not period2_points:
        print("⚠ No period-2 points found. Cannot compute Jacobian.")
        return None
    
    # ─────────────────────────────────────────────────────────────────────
    # Build finite-difference Jacobian
    # ─────────────────────────────────────────────────────────────────────
    
    analyzer = NumericalMultiplierAnalysis(n=n, verbose=False)
    
    # Base multipliers
    base_mults = analyzer.compute_multipliers(period2_points, parameters)
    num_multipliers = len(base_mults)
    
    # Jacobian will be (num_multipliers) × (2n) matrix
    # (since each complex a_k counts as 2 real parameters: Re, Im)
    
    jacobian_fd = np.zeros((num_multipliers, 2*n), dtype=complex)
    
    for k in range(n):
        # Perturb Re(a_k)
        params_perturb = parameters.copy()
        params_perturb[k] = parameters[k] + epsilon
        mults_perturb_re = analyzer.compute_multipliers(period2_points, params_perturb)
        
        # Perturb Im(a_k)
        params_perturb = parameters.copy()
        params_perturb[k] = parameters[k] + 1j*epsilon
        mults_perturb_im = analyzer.compute_multipliers(period2_points, params_perturb)
        
        # Finite differences: ∂m/∂(Re a_k) and ∂m/∂(Im a_k)
        jacobian_fd[:, 2*k] = (mults_perturb_re - base_mults) / epsilon
        jacobian_fd[:, 2*k + 1] = (mults_perturb_im - base_mults) / epsilon
    
    # ─────────────────────────────────────────────────────────────────────
    # Analyze Jacobian
    # ─────────────────────────────────────────────────────────────────────
    
    # SVD to get rank and condition number
    U, singular_values, Vt = np.linalg.svd(jacobian_fd, full_matrices=False)
    
    rank = np.sum(singular_values > 1e-10)
    expected_rank = min(num_multipliers, 2*n)
    
    if singular_values[-1] > 1e-10:
        condition_number = singular_values[0] / singular_values[-1]
    else:
        condition_number = np.inf
    
    # ─────────────────────────────────────────────────────────────────────
    # Results
    # ─────────────────────────────────────────────────────────────────────
    
    phase2_result = {
        'method': 'finite_difference',
        'n': n,
        'epsilon': epsilon,
        'jacobian': jacobian_fd,
        'singular_values': singular_values.tolist(),
        'rank': int(rank),
        'expected_rank': int(expected_rank),
        'condition_number': float(condition_number),
        'is_full_rank': (rank == expected_rank),
    }
    
    print(f"\n{'─'*80}")
    print("JACOBIAN ANALYSIS RESULTS:")
    print(f"{'─'*80}")
    print(f"Jacobian shape:  {jacobian_fd.shape}")
    print(f"Rank:            {rank}")
    print(f"Expected (full): {expected_rank}")
    print(f"Full rank?       {phase2_result['is_full_rank']}")
    print(f"Condition number: {condition_number:.2e}")
    print(f"Smallest singular value: {singular_values[-1]:.2e}")
    
    return phase2_result


# ============================================================================
# PHASE 3: COMBINED CONCLUSION
# ============================================================================

def phase_3_conclusion(phase1_result: Dict, phase2_result: Dict) -> str:
    """
    Phase 3: Synthesize evidence from both approaches.
    
    Returns a publication-ready conclusion.
    """
    
    print("\n" + "="*80)
    print("PHASE 3: SYNTHESIS & CONCLUSION")
    print("="*80)
    
    conclusion_parts = []
    
    # From Phase 1
    if phase1_result['conclusion'] == "✓ SUFFICIENT":
        conclusion_parts.append(
            f"✓ Numerical perturbation test: Multipliers change by "
            f"{phase1_result['sensitivity_1pct']['mean_change']:.4f}% when "
            f"parameters perturbed by 1% (p={phase1_result['statistical_test']['p_value']:.6f})"
        )
    else:
        conclusion_parts.append(
            f"⚠ Numerical perturbation test INCONCLUSIVE (p={phase1_result['statistical_test']['p_value']:.6f})"
        )
    
    # From Phase 2
    if phase2_result is not None:
        if phase2_result.get('is_full_rank'):
            conclusion_parts.append(
                f"✓ Jacobian rank analysis: Full rank ({phase2_result['rank']}/{phase2_result['expected_rank']}) "
                f"with condition number {phase2_result['condition_number']:.2e}"
            )
        else:
            conclusion_parts.append(
                f"⚠ Jacobian rank analysis: Not full rank ({phase2_result['rank']}/{phase2_result['expected_rank']})"
            )
    
    # Final synthesis
    print("\nFINAL EVIDENCE:")
    for i, part in enumerate(conclusion_parts, 1):
        print(f"  {i}. {part}")
    
    if all("✓" in part for part in conclusion_parts):
        final = "✓✓ STRONG EVIDENCE: Multipliers are sensitive, independent coordinates"
    else:
        final = "⚠ MIXED EVIDENCE: Examine individual results above"
    
    print(f"\n{final}")
    
    return final


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    
    print("\n" + "="*80)
    print("MXM RESEARCH: NUMERICAL + JACOBIAN HYBRID APPROACH")
    print("="*80)
    
    # Parameters
    n = 10  # Number of parameters
    num_perturbation_trials = 100
    num_jacobian_samples = 10
    
    # ────────────────────────────────────────────────────────────────────
    # PHASE 1: Quick numerical evidence
    # ────────────────────────────────────────────────────────────────────
    
    phase1 = phase_1_numerical_sensitivity(n=n, num_trials=num_perturbation_trials)
    
    # ────────────────────────────────────────────────────────────────────
    # PHASE 2: Theoretical Jacobian validation
    # ────────────────────────────────────────────────────────────────────
    
    phase2 = phase_2_jacobian_rank(phase1, n=n, num_sample_points=num_jacobian_samples)
    
    # ────────────────────────────────────────────────────────────────────
    # PHASE 3: Synthesize conclusion
    # ────────────────────────────────────────────────────────────────────
    
    conclusion = phase_3_conclusion(phase1, phase2)
    
    # ────────────────────────────────────────────────────────────────────
    # SAVE RESULTS
    # ────────────────────────────────────────────────────────────────────
    
    all_results = {
        'phase1': phase1,
        'phase2': phase2,
        'conclusion': conclusion,
        'parameters': {'n': n, 'trials': num_perturbation_trials}
    }
    
    filename = f'mxm_analysis_n{n}.pkl'
    with open(filename, 'wb') as f:
        pickle.dump(all_results, f)
    
    print(f"\nResults saved to: {filename}")
    print("\n" + "="*80)
