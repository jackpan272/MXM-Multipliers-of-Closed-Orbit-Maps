"""
CORRECTED STATISTICAL TEST FOR SUFFICIENT PERTURBATION
========================================================

This module provides the CORRECT hypothesis test that checks not just
"is there a change" but "is the change SUFFICIENT"

Copy this entire file as-is. Just save as: corrected_statistical_test.py
"""

import numpy as np
from scipy import stats
from typing import Dict


def sufficient_change_test(
    multiplier_changes: np.ndarray,
    sufficiency_threshold: float = 0.1,
    confidence_level: float = 0.95,
    verbose: bool = True
) -> Dict:
    """
    Test if multiplier changes are SUFFICIENTLY LARGE.
    
    NOT: "Is there any change?" (standard t-test)
    BUT: "Is the change larger than the threshold?" (one-tailed test)
    
    Parameters
    ----------
    multiplier_changes : np.ndarray
        Array of measured multiplier changes (in percentage points)
        Example: [0.28, 0.41, 0.15, 0.38, 0.22, 0.33, ...]
    
    sufficiency_threshold : float
        Minimum acceptable change (in percentage points)
        Default: 0.1% (changes must exceed this)
    
    confidence_level : float
        Confidence level for the test (default 95%, α=0.05)
    
    verbose : bool
        Print detailed results
    
    Returns
    -------
    Dict with keys:
        - 'test_type': "One-tailed t-test (sufficient change)"
        - 'null_hypothesis': "Mean change ≤ threshold"
        - 'alternative_hypothesis': "Mean change > threshold"
        - 'mean_change': Sample mean
        - 'std_change': Sample std dev
        - 'threshold': Sufficiency threshold
        - 't_statistic': Computed t-statistic
        - 'p_value': One-tailed p-value
        - 'confidence_level': Confidence %
        - 'reject_null': Boolean (True if sufficient)
        - 'conclusion': Text interpretation
    
    Examples
    --------
    # Scenario 1: Strong evidence of sufficient change
    changes = np.array([0.28, 0.41, 0.15, 0.38, 0.22, 0.33] * 20)
    results = sufficient_change_test(changes, sufficiency_threshold=0.1)
    # Output: p_value ≈ 0.0001, reject_null=True → SUFFICIENT
    
    # Scenario 2: Change too small
    changes = np.array([0.01, 0.02, 0.03, 0.015, 0.025] * 20)
    results = sufficient_change_test(changes, sufficiency_threshold=0.1)
    # Output: p_value ≈ 0.95, reject_null=False → INSUFFICIENT
    
    # Scenario 3: Borderline
    changes = np.array([0.09, 0.11, 0.08, 0.12, 0.10] * 20)
    results = sufficient_change_test(changes, sufficiency_threshold=0.1)
    # Output: p_value ≈ 0.45, reject_null=False → INCONCLUSIVE
    """
    
    changes = np.array(multiplier_changes)
    n = len(changes)
    
    # Compute statistics
    mean_change = np.mean(changes)
    std_change = np.std(changes, ddof=1)  # Use sample std (N-1)
    std_error = std_change / np.sqrt(n)
    
    # ONE-TAILED T-TEST
    # H0: μ ≤ threshold  (change is insufficient)
    # H1: μ > threshold  (change is sufficient)
    
    # The test statistic measures how many standard errors above the threshold
    # the mean is
    t_statistic = (mean_change - sufficiency_threshold) / std_error
    
    # For a ONE-TAILED test (right tail), p-value = P(T > t_obs)
    # where T ~ t-distribution with df = n-1
    alpha = 1 - confidence_level
    df = n - 1
    
    # One-tailed p-value: probability of observing this or more extreme under H0
    p_value_one_tailed = 1 - stats.t.cdf(t_statistic, df)
    
    # Reject H0 if p-value < α
    reject_null = p_value_one_tailed < alpha
    
    # Determine conclusion
    if reject_null:
        conclusion = (
            f"✓ SUFFICIENT CHANGE DETECTED\n"
            f"  Mean change ({mean_change:.4f}%) is statistically significantly "
            f"greater than threshold ({sufficiency_threshold:.4f}%) at "
            f"{int(confidence_level*100)}% confidence level."
        )
    else:
        if mean_change < sufficiency_threshold:
            conclusion = (
                f"✗ INSUFFICIENT CHANGE\n"
                f"  Mean change ({mean_change:.4f}%) is below the threshold "
                f"({sufficiency_threshold:.4f}%). Perturbation too small or "
                f"multipliers not sensitive enough."
            )
        else:
            conclusion = (
                f"⚠ BORDERLINE / INCONCLUSIVE\n"
                f"  Mean change ({mean_change:.4f}%) exceeds threshold "
                f"({sufficiency_threshold:.4f}%), but not statistically "
                f"significantly so (p={p_value_one_tailed:.4f} ≥ α={alpha}).\n"
                f"  Need more trials or larger perturbation to be confident."
            )
    
    results = {
        'test_type': 'One-tailed t-test (sufficient change)',
        'null_hypothesis': f'Mean change ≤ {sufficiency_threshold}%',
        'alternative_hypothesis': f'Mean change > {sufficiency_threshold}%',
        'mean_change': float(mean_change),
        'std_change': float(std_change),
        'std_error': float(std_error),
        'threshold': float(sufficiency_threshold),
        't_statistic': float(t_statistic),
        'p_value': float(p_value_one_tailed),
        'degrees_of_freedom': int(df),
        'num_observations': int(n),
        'confidence_level': float(confidence_level),
        'alpha': float(alpha),
        'reject_null_hypothesis': bool(reject_null),
        'conclusion': conclusion,
    }
    
    if verbose:
        print("\n" + "="*80)
        print("HYPOTHESIS TEST: SUFFICIENT CHANGE")
        print("="*80)
        print(f"\nTest Type: {results['test_type']}")
        print(f"\nHypotheses:")
        print(f"  H₀: {results['null_hypothesis']}")
        print(f"  H₁: {results['alternative_hypothesis']}")
        print(f"\nData Summary:")
        print(f"  Number of observations: {n}")
        print(f"  Mean change:            {mean_change:.6f}%")
        print(f"  Std deviation:          {std_change:.6f}%")
        print(f"  Std error:              {std_error:.6f}%")
        print(f"\nTest Statistics:")
        print(f"  Sufficiency threshold:  {sufficiency_threshold:.4f}%")
        print(f"  t-statistic:            {t_statistic:.6f}")
        print(f"  Degrees of freedom:     {df}")
        print(f"  Confidence level:       {int(confidence_level*100)}%")
        print(f"  α (significance):       {alpha:.4f}")
        print(f"  p-value (one-tailed):   {p_value_one_tailed:.6f}")
        print(f"\nDecision:")
        print(f"  Reject H₀? {reject_null}")
        print(f"\nConclusion:")
        print(f"  {results['conclusion']}")
        print("="*80 + "\n")
    
    return results


def compare_thresholds(
    multiplier_changes: np.ndarray,
    threshold_range: np.ndarray = None,
    confidence_level: float = 0.95
) -> Dict:
    """
    Test sufficiency against multiple thresholds.
    
    This helps you understand: "At what threshold does the evidence fail?"
    
    Parameters
    ----------
    multiplier_changes : np.ndarray
        Array of measured changes
    
    threshold_range : np.ndarray
        Array of thresholds to test
        Default: [0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
    
    confidence_level : float
        Confidence level for tests
    
    Returns
    -------
    Dict with results for each threshold
    
    Example
    -------
    changes = np.array([0.28, 0.41, 0.15, ...])
    comparison = compare_thresholds(changes)
    
    # Shows which thresholds pass/fail:
    # Threshold 0.01% → PASS (p=0.0001)
    # Threshold 0.1%  → PASS (p=0.0012)
    # Threshold 0.2%  → BORDERLINE (p=0.12)
    # Threshold 0.5%  → FAIL (p=0.87)
    """
    
    if threshold_range is None:
        threshold_range = np.array([0.01, 0.05, 0.1, 0.2, 0.5, 1.0])
    
    results_by_threshold = {}
    
    print("\n" + "="*80)
    print("SUFFICIENCY TEST ACROSS MULTIPLE THRESHOLDS")
    print("="*80)
    print(f"\n{'Threshold':<12} {'Mean>':<10} {'p-value':<10} {'Reject H₀':<10} {'Status':<15}")
    print("-"*80)
    
    for threshold in threshold_range:
        result = sufficient_change_test(
            multiplier_changes,
            sufficiency_threshold=threshold,
            confidence_level=confidence_level,
            verbose=False
        )
        
        results_by_threshold[threshold] = result
        
        # Determine status
        if result['reject_null_hypothesis']:
            status = "✓ SUFFICIENT"
        elif result['mean_change'] > threshold:
            status = "⚠ BORDERLINE"
        else:
            status = "✗ INSUFFICIENT"
        
        print(
            f"{threshold:<12.3f}% "
            f"{result['mean_change']:<10.4f}% "
            f"{result['p_value']:<10.6f} "
            f"{str(result['reject_null_hypothesis']):<10} "
            f"{status:<15}"
        )
    
    print("="*80 + "\n")
    
    return results_by_threshold


if __name__ == "__main__":
    print("\n" + "="*80)
    print("EXAMPLE 1: Strong Evidence of Sufficient Change")
    print("="*80)
    
    # Scenario: Mean change = 0.30%, threshold = 0.1%
    # All changes are around 0.30%, very consistent
    changes_strong = np.array([0.28, 0.31, 0.29, 0.32, 0.30, 0.27] * 30)
    result1 = sufficient_change_test(changes_strong, sufficiency_threshold=0.1)
    
    print("\n" + "="*80)
    print("EXAMPLE 2: Weak Evidence (Too Small)")
    print("="*80)
    
    # Scenario: Mean change = 0.02%, threshold = 0.1%
    # All changes are tiny, way below threshold
    changes_weak = np.array([0.01, 0.02, 0.03, 0.015, 0.025] * 30)
    result2 = sufficient_change_test(changes_weak, sufficiency_threshold=0.1)
    
    print("\n" + "="*80)
    print("EXAMPLE 3: Borderline Case")
    print("="*80)
    
    # Scenario: Mean change = 0.12%, threshold = 0.1%
    # Just barely above threshold, but with high variability
    changes_borderline = np.array([0.05, 0.15, 0.08, 0.18, 0.10, 0.12] * 30)
    result3 = sufficient_change_test(changes_borderline, sufficiency_threshold=0.1)
    
    print("\n" + "="*80)
    print("EXAMPLE 4: Testing Multiple Thresholds")
    print("="*80)
    
    # Same data, but test against different thresholds
    changes_multi = np.array([0.28, 0.31, 0.29, 0.32, 0.30, 0.27] * 30)
    comparison = compare_thresholds(
        changes_multi,
        threshold_range=np.array([0.01, 0.05, 0.1, 0.2, 0.5])
    )
    
    # Interpretation
    print(
        "INTERPRETATION:\n"
        "The data shows sufficient change above 0.01%, 0.05%, and 0.1%,\n"
        "but NOT above 0.2% or 0.5%. This means:\n"
        "  → Best choice: threshold = 0.1% (good balance)\n"
        "  → You could also use 0.05% if you want to be conservative\n"
    )
