

"""
MXM Jacobian Scaling Analysis with Grid of Entry Norms (FIXED v3)
================================================================
Runs the numerical Jacobian analysis for n = 2, 3, 4, ...
and generates:
1. Jacobian Rank vs n (line plot)
2. 5×5 Grid of Entry Norms vs n (25 subplots showing evolution of each entry)

WITH DIAGNOSTIC OUTPUT to debug the rank issue
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Headless mode for servers
import pickle

# Import your existing modules
from main_MXM_NUMERICAL import phase_1_numerical_sensitivity, phase_2_jacobian_rank

def run_scaling_analysis(max_n=8):
    print("="*80)
    print(f"JACOBIAN SCALING ANALYSIS (n = 2 to {max_n})")
    print("="*80)

    results_by_n = {}

    # Range of n to test
    n_values = range(2, max_n + 1)

    for n in n_values:
        print(f"\n{'='*80}")
        print(f"Analyzing n = {n}...")
        print(f"{'='*80}")

        # 1. Run Phase 1 to get period-2 points
        phase1 = phase_1_numerical_sensitivity(n=n, num_trials=10)

        # 2. Run Phase 2 to get Jacobian
        if not phase1['period2_points']:
            print(f"⚠ Skipping n={n}: No period-2 points found.")
            continue

        phase2 = phase_2_jacobian_rank(phase1, n=n, num_sample_points=1)

        if phase2:
            # ✓ ADD DIAGNOSTIC OUTPUT HERE
            jacobian = phase2['jacobian']
            num_multipliers = phase2['jacobian'].shape[0]

            print(f"\n[DIAGNOSTIC OUTPUT for n={n}]")
            print(f"  Jacobian shape: {jacobian.shape} (rows × cols)")
            print(f"  Number of multipliers (rows): {num_multipliers}")
            print(f"  Expected (2^(n+1) - 2): {2**(n+1) - 2}")
            print(f"  Number of parameter dimensions (cols): {jacobian.shape[1]}")
            print(f"  Expected (2*n): {2*n}")
            print(f"  Computed rank: {phase2['rank']}")
            print(f"  Expected full rank: {phase2['expected_rank']}")
            print(f"  Condition number: {phase2['condition_number']:.2e}")
            print(f"  Singular values (first 5): {phase2['singular_values'][:5]}")

            results_by_n[n] = phase2

        else:
            print(f"⚠ Phase 2 returned None for n={n}")

    return results_by_n

def plot_results(results_by_n, max_n=8):
    print(f"\n{'='*80}")
    print("Creating plots...")
    print(f"{'='*80}")

    n_values = sorted(results_by_n.keys())

    if not n_values:
        print("No results to plot.")
        return

    # ─────────────────────────────────────────────────────────────────────
    # FIGURE 1: Rank vs n (line plot)
    # ─────────────────────────────────────────────────────────────────────
    fig1, ax1 = plt.subplots(1, 1, figsize=(10, 6))

    ranks = []
    expected_ranks = []
    jacobian_shapes = []

    for n in n_values:
        res = results_by_n[n]
        jacobian = res['jacobian']
        ranks.append(res['rank'])
        expected_ranks.append(res['expected_rank'])
        jacobian_shapes.append(f"{jacobian.shape[0]}×{jacobian.shape[1]}")

    print("\nJacobian shapes by n:")
    for n, shape in zip(n_values, jacobian_shapes):
        print(f"  n={n}: {shape}")

    # Plot Rank vs n
    ax1.plot(n_values, ranks, 'o-', color='blue', linewidth=2.5, markersize=8, label='Computed Rank')
    ax1.plot(n_values, expected_ranks, 's--', color='gray', alpha=0.6, linewidth=2, label='Expected (Full Rank)')
    ax1.set_xlabel('n (Number of Parameters)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Rank', fontsize=12, fontweight='bold')
    ax1.set_title('Jacobian Rank vs n', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(n_values)

    plt.tight_layout()
    plt.savefig('jacobian_rank_vs_n.png', dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved rank plot to jacobian_rank_vs_n.png")
    plt.close()

    # ─────────────────────────────────────────────────────────────────────
    # FIGURE 2: 5×5 Grid of Entry Norms vs n
    # ─────────────────────────────────────────────────────────────────────
    print("\nGenerating 5×5 grid of entry norms...")

    # FIX v3: Use tight_layout=False and constrained_layout=False to have full control
    fig2, axes = plt.subplots(5, 5, figsize=(18, 19))
    
    # Collect entry values for each (row, col) across all n
    entry_data = {}  # (row, col) -> list of magnitudes across n values

    for n in n_values:
        res = results_by_n[n]
        jacobian = res['jacobian']

        # Extract entries (0-4, 0-4) - the 5×5 upper left corner
        for row in range(min(5, jacobian.shape[0])):
            for col in range(min(5, jacobian.shape[1])):
                key = (row, col)
                if key not in entry_data:
                    entry_data[key] = []
                entry_data[key].append(abs(jacobian[row, col]))

    # Plot each entry
    for row in range(5):
        for col in range(5):
            ax = axes[row, col]
            key = (row, col)

            if key in entry_data and len(entry_data[key]) > 0:
                magnitudes = entry_data[key]
                # Pad with NaN if not enough values (for smaller matrices)
                while len(magnitudes) < len(n_values):
                    magnitudes.insert(0, np.nan)
                
                ax.plot(n_values, magnitudes, 'o-', color='darkblue', linewidth=2, markersize=6)
                ax.grid(True, alpha=0.3)
                ax.set_title(f'Entry ({row+1},{col+1})', fontsize=10, fontweight='bold')
                ax.set_xticks(n_values)
                ax.set_ylabel('|∂m/∂a|', fontsize=9)
                
                if row == 4:  # Bottom row
                    ax.set_xlabel('n', fontsize=9)
            else:
                ax.text(0.5, 0.5, f'Entry ({row+1},{col+1})\nN/A', 
                       ha='center', va='center', fontsize=10, color='gray')
                ax.set_xticks([])
                ax.set_yticks([])

    # FIX v3: Use subplots_adjust to manually control spacing
    # This gives the title plenty of room (top=0.94) and adjusts hspace/wspace for subplots
    fig2.subplots_adjust(top=0.94, hspace=0.35, wspace=0.30)
    
    fig2.suptitle('Jacobian Entry Norms vs n (5×5 Grid, Entry (1,1) to (5,5))', 
                  fontsize=16, fontweight='bold', y=0.98)
    
    plt.savefig('jacobian_entry_grid_vs_n.png', dpi=300, bbox_inches='tight')
    print(f"✓ Saved entry grid plot to jacobian_entry_grid_vs_n.png")
    plt.close()

if __name__ == "__main__":
    # Run analysis up to n=8
    data = run_scaling_analysis(max_n=8)

    # Save raw data
    with open('scaling_data.pkl', 'wb') as f:
        pickle.dump(data, f)
    print(f"\n✓ Saved raw data to scaling_data.pkl")

    # Plot
    plot_results(data, max_n=8)

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print("\nGenerated files:")
    print("  1. jacobian_rank_vs_n.png        (rank line plot)")
    print("  2. jacobian_entry_grid_vs_n.png  (5×5 grid of entry evolution)")
    print("  3. scaling_data.pkl              (raw data)")
