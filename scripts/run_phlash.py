#!/usr/bin/env python3
"""
Run Phlash analysis from a PSMCFA file and export results.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pickle

def main():
    # ==================================================
    # Parse arguments
    # ==================================================
    if len(sys.argv) < 2:
        print("\nUsage: python run_phlash.py <psmcfa_file>")
        print("\nExample: python run_phlash.py simulated_validation.psmcfa\n")
        return
    
    psmcfa_file = sys.argv[1]
    
    if not os.path.exists(psmcfa_file):
        print(f"✗ Error: PSMCFA file not found: {psmcfa_file}")
        return

    print("\n" + "="*70)
    print("PHLASH ANALYSIS")
    print("="*70)
    
    # ==================================================
    # Run Phlash
    # ==================================================
    print(f"\n[Step 1] Running Phlash on {psmcfa_file}...")
    
    phlash_success = False
    phlash_times = None
    phlash_median_ne = None
    phlash_lower = None
    phlash_upper = None
    posterior_samples = None
    
    try:
        import phlash
        print(f"✓ Using Phlash from: {phlash.__file__}")
        if os.environ.get("PHLASH_DEBUG") == "1":
             os.environ["PHLASH_DEBUG_MODE"] = "1"
             print("DEBUG: Phlash debug mode enabled")
        
        # Define mutation rate (used for both inference and rescaling)
        mutation_rate = 1.29e-8
        
        # Use phlash.psmc() convenience function for PSMCFA files
        print("  Running Phlash inference (this may take 1-2 minutes)...")
        posterior_samples = phlash.psmc([psmcfa_file], theta=4 * 10000 * mutation_rate)
        print(f"  ✓ Phlash complete ({len(posterior_samples)} posterior samples)")
        
        # Extract results
        times_all = np.array([dm.eta.t[1:] for dm in posterior_samples])
        phlash_times = np.geomspace(times_all.min(), times_all.max(), 500)
        
        # Rescale with mutation rate
        rescaled_samples = [dm.rescale(mutation_rate) for dm in posterior_samples]
        Nes = np.array([dm.eta(phlash_times, Ne=True) for dm in rescaled_samples])
        
        phlash_median_ne = np.median(Nes, axis=0)
        phlash_lower = np.percentile(Nes, 2.5, axis=0)
        phlash_upper = np.percentile(Nes, 97.5, axis=0)
        
        phlash_success = True
        
    except Exception as e:
        import traceback
        print(f"  ⚠ Phlash failed: {str(e)}")
        print(f"  Full traceback:")
        traceback.print_exc()
        return
    
    # ==================================================
    # Visualization (PSMC-like style)
    # ==================================================
    print("\n[Step 2] Creating visualization...")

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    # psmc_plot.pl uses log-x and linear-y with step lines.
    ax.step(phlash_times, phlash_median_ne / 10000.0, where='post', linewidth=2.2,
            label='Phlash median', color='red')
    ax.fill_between(phlash_times, phlash_lower / 10000.0, phlash_upper / 10000.0,
                    alpha=0.2, color='red', label='95% CI', step='post')
    ax.set_xscale('log')
    ax.set_xlabel('Time (generations ago)', fontsize=11)
    ax.set_ylabel('Effective population size (x10^4)', fontsize=11)
    ax.set_title('Phlash Demographic History', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend(fontsize=10, loc='upper right')
    
    plt.tight_layout()
    
    output_file = 'phlash_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")

    # ==================================================
    # Save execution outputs as pickle
    # ==================================================
    pickle_file = 'phlash_output.pkl'
    output_payload = {
        'input_psmcfa_file': psmcfa_file,
        'mutation_rate': mutation_rate,
        'n_posterior_samples': len(posterior_samples),
        'times_generations': phlash_times,
        'median_ne': phlash_median_ne,
        'lower_95_ne': phlash_lower,
        'upper_95_ne': phlash_upper,
    }
    with open(pickle_file, 'wb') as f:
        pickle.dump(output_payload, f)
    print(f"  ✓ Saved: {pickle_file}")
    
    # ==================================================
    # Summary
    # ==================================================
    print("\n" + "="*70)
    print("PHLASH ANALYSIS COMPLETE")
    print("="*70)
    print("\nResults:")
    if phlash_success:
        print(f"  • Phlash: {len(posterior_samples)} posterior samples")
        print(f"  • Phlash Ne range: {phlash_median_ne.min():.0f} - {phlash_median_ne.max():.0f}")
    
    print("\nGenerated files:")
    print(f"  • {output_file}")
    print(f"  • {pickle_file}")
    print("\n" + "="*70)

if __name__ == '__main__':
    main()

