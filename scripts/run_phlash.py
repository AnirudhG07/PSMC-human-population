#!/usr/bin/env python3
"""
Run Phlash analysis: Takes PSMCFA file and compares with PSMC results
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

def parse_psmc(filename):
    """Parse PSMC output file"""
    times, lambdas = [], []
    in_final_rd = False
    with open(filename) as f:
        for line in f:
            if line.startswith('RD'):
                in_final_rd = True
            elif line.startswith('RS') and in_final_rd:
                parts = line.split()
                times.append(float(parts[2]))
                lambdas.append(float(parts[3]))
    return np.array(times), np.array(lambdas)

def main():
    # ==================================================
    # Parse arguments
    # ==================================================
    if len(sys.argv) < 2:
        print("\nUsage: python run_phlash.py <psmcfa_file> [psmc_file]")
        print("\nExample: python run_phlash.py simulated_validation.psmcfa simulated_validation.psmc")
        print("\nIf psmc_file is not provided, comparison plot will be skipped.\n")
        return
    
    psmcfa_file = sys.argv[1]
    psmc_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.exists(psmcfa_file):
        print(f"✗ Error: PSMCFA file not found: {psmcfa_file}")
        return
    
    if psmc_file and not os.path.exists(psmc_file):
        print(f"✗ Error: PSMC file not found: {psmc_file}")
        return

    print("\n" + "="*70)
    print("PHLASH ANALYSIS")
    print("="*70)
    
    # ==================================================
    # Parse PSMC results (if provided)
    # ==================================================
    psmc_times = None
    psmc_lambdas = None
    
    if psmc_file:
        print(f"\n[Step 1] Parsing PSMC results from {psmc_file}...")
        psmc_times, psmc_lambdas = parse_psmc(psmc_file)
        print(f"  ✓ Extracted {len(psmc_times)} time points from PSMC")
    
    # ==================================================
    # Run Phlash
    # ==================================================
    print(f"\n[Step 2] Running Phlash on {psmcfa_file}...")
    
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
    # Visualization
    # ==================================================
    print("\n[Step 3] Creating visualization...")
    
    if psmc_file and psmc_times is not None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    else:
        fig, ax2 = plt.subplots(1, 1, figsize=(10, 6))
        ax1 = None
    
    # Plot 1: PSMC (if available)
    if ax1 is not None:
        ax1.plot(psmc_times, psmc_lambdas, 'o-', linewidth=2, markersize=5,
                 color='blue', label='PSMC', alpha=0.8)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel('Scaled time (2N₀ generations)', fontsize=11)
        ax1.set_ylabel('Scaled population size (λ)', fontsize=11)
        ax1.set_title('PSMC Results', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3, which='both')
        ax1.legend(fontsize=10)
        from matplotlib.ticker import ScalarFormatter
        ax1.yaxis.set_major_formatter(ScalarFormatter())
        ax1.ticklabel_format(style='scientific', axis='y', scilimits=(0,0))
    
    # Plot 2: Phlash
    ax2.plot(phlash_times, phlash_median_ne, linewidth=2.5,
             label='Phlash median', color='green')
    ax2.fill_between(phlash_times, phlash_lower, phlash_upper,
                     alpha=0.2, color='green', label='95% CI')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Time (generations ago)', fontsize=11)
    ax2.set_ylabel('N_e', fontsize=11)
    ax2.set_title('Phlash Results (Bayesian)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend(fontsize=10)
    
    plt.tight_layout()
    
    output_file = 'phlash_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")
    
    # ==================================================
    # Summary
    # ==================================================
    print("\n" + "="*70)
    print("PHLASH ANALYSIS COMPLETE")
    print("="*70)
    print("\nResults:")
    if psmc_times is not None:
        print(f"  • PSMC: {len(psmc_times)} time points")
    if phlash_success:
        print(f"  • Phlash: {len(posterior_samples)} posterior samples")
        print(f"  • Phlash Ne range: {phlash_median_ne.min():.0f} - {phlash_median_ne.max():.0f}")
    
    print("\nGenerated files:")
    print(f"  • {output_file}")
    print("\n" + "="*70)

if __name__ == '__main__':
    main()

