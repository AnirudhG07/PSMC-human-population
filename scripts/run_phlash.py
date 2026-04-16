#!/usr/bin/env python3
"""
Run Phlash analysis: Takes PSMCFA file and compares with PSMC results.
Fixes the scaling bug where rescaled models were evaluated at unscaled time points.
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pickle

def parse_psmc(filename):
    """Parse PSMC output file to extract lambda and time points."""
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
        print("\nUsage: python run_phlash.py <psmcfa_file> [psmc_file] [mutation_rate]")
        print("\nExample (Human): python run_phlash.py data.psmcfa data.psmc 1.29e-8")
        print("Example (Simulated): python run_phlash.py data.psmcfa data.psmc 2e-7\n")
        return
    
    psmcfa_file = sys.argv[1]
    psmc_file = None
    mutation_rate = 1.29e-8 # Default human rate
    
    # Flexible argument parsing
    for arg in sys.argv[2:]:
        if os.path.exists(arg):
            psmc_file = arg
        else:
            try:
                mutation_rate = float(arg)
            except ValueError:
                pass

    if not os.path.exists(psmcfa_file):
        print(f"✗ Error: PSMCFA file not found: {psmcfa_file}")
        return

    print("\n" + "="*70)
    print("PHLASH ANALYSIS")
    print(f"Input: {psmcfa_file}")
    print(f"Mutation Rate: {mutation_rate}")
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
    print(f"\n[Step 2] Running Phlash inference...")
    
    phlash_success = False
    phlash_times = None
    phlash_median_ne = None
    phlash_lower = None
    phlash_upper = None
    posterior_samples = None
    
    try:
        import phlash
        
        # Define prior theta (4 * N0 * mu) assuming N0=10000
        theta_prior = 4 * 10000 * mutation_rate
        
        # Use phlash.psmc() convenience function
        print(f"  Running inference with prior theta={theta_prior:.6f}...")
        posterior_samples = phlash.psmc([psmcfa_file], theta=theta_prior)
        print(f"  ✓ Phlash complete ({len(posterior_samples)} posterior samples)")
        
        # FIX: Rescale models to generations/Ne units BEFORE evaluation
        # The 'theta' in the demographic model is per-window, so we must 
        # rescale using the per-window mutation rate.
        window_size = 100 # Default used in phlash.psmc
        rescaled_samples = [dm.rescale(mutation_rate * window_size) for dm in posterior_samples]
        
        # FIX: Extract results in rescaled units (generations)
        # Use a geometric grid over the inferred time range
        times_all = np.array([dm.eta.t[1:] for dm in rescaled_samples])
        phlash_times = np.geomspace(times_all.min(), times_all.max(), 500)
        
        # Evaluate rescaled models at rescaled time points
        Nes = np.array([dm.eta(phlash_times, Ne=True) for dm in rescaled_samples])
        
        phlash_median_ne = np.median(Nes, axis=0)
        phlash_lower = np.percentile(Nes, 2.5, axis=0)
        phlash_upper = np.percentile(Nes, 97.5, axis=0)
        
        phlash_success = True
        
    except Exception as e:
        import traceback
        print(f"  ⚠ Phlash failed: {str(e)}")
        traceback.print_exc()
        return
    
    # ==================================================
    # Save results
    # ==================================================
    print("\n[Step 3] Saving results for downstream analysis...")
    
    # Save results for downstream analysis
    pickle_file = 'phlash_output.pkl'
    with open(pickle_file, 'wb') as f:
        pickle.dump({
            'times': phlash_times,
            'median_ne': phlash_median_ne,
            'lower_ne': phlash_lower,
            'upper_ne': phlash_upper,
            'mutation_rate': mutation_rate,
            'window_size': window_size
        }, f)
    print(f"  ✓ Saved: {pickle_file}")
    print("  (Plotting is now handled by plot_phlash.py)")
    
    print("\n" + "="*70)
    print("PHLASH ANALYSIS COMPLETE")
    print("="*70)

if __name__ == '__main__':
    main()
