#!/usr/bin/env python3
"""
Generate simulated data: Creates PSMCFA and PSMC files
"""
import msprime
import numpy as np
import subprocess
import os
    
def main():

    print("\n" + "="*70)
    print("SIMULATION: Generate PSMCFA and PSMC Files")
    print("="*70)
    
    # Paths
    PSMC_BIN = "/home/anirudhgupta/bigdatabio/psmc/psmc"
    
    # ==================================================
    # Step 1: Generate simulated data with msprime
    # ==================================================
    print("\n[Step 1] Generating simulated data...")
    
    ts = msprime.sim_ancestry(
        samples=2,
        sequence_length=1e6,
        recombination_rate=1e-8,
        population_size=10000,
        random_seed=42
    )
    ts = msprime.sim_mutations(ts, rate=1.29e-8, random_seed=42)
    print(f"  ✓ Simulated {ts.num_sites} SNPs from 1 Mb sequence")
    
    # ==================================================
    # Step 2: Create PSMCFA file
    # ==================================================
    print("\n[Step 2] Creating PSMCFA file...")
    
    psmcfa_file = "simulated_validation.psmcfa"
    with open(psmcfa_file, 'w') as f:
        f.write(">simulated_1Mb\n")
        # Simple encoding: N for no het, T for het
        seq = "N" * ts.num_sites + "T" * (10000 - ts.num_sites)
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
    
    print(f"  ✓ Created {psmcfa_file} ({os.path.getsize(psmcfa_file)/1024:.1f} KB)")
    
    # ==================================================
    # Step 3: Run PSMC
    # ==================================================
    print("\n[Step 3] Running PSMC...")
    
    psmc_file = "simulated_validation.psmc"
    cmd = f"{PSMC_BIN} -N25 -t15 -r5 -p '4+25*2+4+6' -o {psmc_file} {psmcfa_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if os.path.exists(psmc_file) and os.path.getsize(psmc_file) > 1000:
        lines = sum(1 for _ in open(psmc_file))
        print(f"  ✓ PSMC complete ({lines} lines)")
    else:
        print(f"  ✗ PSMC failed")
        return
    
    # ==================================================
    # Summary
    # ==================================================
    print("\n" + "="*70)
    print("SIMULATION COMPLETE")
    print("="*70)
    print("\nGenerated files:")
    print(f"  • {psmcfa_file}")
    print(f"  • {psmc_file}")
    print("\n" + "="*70)

if __name__ == '__main__':
    main()

