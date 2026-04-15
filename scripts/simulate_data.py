#!/usr/bin/env python3
__about__ = """
PSMC Simulated data
- Variable demography
- Proper SNP → sequence mapping
- 100 bp binning (PSMC style)
- Unique output per run (seed-based tag)
"""

import msprime
import numpy as np
import subprocess
import os
from pathlib import Path
import random

ROOT_DIR = Path(__file__).parent.parent
DATASETS_DIR = ROOT_DIR / "datasets"

PSMC_BIN = ROOT_DIR / "psmc" / "psmc"

SEQ_LENGTH = int(5e6) 
RECOMB_RATE = 1e-8
MUT_RATE = 2e-7
BIN_SIZE = 10


# Demography
def build_demography():
    demography = msprime.Demography()
    demography.add_population(name="pop", initial_size=10000)

    # Strong oscillations (amplified for small data)
    demography.add_population_parameters_change(time=1000, initial_size=1000, population="pop")
    demography.add_population_parameters_change(time=3000, initial_size=50000, population="pop")
    demography.add_population_parameters_change(time=8000, initial_size=2000, population="pop")
    demography.add_population_parameters_change(time=20000, initial_size=80000, population="pop")
    demography.add_population_parameters_change(time=40000, initial_size=3000, population="pop")
    demography.add_population_parameters_change(time=100000, initial_size=60000, population="pop")

    return demography

# Convert tree sequence → PSMCFA
def ts_to_psmcfa(ts, output_file):
    L = int(ts.sequence_length)

    print(f"  → Building sequence of length {L:,}")

    # initialize
    seq = np.zeros(L, dtype=np.uint8)  # 0 = homozygous, 1 = heterozygous

    # mark heterozygous sites
    for var in ts.variants():
        g = var.genotypes
        if g[0] != g[1]:
            pos = int(var.site.position)
            if pos < L:
                seq[pos] = 1

    # bin into 100bp windows
    print(f"  → Binning into {BIN_SIZE} bp windows...")
    num_bins = L // BIN_SIZE
    binned = []

    for i in range(num_bins):
        chunk = seq[i * BIN_SIZE:(i + 1) * BIN_SIZE]
        if np.any(chunk):
            binned.append("T")
        else:
            binned.append("N")

    print(f"  → Writing PSMCFA ({len(binned):,} bins)...")

    with open(output_file, "w") as f:
        f.write(">simulated\n")
        for i in range(0, len(binned), 60):
            f.write("".join(binned[i:i+60]) + "\n")


# ==================================================
# Main
# ==================================================
def main():

    print("\n" + "="*70)
    print("SIMULATION: Realistic PSMC Data")
    print("="*70)

    # Ensure output dir
    DATASETS_DIR.mkdir(exist_ok=True)

    # Random seed + tag
    seed = random.randint(1, 10_000_000)
    tag = f"{seed}"

    print(f"\n[Info] Seed: {seed}")
    print(f"[Info] Tag: {tag}")

    # Output files
    psmcfa_file = DATASETS_DIR / f"simulated_data_{tag}.psmcfa"
    psmc_file = DATASETS_DIR / f"simulated_data_{tag}.psmc"

    # Simulate ancestry
    print("\n[Step 1] Simulating ancestry...")

    demography = build_demography()

    ts = msprime.sim_ancestry(
        samples=1,
        ploidy=2,
        sequence_length=SEQ_LENGTH,
        recombination_rate=RECOMB_RATE,
        demography=demography,
        random_seed=seed
    )

    ts = msprime.sim_mutations(ts, rate=MUT_RATE, random_seed=seed)

    print(f"  ✓ Sites: {ts.num_sites:,}")

    # Convert to PSMCFA
    print("\n[Step 2] Creating PSMCFA...")

    ts_to_psmcfa(ts, psmcfa_file)

    print(f"  ✓ {psmcfa_file}")

    # Run PSMC
    print("\n[Step 3] Running PSMC...")

    cmd = [
        str(PSMC_BIN),
        "-N25",
        "-t15",
        "-r5",
        "-p", "4+10*2+4+6",
        "-o", str(psmc_file),
        str(psmcfa_file)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if os.path.exists(psmc_file) and os.path.getsize(psmc_file) > 1000:
        lines = sum(1 for _ in open(psmc_file))
        print(f"  ✓ PSMC complete ({lines} lines)")
    else:
        print("  ✗ PSMC failed")
        print(result.stderr)
        return

    print("\n" + "="*70)
    print("DONE")
    print("="*70)

    print("\nGenerated:")
    print(f"- {psmcfa_file}")
    print(f"- {psmc_file}")

if __name__ == "__main__":
    main()