#!/usr/bin/env python3
import msprime
import numpy as np
import subprocess
import os
from pathlib import Path
import random
import json

ROOT_DIR = Path(__file__).parent.parent
DATASETS_DIR = ROOT_DIR / "datasets"

PSMC_BIN = ROOT_DIR / "psmc" / "psmc"

# Simulation parameters: Increased segments for better deep-time signal
NUM_SEGMENTS = 20
SEGMENT_LENGTH = int(5e6) 
RECOMB_RATE = 1e-8
MUT_RATE = 2e-7
BIN_SIZE = 100

# Demography: Bottleneck and ancient structure
def build_demography():
    demography = msprime.Demography()
    demography.add_population(name="pop", initial_size=10000)
    
    # Bottleneck
    demography.add_population_parameters_change(time=10000, initial_size=1200, population="pop")
    # Recovery
    demography.add_population_parameters_change(time=25000, initial_size=12000, population="pop")
    # Ancient expansion
    demography.add_population_parameters_change(time=100000, initial_size=8000, population="pop")

    return demography

def save_truth(demography, filename):
    """Saves the demographic steps for plotting."""
    # Since we know our model, we hardcode for plotting clarity
    history = [
        {"time": 0, "Ne": 10000},
        {"time": 10000, "Ne": 1200},
        {"time": 25000, "Ne": 12000},
        {"time": 100000, "Ne": 8000},
        {"time": 1000000, "Ne": 8000}
    ]
    with open(filename, 'w') as f:
        json.dump(history, f)

# Convert multiple tree sequences → single PSMCFA
def ts_to_psmcfa(ts_list, output_file):
    with open(output_file, "w") as f:
        for idx, ts in enumerate(ts_list):
            L = int(ts.sequence_length)
            seq = np.zeros(L, dtype=np.uint8)
            for var in ts.variants():
                if var.genotypes[0] != var.genotypes[1]:
                    seq[int(var.site.position)] = 1
            
            num_bins = L // BIN_SIZE
            binned = []
            for i in range(num_bins):
                binned.append("K" if np.any(seq[i*BIN_SIZE:(i+1)*BIN_SIZE]) else "T")
            
            f.write(f">scaffold_{idx}\n")
            for i in range(0, len(binned), 60):
                f.write("".join(binned[i:i+60]) + "\n")

def main():
    DATASETS_DIR.mkdir(exist_ok=True)
    seed = random.randint(1, 10_000_000)
    tag = f"{seed}"
    
    psmcfa_file = DATASETS_DIR / f"sim_data_{tag}.psmcfa"
    psmc_file = DATASETS_DIR / f"sim_data_{tag}.psmc"
    truth_file = DATASETS_DIR / f"sim_data_{tag}_truth.json"

    print(f"\n[Info] Seed: {seed} | Tag: {tag}")
    
    demography = build_demography()
    save_truth(demography, truth_file)

    ts_list = []
    print(f"[Step 1] Simulating {NUM_SEGMENTS} segments of {SEGMENT_LENGTH/1e6}Mb...")
    for i in range(NUM_SEGMENTS):
        ts = msprime.sim_ancestry(
            samples=1, ploidy=2, sequence_length=SEGMENT_LENGTH,
            recombination_rate=RECOMB_RATE, demography=demography, random_seed=seed + i
        )
        ts = msprime.sim_mutations(ts, rate=MUT_RATE, random_seed=seed + i)
        ts_list.append(ts)

    print("[Step 2] Creating multi-scaffold PSMCFA...")
    ts_to_psmcfa(ts_list, psmcfa_file)

    print("[Step 3] Running PSMC...")
    # -t30 for deep reach, -p "2+28*2+6" for high recent resolution
    subprocess.run([str(PSMC_BIN), "-N25", "-t30", "-r5", "-p", "2+28*2+6", "-o", str(psmc_file), str(psmcfa_file)])
    print(f"  ✓ Done. Files: {psmc_file.name}, {truth_file.name}")

if __name__ == "__main__":
    main()
