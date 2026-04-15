#!/usr/bin/env python3
"""
For Human Genome Data

Optimized PSMC pipeline for demographic analysis.
Supports running specific samples via command line arguments.
Example: uv run scripts/pipeline_psmc.py NA18561.bam NA12878.bam
"""

import os
import sys
import subprocess
import multiprocessing as mp
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import time

# Base Paths
ROOT_DIR = "/home/anirudhgupta/project/BigDataBio"
PSMC_BIN = os.path.join(ROOT_DIR, "psmc/psmc")
PSMC_UTILS = os.path.join(ROOT_DIR, "psmc/utils")
DATASETS_DIR = os.path.join(ROOT_DIR, "datasets")
REF_GENOME = os.path.join(DATASETS_DIR, "hs37d5.fa.gz")

def log(msg, level="INFO"):
    """Pretty logging"""
    prefix = f"[{level}]" if level != "SUCCESS" else "[\u2713]"
    print(f"{prefix} {msg}", flush=True)

def file_valid(path, min_size=100):
    """Check if file exists and has minimum size"""
    if not os.path.exists(path):
        return False
    size = os.path.getsize(path)
    return size >= min_size

def run_cmd(cmd, label=""):
    """Run command and handle errors"""
    try:
        log(f"DEBUG: Executing command for {label}: {cmd}")
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=14400  # 4hr timeout per step
        )
        if result.returncode != 0:
            log(f"FAILED: {label}", "ERROR")
            log(f"  exit code: {result.returncode}", "ERROR")
            log(f"  stderr: {result.stderr[:1000]}", "ERROR")
            return False
        return True
    except subprocess.TimeoutExpired:
        log(f"TIMEOUT: {label}", "ERROR")
        return False
    except Exception as e:
        log(f"ERROR: {label} - {str(e)}", "ERROR")
        return False

def ensure_bam_index(bam_file):
    """Ensure BAM index exists, generate if missing"""
    possible_indices = [f"{bam_file}.bai", bam_file.replace(".bam", ".bai")]
    for idx in possible_indices:
        if os.path.exists(idx):
            return True
    
    log(f"Index missing for {bam_file}, generating...")
    if not run_cmd(f"samtools index -@ 4 {bam_file}", f"Indexing {bam_file}"):
        return False
    return True

def get_avg_depth(bam_file):
    """Estimate average depth by sampling multiple 100kb regions across chromosomes 1, 2, and 3"""
    regions = ["1:10000000-10100000", "2:10000000-10100000", "3:10000000-10100000"]
    depths = []
    for reg in regions:
        try:
            # -a includes zero-coverage positions. Window size is 100,000 bp.
            cmd = f"samtools depth -a -r {reg} {bam_file} | awk '{{sum+=$3}} END {{print sum/100000}}'"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            val = float(result.stdout.strip()) if result.stdout.strip() else 0
            depths.append(val)
        except:
            continue
    
    return sum(depths) / len(depths) if depths else 0

def generate_fq(sample_name, bam_file):
    """BAM -> FQ.GZ following README exactly"""
    fq_file = f"{sample_name}.fq.gz"
    
    if os.path.exists(fq_file) and os.path.getsize(fq_file) < 1000:
        os.remove(fq_file)

    if file_valid(fq_file, min_size=10*1024*1024):
        log(f"FASTQ exists for {sample_name}, skipping", "SKIP")
        return fq_file
    
    log(f"Generating FASTQ for {sample_name}...")
    
    # Estimate depth to set -d threshold (1/3 of average depth is recommended for low coverage)
    avg_depth = get_avg_depth(bam_file)
    min_depth = max(2, int(avg_depth / 3))
    # PSMC README recommends -d = 1/3 and -D = 2x average depth
    max_depth = max(50, int(avg_depth * 2.0)) 
    
    log(f"Detected true avg depth: {avg_depth:.2f}x. Using filters: -d {min_depth} -D {max_depth}")
    
    cmd = (f"bcftools mpileup --threads 4 -C50 -f {REF_GENOME} {bam_file} | "
           f"bcftools call --threads 4 -c -V indels | "
           f"vcfutils.pl vcf2fq -d {min_depth} -D {max_depth} | gzip > {fq_file}")
    
    if not run_cmd(cmd, f"FASTQ generation: {sample_name}"):
        return None
    return fq_file

def generate_psmcfa(sample_name):
    """FQ.GZ -> PSMCFA following README exactly"""
    fq_file = f"{sample_name}.fq.gz"
    psmcfa_file = f"{sample_name}.psmcfa"
    
    # Check if PSMCFA already valid (1MB min for human genome)
    if file_valid(psmcfa_file, min_size=1*1024*1024):
        log(f"PSMCFA exists for {sample_name}, skipping", "SKIP")
        return psmcfa_file
    
    if os.path.exists(psmcfa_file):
        log(f"DEBUG: Found corrupted/empty {psmcfa_file}, deleting...")
        os.remove(psmcfa_file)
    
    log(f"Generating PSMCFA for {sample_name}...")
    cmd = f"{PSMC_UTILS}/fq2psmcfa -q20 {fq_file} > {psmcfa_file}"
    if not run_cmd(cmd, f"PSMCFA generation: {sample_name}"):
        return None
    return psmcfa_file

def generate_psmc(sample_name):
    """PSMCFA -> PSMC following README exactly"""
    psmcfa_file = f"{sample_name}.psmcfa"
    psmc_file = f"{sample_name}.psmc"
    
    if file_valid(psmc_file, min_size=10000):
        log(f"PSMC exists for {sample_name}, skipping", "SKIP")
        return psmc_file
    
    log(f"Running PSMC for {sample_name}...")
    cmd = f"{PSMC_BIN} -N25 -t15 -r5 -p '4+25*2+4+6' -o {psmc_file} {psmcfa_file}"
    if not run_cmd(cmd, f"PSMC: {sample_name}"):
        return None
    return psmc_file

def process_sample(sample_name, bam_file):
    """Full pipeline for a single sample"""
    log(f"--- Processing {sample_name} ---")
    os.chdir(DATASETS_DIR)
    
    if not ensure_bam_index(bam_file):
        return False
    if not generate_fq(sample_name, bam_file):
        return False
    if not generate_psmcfa(sample_name):
        return False
    if not generate_psmc(sample_name):
        return False
    
    log(f"✓ Completed {sample_name}", "SUCCESS")
    return True

def main():
    args = sys.argv[1:]
    if not args:
        # Default to all BAM files in datasets directory if no args provided
        args = [str(p.name) for p in Path(DATASETS_DIR).glob("*.bam")]
    
    target_samples = {}
    for arg in args:
        bam_path = Path(arg)
        # Try finding the file directly or in DATASETS_DIR
        if not bam_path.exists():
            bam_path = Path(DATASETS_DIR) / arg
        
        if not bam_path.exists():
            log(f"BAM file not found: {arg}", "ERROR")
            continue
            
        bam_path = bam_path.absolute()
        # Derive sample name from filename (e.g., NA12878.bam -> NA12878)
        sample_name = bam_path.name.split('.')[0]
        target_samples[sample_name] = str(bam_path)

    if not target_samples:
        log("No valid samples to process. Example: python scripts/pipeline.py NA18561.bam", "ERROR")
        return

    log(f"Starting pipeline for samples: {', '.join(target_samples.keys())}")
    
    # Ensure max_workers is at least 1 and capped at 3
    num_workers = max(1, min(len(target_samples), 3))
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = {s: executor.submit(process_sample, s, b) for s, b in target_samples.items()}
        for s, f in futures.items():
            try:
                if not f.result(): log(f"Sample {s} failed", "ERROR")
            except Exception as e:
                log(f"Error processing {s}: {str(e)}", "ERROR")
    
    log("=" * 60)
    log("Pipeline finished. Use 'python scripts/plot_psmc.py <psmc_files>' to generate plots.")
    log("=" * 60)

if __name__ == '__main__':
    main()
