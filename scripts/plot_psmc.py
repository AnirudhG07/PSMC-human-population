#!/usr/bin/env python3
"""
PSMC Plotter (fixed + enhanced)

Features:
- Proper Ne scaling (λ → Ne)
- Supports multiple input files (overlaid)
- Optional iteration gradient (-it flag)
- Smoothing for noisy small datasets
- Ground Truth comparison (clipped to PSMC range)
"""

import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import os
import json

# ==================================================
# Truth Plotter
# ==================================================
def plot_truth(ax, truth_file, max_time):
    if not os.path.exists(truth_file):
        print(f"⚠️ Truth file {truth_file} not found")
        return

    with open(truth_file) as f:
        history = json.load(f)
    
    times = [h["time"] for h in history]
    nes = [h["Ne"] for h in history]
    
    # Filter to match PSMC reach
    filtered_times = []
    filtered_nes = []
    for t, n in zip(times, nes):
        if t <= max_time:
            filtered_times.append(t)
            filtered_nes.append(n)
        else:
            # Add one last point at the boundary
            filtered_times.append(max_time)
            filtered_nes.append(n)
            break
            
    ax.step(filtered_times, filtered_nes, color="black", linestyle="--", linewidth=1.5, 
            label="Ground Truth", where='post', alpha=0.7)


# ==================================================
# Parse PSMC
# ==================================================
def parse_psmc_iterations(filename):
    iterations = []
    current_t, current_l = [], []
    in_rd = False
    theta_0 = None

    with open(filename) as f:
        for line in f:
            if line.startswith('TR'):
                theta_0 = float(line.split()[1])
            elif line.startswith('RD'):
                if current_t:
                    iterations.append((np.array(current_t), np.array(current_l), theta_0))
                current_t, current_l = [], []
                in_rd = True

            elif line.startswith('RS') and in_rd:
                p = line.split()
                t = float(p[2])
                l = float(p[3])

                if t > 0 and l > 0:
                    current_t.append(t)
                    current_l.append(l)

        if current_t:
            iterations.append((np.array(current_t), np.array(current_l), theta_0))

    return iterations


# ==================================================
# Smoothing (important for small Mb data)
# ==================================================
def smooth(y, window=5):
    if len(y) < window:
        return y
    return np.convolve(y, np.ones(window)/window, mode='same')


# ==================================================
# Plotting
# ==================================================
def plot_sample(ax, iterations, label, color, show_iterations=False, no_smooth=False, mu=2e-7, s=100):

    if not iterations:
        return None

    def process_it(it):
        t, l, theta_0 = it
        mask = (t > 0) & (l > 0)
        t, l = t[mask], l[mask]
        if len(t) == 0:
            return None, None
        
        N_0 = theta_0 / (4 * mu * s)
        Ne = l * N_0
        T = 2 * N_0 * t
        
        if not no_smooth:
            Ne = smooth(Ne)
        return T, Ne

    max_t_inferred = 0

    if not show_iterations:
        # Only final iteration
        T, Ne = process_it(iterations[-1])
        if T is not None:
            ax.step(T, Ne, color=color, linewidth=2.5, label=f"{label} (PSMC)", where='post')
            max_t_inferred = T[-1]
    else:
        # Plot all iterations with gradient
        n = len(iterations)
        for i, it in enumerate(iterations):
            T, Ne = process_it(it)
            if T is None:
                continue

            alpha = min(1.0, 0.1 + 0.9 * (i + 1) / n)
            lw = 0.5 + 2.0 * (i + 1) / n

            ax.step(
                T,
                Ne,
                color=color,
                alpha=alpha,
                linewidth=lw,
                label=f"{label} (Final)" if i == n - 1 else "",
                where='post'
            )
            if i == n - 1:
                max_t_inferred = T[-1]
    
    return max_t_inferred


# ==================================================
# Main
# ==================================================
def main():
    parser = argparse.ArgumentParser(description="PSMC Plotter (fixed scaling)")
    parser.add_argument("psmc_files", nargs="+", help="PSMC output files")
    parser.add_argument("-o", "--out", default="psmc_plot.png", help="Output file")
    parser.add_argument("-it", "--iterations", action="store_true",
                        help="Show all iterations with gradient")
    parser.add_argument("--no-smooth", action="store_true",
                        help="Disable smoothing")
    parser.add_argument("--mu", type=float, default=2e-7, help="Mutation rate")
    parser.add_argument("--truth", help="JSON file containing true demography")
    args = parser.parse_args()

    # ==================================================
    # Plot
    # ==================================================
    fig, ax = plt.subplots(figsize=(10, 7))

    colors = plt.cm.tab10.colors
    overall_max_t = 0

    # First pass: Plot samples and find max inferred time
    for i, f in enumerate(args.psmc_files):
        iterations = parse_psmc_iterations(f)

        if not iterations:
            print(f"⚠️ Skipping {f} (no data)")
            continue

        label = Path(f).stem
        color = colors[i % len(colors)]

        m_t = plot_sample(
            ax,
            iterations,
            label,
            color,
            show_iterations=args.iterations,
            no_smooth=args.no_smooth,
            mu=args.mu
        )
        if m_t:
            overall_max_t = max(overall_max_t, m_t)

    # Second pass: Plot Truth clipped to max inferred time
    if args.truth and overall_max_t > 0:
        plot_truth(ax, args.truth, overall_max_t)

    # ==================================================
    # Axes
    # ==================================================
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('Time (generations)')
    ax.set_ylabel('Effective population size (Ne)')

    ax.set_title('PSMC Demographic Inference vs Ground Truth')

    ax.grid(True, which='both', alpha=0.3)

    # Auto-adjust limits
    ax.set_ylim(bottom=1e2)
    # Set x-limits based on data
    ax.set_xlim(left=1e3, right=overall_max_t * 1.1 if overall_max_t > 0 else 1e6)

    ax.legend(frameon=True, facecolor='white', framealpha=0.9)

    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    plt.close()

    print(f"✓ Plot saved to {args.out}")


if __name__ == "__main__":
    main()
