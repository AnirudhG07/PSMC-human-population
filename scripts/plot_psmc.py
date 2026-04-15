#!/usr/bin/env python3
"""
PSMC Plotter (fixed + enhanced)

Features:
- Proper Ne scaling (λ → Ne)
- Supports multiple input files (overlaid)
- Optional iteration gradient (-it flag)
- Smoothing for noisy small datasets
"""

import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import os

# ==================================================
# Parse PSMC
# ==================================================
def parse_psmc_iterations(filename):
    iterations = []
    current_t, current_l = [], []
    in_rd = False

    with open(filename) as f:
        for line in f:
            if line.startswith('RD'):
                if current_t:
                    iterations.append((np.array(current_t), np.array(current_l)))
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
            iterations.append((np.array(current_t), np.array(current_l)))

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
def plot_sample(ax, iterations, label, color, show_iterations=False, mu=1.25e-8):

    if not iterations:
        return

    if not show_iterations:
        # Only final iteration
        t, l = iterations[-1]

        mask = (t > 0) & (l > 0)
        t, l = t[mask], l[mask]

        if len(t) == 0:
            return

        Ne = l / (2 * mu)
        Ne = smooth(Ne)

        ax.plot(t, Ne, color=color, linewidth=2.5, label=label)

    else:
        # Plot all iterations with gradient
        n = len(iterations)

        for i, (t, l) in enumerate(iterations):
            mask = (t > 0) & (l > 0)
            t, l = t[mask], l[mask]

            if len(t) == 0:
                continue

            Ne = l / (2 * mu)
            Ne = smooth(Ne)

            alpha = 0.1 + 0.9 * (i + 1) / n
            lw = 0.5 + 2.0 * (i + 1) / n

            ax.plot(
                t,
                Ne,
                color=color,
                alpha=alpha,
                linewidth=lw,
                label=label if i == n - 1 else ""
            )


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
    args = parser.parse_args()

    # Override smoothing if disabled
    global smooth
    if args.no_smooth:
        smooth = lambda x: x

    # ==================================================
    # Plot
    # ==================================================
    fig, ax = plt.subplots(figsize=(10, 7))

    colors = plt.cm.tab10.colors

    for i, f in enumerate(args.psmc_files):
        iterations = parse_psmc_iterations(f)

        if not iterations:
            print(f"⚠️ Skipping {f} (no data)")
            continue

        label = Path(f).stem
        color = colors[i % len(colors)]

        plot_sample(
            ax,
            iterations,
            label,
            color,
            show_iterations=args.iterations
        )

    # ==================================================
    # Axes
    # ==================================================
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('Time (coalescent units)')
    ax.set_ylabel('Effective population size (Ne)')

    ax.set_title('PSMC Demographic Inference')

    ax.grid(True, which='both', alpha=0.3)

    ax.set_ylim(bottom=1e2)

    ax.legend()

    plt.tight_layout()
    plt.savefig(args.out, dpi=300)
    plt.close()

    print(f"✓ Plot saved to {args.out}")


if __name__ == "__main__":
    main()