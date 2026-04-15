import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from pathlib import Path
import argparse

def parse_psmc_iterations(filename):
    """Parse all iterations from PSMC output file"""
    iterations = []
    current_t, current_l = [], []
    in_rd = False
    try:
        with open(filename) as f:
            for line in f:
                if line.startswith('RD'):
                    if current_t:
                        iterations.append((np.array(current_t), np.array(current_l)))
                    current_t, current_l = [], []
                    in_rd = True
                elif line.startswith('RS') and in_rd:
                    p = line.split()
                    current_t.append(float(p[2]))
                    current_l.append(float(p[3]))
            if current_t:
                iterations.append((np.array(current_t), np.array(current_l)))
    except Exception as e:
        print(f"Error parsing {filename}: {e}")
    return iterations

def plot_sample(ax, iterations, label, base_color):
    """Plot iterations with a gradient from light to dark"""
    n = len(iterations)
    if n == 0: return
    
    for i, (t, l) in enumerate(iterations):
        # Final iteration is darkest and thickest
        is_final = (i == n - 1)
        alpha = 0.1 + 0.8 * (i + 1) / n
        linewidth = 0.5 + 1.5 * (i + 1) / n
        
        ax.plot(t, l, color=base_color, alpha=alpha, linewidth=linewidth, 
                label=label if is_final else "")

def main():
    parser = argparse.ArgumentParser(description="Plot PSMC results with iteration gradients")
    parser.add_argument("psmc_files", nargs="+", help="PSMC output files to plot")
    parser.add_argument("--out-dir", default=".", help="Directory to save the plot")
    parser.add_argument("--name", help="Custom name for the output file")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    
    # Group files into Human and Simulated
    human_files = [f for f in args.psmc_files if Path(f).name.startswith('NA')]
    simulated_files = [f for f in args.psmc_files if 'simulated' in Path(f).name.lower()]
    other_files = [f for f in args.psmc_files if f not in human_files and f not in simulated_files]

    groups = []
    if len(args.psmc_files) == 1:
        # If only one file, just plot it
        groups = [("single", args.psmc_files)]
    else:
        # Split into groups as requested
        if human_files: groups.append(("human", human_files))
        if simulated_files: groups.append(("simulated", simulated_files))
        if other_files: groups.append(("other", other_files))

    # Standard colors for plotting
    colors = plt.cm.tab10.colors

    for group_name, files in groups:
        fig, ax = plt.subplots(figsize=(11, 7))
        
        for i, f_path in enumerate(files):
            iterations = parse_psmc_iterations(f_path)
            if iterations:
                label = Path(f_path).stem
                color = colors[i % len(colors)]
                plot_sample(ax, iterations, label, color)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Scaled time (2N₀ generations)')
        ax.set_ylabel('Scaled population size (λ)')
        ax.set_title(f'PSMC: {group_name.capitalize()} Demographic Inference')
        ax.legend()
        ax.grid(True, alpha=0.3, which='both')
        
        if args.name and len(groups) == 1:
            output_name = args.name if args.name.endswith('.png') else f"{args.name}.png"
        else:
            if group_name == "single":
                output_name = f"{Path(files[0]).stem}_plot.png"
            else:
                output_name = f"psmc_{group_name}_combined.png"
        
        output_path = os.path.join(args.out_dir, output_name)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
        plt.close()
        print(f"✓ {group_name.capitalize()} plot saved to {output_path}")

if __name__ == '__main__':
    main()
