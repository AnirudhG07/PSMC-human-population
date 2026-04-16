#!/usr/bin/env python3
import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
import json
import argparse

def plot_truth(ax, truth_file, max_time):
    if not os.path.exists(truth_file):
        return
    with open(truth_file, 'r') as f:
        truth_data = json.load(f)
    truth_data.sort(key=lambda x: x['time'])
    times = [x['time'] for x in truth_data]
    nes = [x['Ne'] for x in truth_data]
    
    # Extend to max_time
    plot_times = times + [max_time]
    plot_nes = nes + [nes[-1]]
    
    ax.step(plot_times, plot_nes, where='post', 
             color='black', linestyle='--', linewidth=2.0, label='Ground Truth', zorder=5, alpha=0.7)

def main():
    parser = argparse.ArgumentParser(description="Plot Phlash results from pickle output")
    parser.add_argument("--pkl", default="phlash_output.pkl", help="Phlash pickle file")
    parser.add_argument("--truth", help="Ground truth JSON file")
    parser.add_argument("--out", help="Output filename prefix")
    parser.add_argument("--zoomed", action="store_true", help="Generate only a zoomed plot")
    parser.add_argument("--full", action="store_true", help="Generate only a full range plot")
    parser.add_argument("--fix-scaling", action="store_true", help="Apply 100x scaling fix for old pkl files")
    
    args = parser.parse_args()

    if not os.path.exists(args.pkl):
        print(f"Error: {args.pkl} not found.")
        return

    print(f"Loading data from {args.pkl}...")
    with open(args.pkl, 'rb') as f:
        data = pickle.load(f)

    times = data['times']
    median_ne = data['median_ne']
    lower = data['lower_ne']
    upper = data['upper_ne']

    if args.fix_scaling:
        print("Applying 100x scaling fix...")
        times /= 100.0
        median_ne /= 100.0
        lower /= 100.0
        upper /= 100.0

    prefix = args.out if args.out else "phlash_analysis"
    
    # helper for single plot
    def create_plot(filename, y_lim=None, title_suffix=""):
        plt.figure(figsize=(10, 6))
        
        if args.truth:
            plot_truth(plt.gca(), args.truth, times.max())

        plt.step(times, median_ne, linewidth=2.5,
                 label='Phlash median', color='red', where='post')
        plt.fill_between(times, lower, upper,
                         alpha=0.2, color='red', label='95% CI', step='post')
        
        plt.xscale('log')
        plt.yscale('log')
        
        if y_lim:
            plt.ylim(y_lim)
        
        plt.xlabel('Time (generations ago)', fontsize=11)
        plt.ylabel('Effective Population Size ($N_e$)', fontsize=11)
        plt.title(f'Phlash Demographic History {title_suffix}', fontsize=12, fontweight='bold')
        plt.grid(True, alpha=0.3, which='both')
        plt.legend(fontsize=10)
        
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"✓ Saved plot to: {filename}")
        plt.close()

    # Determine what to plot
    if not args.zoomed and not args.full:
        # Plot both by default
        create_plot(f"{prefix}_full.png", title_suffix="(Full Range)")
        
        y_min = max(100, median_ne.min() / 5)
        y_max = median_ne.max() * 5
        create_plot(f"{prefix}_zoomed.png", y_lim=(y_min, y_max), title_suffix="(Zoomed)")
    
    elif args.full:
        create_plot(f"{prefix}_full.png", title_suffix="(Full Range)")
        
    elif args.zoomed:
        y_min = max(100, median_ne.min() / 5)
        y_max = median_ne.max() * 5
        create_plot(f"{prefix}_zoomed.png", y_lim=(y_min, y_max), title_suffix="(Zoomed)")

if __name__ == '__main__':
    main()
