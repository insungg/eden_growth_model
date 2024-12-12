"""Utility functions for data processing and visualization for analysis."""
import numpy as np
import itertools
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib import gridspec
from matplotlib.colors import Normalize

from source.draw_lattice import draw_lattice_heatmap
from source.decoder import get_landscape

def calculate_survival_probabilities(landscape, nu_list, s_list, species=2):
    survival_probabilities = {}

    for nu, s in itertools.product(nu_list, s_list):
        phase_params = (nu, s)  # Create a unique key for each combination
        
        if phase_params in landscape:
            snapshots = landscape[phase_params]["snapshots"]
            
            if snapshots:  # Ensure there is data
                # Convert snapshots to a numpy array
                # Shape: (number of snapshots, rows, cols)
                snapshots_array = np.array(snapshots)
                
                # Calculate survival probability:
                # Count of mutants (value == 2) divided by the number of snapshots
                survival_probability = np.sum(snapshots_array == species, axis=0) / len(snapshots)
                
                # Store the result
                survival_probabilities[phase_params] = survival_probability
            else:
                print(f"No snapshots available for phase_params (nu={nu}, s={s}).")
                survival_probabilities[phase_params] = None
        else:
            print(f"No data found for phase_params (nu={nu}, s={s}).")
            survival_probabilities[phase_params] = None
    
    return survival_probabilities

def calculate_mutant_frequency(snapshot):
    return np.sum(snapshot == 2) / snapshot.size

def calculate_mean_mutant_frequency(data_path, prefix, verbose=True):
    # Initialize dictionary to store aggregated results for (nu, s)
    results = {}
    # Iterate over all landscapes with given prefix
    for filename in os.listdir(data_path):
        if filename.startswith(prefix) and filename.endswith(".bin"):
            file_path = os.path.join(data_path, filename)
            if verbose: print(f"Processing file: {filename}")

            landscape_results = get_landscape(data_path, prefix=filename, verbose=False)
            for phase_params, data in landscape_results.items():
                snapshots = data["snapshots"]
                # Initialize storage if this is the first file for this (nu, s)
                if phase_params not in results:
                    results[phase_params] = []
                # Compute mutant frequencies for all snapshots in the current file
                for snapshot in snapshots:
                    results[phase_params].append(calculate_mutant_frequency(snapshot))
                del snapshots
            del landscape_results

    # Compute the final average mutant frequency for each (nu, s)
    final_results = {}
    for phase_params, frequencies in results.items():
        final_results[phase_params] = np.mean(frequencies)

    # Print the results
    if verbose:
        print("Final Aggregated Results (sorted by phase_params):")
        sorted_phase_params = sorted(final_results.keys(), key=lambda x: (x[0], x[1]))
        for phase_params in sorted_phase_params:
            nu, s = phase_params
            avg_frequency = final_results[phase_params]
            print(f"Phase Params: (nu={nu}, s={s}) -> Avg Mutant Frequency: {avg_frequency:.4f}")

    return final_results

def create_probability_pairplot(probabilities, nus, ss, hotspots=None, title=None, fig_name = None, color_map="turbo", vmin=0, vmax=0.5):
    n_nu = len(nus)
    n_s = len(ss)

    # spacing fine-tuning
    fig = plt.figure(figsize=(5 * n_nu + 2, 5 * n_s))  # Add extra width for colorbar
    gs = gridspec.GridSpec(n_s, n_nu, figure=fig, wspace=0.01, hspace=0.01)
    # plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.005, hspace=0.005)

    for i, s in enumerate(ss[::-1]): 
        for j, nu in enumerate(nus):
            ax = fig.add_subplot(gs[i, j])
            density = probabilities[(nu, s)]
            draw_lattice_heatmap(density, hotspots=hotspots,ax=ax, color_map=color_map, vmin=vmin, vmax=vmax)
            ax.set_xticks([])
            ax.set_yticks([])

    # Add global axes
    global_ax = fig.add_subplot(111, frame_on=False)  # Add a new axis on top of all subplots
    global_ax.set_xticks([i + 0.5 for i in range(n_nu)])  # Centered x-ticks
    global_ax.set_yticks([i + 0.5 for i in range(n_s)])   # Centered y-ticks
    global_ax.set_xticklabels([f"{nu}" for nu in nus], fontsize=26)  # Larger x-axis labels
    global_ax.set_yticklabels([f"{s:.2f}" for s in ss], fontsize=26)  # Larger y-axis labels
    global_ax.tick_params(labelsize=26)
    global_ax.set_xlabel(r"Intensity $\nu$", fontsize=35, labelpad=25)
    global_ax.set_ylabel(r"Selection $s$", fontsize=35, labelpad=25)
    global_ax.xaxis.set_label_position("bottom")
    global_ax.yaxis.set_label_position("left")
    global_ax.set_xlim(0, n_nu)
    global_ax.set_ylim(0, n_s)
    global_ax.xaxis.tick_bottom()
    global_ax.yaxis.tick_left()
    global_ax.axhline(y=0, color="black", linestyle="-", linewidth=2)
    global_ax.axvline(x=0, color="black", linestyle="-", linewidth=2)

    # Add shared colorbar
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(color_map)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Required for ScalarMappable
    cbar = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # Adjusted position of the colorbar
    cbar = fig.colorbar(sm, cax=cbar)
    cbar.set_label(r"$M(x,y)$", fontsize=35)
    cbar.ax.tick_params(labelsize=26)

    if title:
        fig.suptitle(title, fontsize=35, y=0.95)
    if fig_name:
        fig.savefig(fig_name, dpi=300, bbox_inches="tight") 

    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.show()

def create_frequency_heatmap(frequency, nu_list=None, s_list=None, title=None, fig_name=None, ax=None):
    # Use provided nu_list and s_list or determine them from the frequency dictionary
    if nu_list is None:
        nu_list = sorted(set(key[0] for key in frequency.keys()))
    if s_list is None:
        s_list = sorted(set(key[1] for key in frequency.keys()))
    # Validate the filtered lists
    if not nu_list or not s_list:
        raise ValueError("Filtered nu_list or s_list results in an empty plot. Please check the provided lists.")

    # Create the grid based on nu_list and s_list
    fm_grid = np.zeros((len(s_list), len(nu_list)))  # Rows: s, Columns: ν
    for i, s in enumerate(s_list):
        for j, nu in enumerate(nu_list):
            # Use frequency[(nu, s)] if it exists; otherwise, default to 0
            fm_grid[i, j] = frequency.get((nu, s), 0)

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 12))

    im = ax.imshow(fm_grid, extent=[min(nu_list), max(nu_list), min(s_list), max(s_list)],
                   origin='lower', aspect='auto', cmap='coolwarm')
    ax.set_xlabel(r'Intensity, $\nu$', fontsize=20)
    ax.set_ylabel(r'Selection, $s$', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=17)

    # Add colorbar to the current figure
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r'$f_m$', fontsize=20, rotation=0, labelpad=20)
    cbar.ax.tick_params(labelsize=17)
    
    if title: ax.set_title(title, fontsize=30)
    if fig is not None:
        plt.tight_layout()
        plt.show()
    if fig_name: plt.savefig(fig_name, dpi=300)