#!/usr/bin/env python3

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# ---------------------------
# Read command-line arguments
# ---------------------------
epsilon = sys.argv[1]
rmin = sys.argv[2]
r1 = sys.argv[3]
rc = sys.argv[4]
bFC = sys.argv[5]
aFC = sys.argv[6]
Friction = sys.argv[7]
NFE = sys.argv[8]

skiprow = 26

# -----------------------------
# Function to plot radial density
# -----------------------------
def radial_density(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE, Fitness):
    base_filename = f"{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}"
    rdf_files = {
        "T-T": f"g_T_T_{base_filename}.xvg",
        "H-H": f"g_H_H_{base_filename}.xvg",
        "H-T": f"g_T_H_{base_filename}.xvg"
    }

    # Check if at least one RDF file exists
    if not os.path.isfile(rdf_files["T-T"]):
        print(f"Simulation with parameters {base_filename} encountered an error. RDF cannot be plotted.")
        return None

    fig, ax = plt.subplots(figsize=(6, 4))

    for label, filename in rdf_files.items():
        if os.path.isfile(filename):
            data = np.loadtxt(filename, skiprows=skiprow)
            ax.plot(data[:, 0], data[:, 1], '-o', label=label)
        else:
            print(f"Warning: {filename} not found. Skipping.")

    ax.set_xlabel('Distance (Å)', fontsize=12)
    ax.set_ylabel('g(r)', fontsize=12)
    ax.legend(loc='upper right')
    plt.tight_layout()

    output_file = f"RDF_{base_filename}_{Fitness}.png"
    plt.savefig(output_file)
    print(f"RDF plot saved as {output_file}")

# -----------------------------
# Run the function
# -----------------------------
radial_density(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE, Fitness)
