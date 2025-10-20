#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
import statistics
import Conversion

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

Sample_Num = 25
Sample_Pass = 1680000

# -----------------------------
# Utility functions
# -----------------------------
def average(lst):
    return np.mean(lst)

def read_column(filename, col_index=1, step=1000):
    """Read the specified column from a GROMACS .xvg file with optional step sampling."""
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"{filename} not found.")
    
    data = []
    with open(filename) as f:
        for i, line in enumerate(f):
            if i >= Sample_Num and i < Sample_Pass:
                values = line.split()
                if len(values) > col_index:
                    data.append(float(values[col_index]))
    return data[::step]

# -----------------------------
# Core functions
# -----------------------------
def Time(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE):
    filename = f"Pxx_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.xvg"
    time_data = read_column(filename, col_index=0)
    return [t for t in time_data]

def Pressure_Components(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE):
    Pxx = read_column(f"Pxx_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.xvg")
    Pyy = read_column(f"Pyy_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.xvg")
    Pzz = read_column(f"Pzz_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.xvg")
    return Pxx, Pyy, Pzz

def Box_Area(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE):
    LX = 20
    LZ = 15
    return LX, LZ  # Non-periodic dimensions (surrounded by water)

def Line_Tension(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE):
    Pxx, Pyy, Pzz = Pressure_Components(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE)
    LX, LZ = Box_Area(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE)

    LT = [
        Conversion.LineTension_Unit_Conversion(
            -0.5 * LX * LZ * (Pyy[i] - 0.5 * (Pxx[i] + Pzz[i]))
        ) for i in range(len(Pxx))
    ]

    mean_LT = average(LT)
    with open("Mean_LT.dat", "w") as fo:
        fo.write(f"{mean_LT}\n")
    return LT, mean_LT

def Plot_Lt(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE):
    LT, _ = Line_Tension(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE)
    LT = LT[-680:]
    T = Time(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE)
    T_ns = [t / 1000 for t in T[-680:]]  # Convert to ns

    mean_LT = statistics.mean(LT)
    sd_LT = statistics.stdev(LT)
    LT_min, LT_max = min(LT), max(LT)
    T_min, T_max = min(T_ns), max(T_ns)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Time series plot
    axes[0].plot(T_ns, LT, 'o-')
    axes[0].set_xlabel("Time (ns)")
    axes[0].set_ylabel("Line Tension (pN)")
    axes[0].set_xticks(np.arange(T_min, T_max + 5, step=20))
    step = max((LT_max - LT_min) // 10, 1)
    axes[0].set_yticks(np.arange(LT_min, LT_max + 1, step=step))

    # Histogram / Gaussian fit
    axes[1].plot(LT, norm.pdf(LT, mean_LT, sd_LT), 'bs')
    axes[1].axvline(mean_LT, color='grey', linestyle='--', linewidth=2, label=f"Mean: {mean_LT:.2f}, std: {sd_LT:.2f}")
    axes[1].axvline(mean_LT - sd_LT, color='black', linestyle='--', linewidth=2)
    axes[1].axvline(mean_LT + sd_LT, color='black', linestyle='--', linewidth=2)
    axes[1].legend(loc='upper right')

    plt.tight_layout()
    plt.savefig(f"LT_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.png")
    plt.show()

def AAA(Sample_Num, Sample_Pass, epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE):
    filename = f"Pxx_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.xvg"
    if os.path.isfile(filename):
        with open(filename) as f0:
            last_line = f0.readlines()[-1]
            if "252000.000000" in last_line:
                Plot_Lt(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE)
                _, LT_mean = Line_Tension(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE)
                return LT_mean
            else:
                return 10000
    else:
        return 20000

# -----------------------------
# Run
# -----------------------------
LT = AAA(Sample_Num, Sample_Pass, epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE)
print(LT)
