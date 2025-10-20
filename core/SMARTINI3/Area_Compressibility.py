import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statistics
from itertools import islice
from scipy.stats import norm
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

Sample_Num = 24
Sample_Pass = 1691 - 24   # Number of samples
Lipid_Num = 816           # Number of lipids in one leaflet

def Average(lst):
    return np.mean(lst)

def read_xvg(column=0):
    filename = f"LX_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.xvg"
    if not os.path.isfile(filename):
        return None, 20000
    with open(filename) as f:
        lines = f.read().splitlines()
        if "252000.000000" not in lines[-1]:
            return None, 10000
    with open(filename) as f:
        next(islice(f, Sample_Num))  # skip initial lines
        data_lines = list(islice(f, Sample_Pass))
        data_array = np.array([line.split() for line in data_lines], dtype=float)
        return data_array[:, column].tolist(), None

def Time(Sample_Num, Sample_Pass):
    time_list, error_code = read_xvg(0)
    if error_code:
        return error_code
    return [t / 1000 for t in time_list]

def Membrane_BOX_LX(Sample_Num, Sample_Pass):
    lx, error_code = read_xvg(1)
    return lx if lx is not None else error_code

def Membrane_BOX_LY(Sample_Num, Sample_Pass):
    ly, error_code = read_xvg(2)
    return ly if ly is not None else error_code

def Membrane_Area(Sample_Num, Sample_Pass):
    LX = Membrane_BOX_LX(Sample_Num, Sample_Pass)
    LY = Membrane_BOX_LY(Sample_Num, Sample_Pass)
    if isinstance(LX, list) and isinstance(LY, list):
        area = np.array(LX) * np.array(LY)
        apl = area / Lipid_Num
        return area.tolist(), apl.tolist()
    return LX, LX  # propagate error codes

def Plot_APL(Sample_Num, Sample_Pass, epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE):
    MEM_BOX_LX = Membrane_BOX_LX(Sample_Num, Sample_Pass)
    MEM_BOX_LY = Membrane_BOX_LY(Sample_Num, Sample_Pass)
    
    if MEM_BOX_LX in [10000, 20000] or MEM_BOX_LY in [10000, 20000]:
        return MEM_BOX_LX
    
    _, APL = Membrane_Area(Sample_Num, Sample_Pass)
    APL = APL[-667:]
    T_ns = Time(Sample_Num, Sample_Pass)[-667:]
    
    mean_apl = statistics.mean(APL)
    sd_apl = statistics.stdev(APL)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Time series
    axes[0].plot(T_ns, APL, 'bo')
    axes[0].set_xlabel('Time / ns')
    axes[0].set_ylabel('APL ($nm^2$)')
    axes[0].legend(loc='upper right')
    
    # Histogram
    axes[1].plot(APL, norm.pdf(APL, mean_apl, sd_apl), 'bo')
    axes[1].axvline(mean_apl, color='grey', linestyle='--', linewidth=2,
                    label=f"Mean: {mean_apl:.3f}, std: {sd_apl:.3f}")
    axes[1].axvline(mean_apl - sd_apl, color='black', linestyle='--', linewidth=2)
    axes[1].axvline(mean_apl + sd_apl, color='black', linestyle='--', linewidth=2)
    axes[1].legend(loc='upper right')
    
    plt.savefig(f"APL_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.png")
    plt.close()
    
    return mean_apl

# Example usage
A = Plot_APL(Sample_Num, Sample_Pass, epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE)
print(A)

# Placeholder for area compressibility (functions can be refactored similarly)
B = 0
print(B)
