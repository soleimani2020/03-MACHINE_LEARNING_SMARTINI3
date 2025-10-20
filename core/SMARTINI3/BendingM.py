#!/usr/bin/env python3
import os
import sys
import numpy as np
import pandas as pd
from itertools import islice
from scipy.optimize import curve_fit
from scipy import ndimage
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy.fft as FF

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

# ---------------------------
# Utility functions
# ---------------------------
def average(lst):
    return np.mean(lst)

def read_gro_box(it):
    """Read simulation box dimensions from GRO file."""
    filename = f"conf_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}_{it}.gro"
    with open(filename) as f:
        lx, ly = map(float, f.readlines()[-1].split())
    return lx, ly

def box_matrix(final_TS, initial_TS=0, Interval=1):
    """Build array of box dimensions for all frames."""
    result = np.zeros([final_TS, 2], dtype=float)
    for it in range(initial_TS, final_TS, Interval):
        lx, ly = read_gro_box(it)
        result[it, :] = [lx, ly]
    return result, np.array(result)

def read_snapshot1(it, num_particles=4896):
    """Read particle coordinates from snapshot (first part)."""
    filename = f"conf_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}_{it}.gro"
    with open(filename) as f:
        lines = list(islice(f, 2, 2 + num_particles))
        data = [line.split() for line in lines]
        df = pd.DataFrame(data, columns=range(len(data[0])))
        df = df.astype({3:'float64', 4:'float64', 5:'float64'})
    return df[3].values, df[4].values, df[5].values

def trajectory_matrix(final_TS, initial_TS=0, Interval=1, num_particles=4896):
    """Build trajectory array for all frames."""
    result = np.zeros([final_TS, num_particles, 3], dtype=float)
    for it in range(initial_TS, final_TS, Interval):
        x, y, z = read_snapshot1(it)
        result[it, :, 0] = x
        result[it, :, 1] = y
        result[it, :, 2] = z
    return result, np.array(result)

# ---------------------------
# Planar modal analysis
# ---------------------------
def planar_modal_analysis(traj, box_dim, n_modes, interp_method="linear"):
    n_steps, n_particles = traj.shape[:2]
    qL_max = np.pi * (n_modes - 1)
    dqL = qL_max / ((n_modes - 1) // 2)

    sampled_qL, sampled_hpow, sampled_hpow_A, sampled_h, sampled_abs = [], [], [], [], []
    Result0 = np.zeros([n_modes, n_modes, n_steps], dtype=float)

    for i in range(n_steps):
        pos = traj[i]
        l_x, l_y = box_dim[i, :2]
        xx, yy = np.linspace(0, l_x, n_modes), np.linspace(0, l_y, n_modes)
        x_2D, y_2D = np.meshgrid(xx, yy)

        # Interpolate height field
        x_, y_, h_ = pos[:, 0], pos[:, 1], pos[:, 2]
        h_2D = griddata((x_, y_), h_, (x_2D, y_2D), method=interp_method)
        h_2D -= np.mean(h_2D)

        # FFT
        hff = FF.fftshift(FF.fft2(h_2D, [n_modes, n_modes])) / (n_modes**2)
        h_pow = np.abs(hff)**2

        freq_x = FF.fftshift(FF.fftfreq(n_modes, 1.0)) * n_modes
        freq_y = FF.fftshift(FF.fftfreq(n_modes, 1.0)) * n_modes
        freq_x_2D, freq_y_2D = np.meshgrid(freq_x, freq_y)
        qL = 2 * np.pi * np.sqrt((freq_x_2D / l_x) ** 2 + (freq_y_2D / l_y) ** 2)

        _ind = np.clip(qL, 0, 100)
        index = np.unique(qL)

        sampled_qL.append(ndimage.mean(qL, _ind, index)[1:-1])
        sampled_h.append(ndimage.mean(np.real(hff), _ind, index)[1:-1])
        sampled_hpow.append(ndimage.mean(np.real(h_pow), _ind, index)[1:-1])
        sampled_hpow_A.append(ndimage.mean(np.real(h_pow) / (l_x*l_y), _ind, index)[1:-1])
        sampled_abs.append(ndimage.mean(h_pow, _ind, index)[1:-1])

    return Result0, np.transpose(np.array(sampled_qL)), np.transpose(np.array(sampled_h)), \
           np.transpose(np.array(sampled_hpow)), np.transpose(np.array(sampled_hpow_A)), \
           np.transpose(np.array(sampled_abs))

# ---------------------------
# Bending rigidity calculation
# ---------------------------
def fluctuation_spectrum_kappa():
    xvg_file = f"LX_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.xvg"
    if not os.path.isfile(xvg_file):
        return 300000

    with open(xvg_file) as f:
        if "252000.000000" not in f.readlines()[-1]:
            return 200000

    # Parameters
    initial_TS, final_TS, Interval, n_modes = 0, 1680, 1, 200
    KT = 1
    a0 = 0.66
    Num_Bead = 816

    # Run analysis
    traj = trajectory_matrix(final_TS, initial_TS, Interval)[0]
    box_dim = box_matrix(final_TS, initial_TS, Interval)[0]
    _, sampled_qL, _, _, _, sampled_abs = planar_modal_analysis(traj, box_dim, n_modes)

    Myave_q = np.mean(sampled_qL, axis=1)
    Myave_abs = np.mean(sampled_abs, axis=1)
    if Myave_abs[-1] >= 0.0009:
        return 100000

    # Fit q^-4 to determine Kappa
    def model(q, Kappa):
        return KT / (Num_Bead * a0 * (Kappa * q**4))

    fit_Kappa, _ = curve_fit(model, Myave_q, Myave_abs, maxfev=10000)
    Kappa = fit_Kappa[0]

    # Plot
    plt.figure(figsize=(6,4))
    plt.plot(Myave_q, Myave_abs, 'bo', label=f"Kappa: {Kappa:.3f}")
    t = np.linspace(Myave_q.min(), Myave_q.max())
    plt.plot(t, model(t, Kappa))
    plt.xlabel('q [1/nm]')
    plt.ylabel('<|u(q)|^2>')
    plt.legend(loc='upper right')
    plt.savefig(f"BR_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.png")

    plt.figure(figsize=(6,4))
    plt.plot(Myave_q, Myave_abs, 'bo', label=f"Kappa: {Kappa:.3f}")
    plt.loglog(t, model(t, Kappa))
    plt.xlabel('q [1/nm]')
    plt.ylabel('<|u(q)|^2>')
    plt.legend(loc='upper right')
    plt.savefig(f"BR_loglog_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.png")

    return Kappa

# ---------------------------
# Main
# ---------------------------
A = 0  # Can call: A = fluctuation_spectrum_kappa()
print(A)
