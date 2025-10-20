#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import islice
from scipy.optimize import curve_fit
import subprocess

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


TEMP1, TEMPN, INTERVAL = 200, 420, 20

# -----------------------------
# Directory setup
# -----------------------------
base_dir = '/home/uni08/soleimani/RUNS/Thesis/Membrane_Project/PAPER/4000/'
run_dir = os.path.join(base_dir, 'RUNS0')
param_str = f"R-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE}"
param_dir = os.path.join(run_dir, param_str)
phase_dir = os.path.join(param_dir, 'Phase_Transition')
phase_temp_dir = os.path.join(phase_dir, 'T-')

# -----------------------------
# Helper Functions
# -----------------------------
def safe_number(values):
    return [float(v) if v not in ['nan', '-nan'] else 1e15 for v in values]


def read_enthalpy(folder, filename, last_line_value):
    enthalpy_path = os.path.join(folder, filename)
    log_path = os.path.join(folder, 'log_Enthelpy.txt')
    H_list = []

    if os.path.isfile(enthalpy_path):
        with open(enthalpy_path) as f:
            lines = f.read().splitlines()
            last_line = lines[-1]
            if str(last_line_value) in last_line and 'nan' not in last_line:
                if os.path.isfile(log_path):
                    with open(log_path) as f_log:
                        log_lines = list(islice(f_log, 6))
                        val = np.char.split(np.array(log_lines)[1])[1]
                        H_list.append(val[0])
                    os.rename(enthalpy_path, os.path.join(folder, f"{filename.split('.')[0]}_{folder.split('/')[-1]}.xvg"))
                    subprocess.run(f"cp {filename.split('.')[0]}_{folder.split('/')[-1]}.xvg ../", shell=True, cwd=folder)
                    return safe_number(H_list)
    with open("Report.dat", "a") as f:
        f.write(f"Problem with {filename} in folder {folder}\n")
    return [1e15] * ((TEMPN - TEMP1) // INTERVAL)


def collect_enthalpy(filename, last_line_value):
    folders = [os.path.join(phase_temp_dir, str(T)) for T in range(TEMP1, TEMPN, INTERVAL)]
    all_H = []
    for folder in folders:
        all_H.extend(read_enthalpy(folder, filename, last_line_value))
    return all_H


def save_list_to_file(filename, data):
    with open(filename, "w") as f:
        for d in data:
            f.write(f"{d}\n")


def snapshot_data():
    df = pd.read_csv("Result.txt", sep="\s+", header=None)
    return df[0].values, df[1].values


def compute_derivative(x, y):
    dy = np.diff(y) / np.diff(x)
    dx = [(x[i] + x[i+1]) / 2 for i in range(len(dy))]
    return np.array(dx), np.array(dy)


def find_transition_zone(deriv, x):
    if all(v == 0 for v in deriv[:4]):
        return 1000, [1000, 1010]
    idx = np.argmax(deriv)
    return (x[idx] + x[idx + 1]) / 2, [x[idx], x[idx + 1]]


def fit_phase_transition(x, H, expected_Tm):
    if H[0] == 1e15:
        Tm = 0
    elif expected_Tm == 210.1234:
        Tm = expected_Tm
    else:
        try:
            def model(x, A1, DeltaA, x0, k, c):
                return A1 + c * x + DeltaA / (1 + np.exp(k * (x - x0)))

            p0 = [np.mean(H), np.diff(H)[0], expected_Tm, 0, 0]
            popt, _ = curve_fit(model, x, H, p0=p0, maxfev=100000)
            x0 = popt[2]

            plt.figure(figsize=(6, 4))
            plt.plot(x, H, 'bo')
            x_model = np.arange(TEMP1, TEMPN, 0.1)
            plt.plot(x_model, model(x_model, *popt), 'r', label=f"T_M: {x0:.2f}")
            plt.axvline(x0, color='k', linestyle='--', linewidth=1)
            plt.xlabel('Temperature (K)')
            plt.ylabel('Enthalpy (kJ/mol)')
            plt.legend(loc='lower right')
            plt.savefig("Phase_Transition_Temperature.png")

            Tm = x0 if TEMP1 < x0 < TEMPN else expected_Tm
        except Exception:
            Tm = 0
    with open("T_m.dat", "w") as f:
        f.write(f"{Tm}\n")
    return Tm

# -----------------------------
# Main Workflow
# -----------------------------
H0 = collect_enthalpy("Enthalpy0.xvg", 36000)
H = collect_enthalpy("Enthalpy.xvg", 34950)
save_list_to_file("Enthalpy.dat", H)

T_list = list(range(TEMP1, TEMPN, INTERVAL))
save_list_to_file("Temp.dat", T_list)

with open('Result.txt', 'w') as writer:
    for t, h in zip(T_list, H):
        writer.write(f"{t}     {h}\n")

x, Enthalpy = snapshot_data()
xprime, yprime = compute_derivative(x, Enthalpy)

plt.figure(figsize=(6, 4))
plt.plot(xprime, yprime, "-o", label="1st-Derivative")
plt.xlabel('Temperature (K)')
plt.ylabel('dH/dT')
plt.show()

Expected_Tm, Expected_Tm_zone = find_transition_zone(yprime, x)
save_list_to_file("Expected_Tm.dat", [Expected_Tm])
save_list_to_file("Expected_Tm_Zone.dat", Expected_Tm_zone)

Tm = fit_phase_transition(x, Enthalpy, Expected_Tm)

subprocess.run(f"cp T_m.dat ../../", shell=True)
result_file = os.path.join(destination_folder2, f"RESULTS_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.txt")
subprocess.run(f"cat T_m.dat >> {result_file}", shell=True)

subprocess.run("python3 ./Plot_ALL_H.py", shell=True, cwd=phase_dir)
