#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import islice
import sys

# -----------------------------
# Core functions for potential/force
# -----------------------------

def get_abc(alpha, r1, rc):
    A = -alpha * ((alpha + 4) * rc - (alpha + 1) * r1) / (rc**(alpha + 2) * (rc - r1)**2)
    B = alpha * ((alpha + 3) * rc - (alpha + 1) * r1) / (rc**(alpha + 2) * (rc - r1)**3)
    C = 1 / rc**alpha - (A / 3) * (rc - r1)**3 - (B / 4) * (rc - r1)**4
    return A, B, C


def get_s(r, alpha, r1, rc):
    A, B, _ = get_abc(alpha, r1, rc)
    return A * (r - r1)**2 + B * (r - r1)**3


def get_switched_force(r, alpha, r1, rc):
    unswitched = alpha / r**(alpha + 1)
    switched = unswitched + get_s(r, alpha, r1, rc)
    switched[r < r1] = unswitched[r < r1]
    switched[rc <= r] = 0.0
    switched[r < 0.04] = 0.0
    return switched


def get_switched_potential(r, alpha, r1, rc):
    unswitched = 1 / r**alpha
    A, B, C = get_abc(alpha, r1, rc)
    result = unswitched - (A / 3) * (r - r1)**3 - (B / 4) * (r - r1)**4 - C
    result[r > rc] = 0.0
    result[r < r1] = unswitched[r < r1] - C
    result[r < 0.04] = 0.0
    return result


def v_vdw(r, A, C):
    return A / r**6 - C / r**0.5


def F_vdw(r, A, C):
    return 6 * A / r**7 - 0.5 * C / r**1.5


def v_vdw_switch(r, r1, rc, A, C):
    vs6 = get_switched_potential(r, alpha=6, r1=r1, rc=rc)
    vshalf = get_switched_potential(r, alpha=0.5, r1=r1, rc=rc)
    return A * vs6, C * vshalf, A * vs6 - C * vshalf


def F_vdw_switch(r, r1, rc, A, C):
    Fs6 = get_switched_force(r, alpha=6, r1=r1, rc=rc)
    Fshalf = get_switched_force(r, alpha=0.5, r1=r1, rc=rc)
    return Fs6 * A, Fshalf * C, A * Fs6 - C * Fshalf


def Parameter_TT(epsilon, rmin):
    A = (epsilon / 5.5) * 0.5 * rmin**6
    C = (epsilon / 5.5) * 6 * rmin**0.5
    return A, C


def v_vdw_HH_switch(r, r1, rc, A):
    vs6 = get_switched_potential(r, alpha=6, r1=r1, rc=rc)
    return 0.4 * A * vs6


def F_vdw_HH_switch(r, r1, rc, A):
    Fs6 = get_switched_force(r, alpha=6, r1=r1, rc=rc)
    return 0.4 * A * Fs6


def Parameter_HH(epsilon, rmin):
    return 0.5 * (epsilon / 5.5) * rmin**6


# -----------------------------
# Table generation functions
# -----------------------------

def Tables(epsilon, rmin, Max_distance, r1, rc, bFC, aFC, Friction, NFE):
    A, C = Parameter_TT(epsilon, rmin)
    r = np.arange(0.0, Max_distance, 0.002)

    # Potentials and forces
    analytic_noswitch_potential = v_vdw(r, A, C)
    H, G, analytic_switch_potential = v_vdw_switch(r, r1, rc, A, C)
    H_force, G_force, analytic_switch_force = F_vdw_switch(r, r1, rc, A, C)
    noswitch_force = F_vdw(r, A, C)

    # Plot potentials
    fig, ax = plt.subplots()
    ax.plot(r, analytic_noswitch_potential, 'k.-', label='Normal potential')
    ax.plot(r, analytic_switch_potential, 'g--', linewidth=3, label='Shifted potential')
    ax.axvline(r1, color='black', linestyle='--', linewidth=2)
    ax.axvline(rc, color='black', linestyle='--', linewidth=2)
    ax.axvline(rmin, color='red', linestyle='--', linewidth=2)
    ax.axhline(0, color='black', linestyle='--', linewidth=2)
    ax.set_xlim(0, 3)
    ax.set_ylim(-8, 4)
    ax.set_xlabel('r (nm)')
    ax.set_ylabel('$V_{LJ}$')
    ax.legend(fancybox=True)
    plt.grid(True)
    plt.savefig(f"Plot_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.png")

    # Zoomed plot
    fig, ax = plt.subplots()
    ax.plot(r, analytic_switch_potential, 'g--', linewidth=3, label='Shifted potential')
    ax.axvline(r1, color='black', linestyle='--', linewidth=2)
    ax.axvline(rc, color='black', linestyle='--', linewidth=2)
    ax.set_xlim(r1 - 0.2, rc + 0.2)
    ax.set_ylim(-0.25, 0.25)
    ax.set_xlabel('r (nm)')
    ax.set_ylabel('$V_{LJ}$')
    ax.legend(fancybox=True)
    plt.grid(True)
    plt.savefig(f"Plot_Zoom_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.png")

    # Electrostatic interactions table
    Distance = np.arange(0.0, Max_distance, 0.002)
    F = [0 if d < 0.04 else 1 / d for d in Distance]
    Minus_F_Prime = [0 if d < 0.04 else -(1 / d)**2 for d in Distance]
    G_vals = [0] * len(Distance)
    Minus_G_Prime = [0] * len(Distance)
    A_HH = Parameter_HH(epsilon, rmin)
    H_vals = v_vdw_HH_switch(Distance, r1, rc, A_HH)
    Minus_H_Prime_vals = F_vdw_HH_switch(Distance, r1, rc, A_HH)

    table_data = np.column_stack((Distance, F, Minus_F_Prime, G_vals, Minus_G_Prime, H_vals, Minus_H_Prime_vals))
    np.savetxt(
        f'table-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE}.xvg',
        table_data, delimiter=' ', fmt='%2f %8f %13f %12f %12f %9f %8f'
    )


def Tables_p(epsilon, rmin, Max_distance, r1, rc, bFC, aFC, Friction, NFE):
    Distance = np.arange(0.0, Max_distance, 0.002)
    zeros = [0] * len(Distance)
    F_vals = [0 if d < 0.04 else 1 / d for d in Distance]
    Minus_F_Prime = [0 if d < 0.04 else -(1 / d)**2 for d in Distance]

    table_data = np.column_stack((Distance, F_vals, Minus_F_Prime, zeros, zeros, zeros, zeros))
    np.savetxt(
        f'Ptable-{epsilon}-{rmin}-{r1}-{rc}-{bFC}-{aFC}-{Friction}-{NFE}.xvg',
        table_data, delimiter=' ', fmt='%2f %8f %13f %12f %12f %9f %8f'
    )
