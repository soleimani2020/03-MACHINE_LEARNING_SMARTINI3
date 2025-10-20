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

SKIPROW = 26

def Density_Profile(
    epsilon: str, rmin: str, r1: str, rc: str,
    bFC: str, aFC: str, Friction: str, NFE: str, Fitness: str
) -> np.ndarray | None:
    """
    Plot density profiles from GROMACS density xvg file.

    Saves a figure 'DP_<params>.png' if file exists, otherwise prints an error.
    Returns the loaded density array if successful.
    """
    filename = f"density_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}.xvg"

    if not os.path.isfile(filename):
        print(f"Simulation with parameters {epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE} "
              f"faced error: density file not found.")
        return None

    # Load data once
    data = np.loadtxt(filename, skiprows=SKIPROW)

    # Plot
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(data[:, 0], data[:, 1], '-o', label="System")
    ax.plot(data[:, 0], data[:, 2], '-o', label="T")
    ax.plot(data[:, 0], data[:, 3], '-o', label="H")

    ax.set_xlabel('Distance (Å)', fontsize=12)
    ax.set_ylabel('g(r)', fontsize=12)
    ax.legend(loc='upper right')
    ax.set_title(f"Density Profile (Fitness={Fitness})")

    out_png = f"DP_{epsilon}_{rmin}_{r1}_{rc}_{bFC}_{aFC}_{Friction}_{NFE}_{Fitness}.png"
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Density profile saved to {out_png}")

    return data


if __name__ == "__main__":
    dp_data = Density_Profile(epsilon, rmin, r1, rc, bFC, aFC, Friction, NFE, Fitness)
    if dp_data is not None:
        print(f"Loaded density data shape: {dp_data.shape}")
