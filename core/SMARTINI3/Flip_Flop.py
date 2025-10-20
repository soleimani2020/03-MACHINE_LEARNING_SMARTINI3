import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from typing import List, Tuple


def count_lipids_cylindrical_x(
    universe: mda.Universe,
    frame_number: int,
    lipid_resnames: List[str] = ["DLPC"],
    headgroup_bead: str = "PO4",
    cylinder_axis: int = 0
) -> Tuple[int, int, np.ndarray, np.ndarray]:
    """
    Count inner and outer leaflet lipids in a cylindrical membrane segment.
    
    Returns:
        inner_count, outer_count, inner_indices, outer_indices
    """
    lipids = universe.select_atoms(f"resname {' '.join(lipid_resnames)}")
    headgroups = lipids.select_atoms(f"name {headgroup_bead}")

    # Coordinates in YZ plane (perpendicular to cylinder_axis)
    coords = headgroups.positions[:, [i for i in range(3) if i != cylinder_axis]]

    # Segment center of mass
    com = lipids.center_of_mass()

    # Compute radial distances in YZ plane
    distances = np.linalg.norm(coords - com[[1, 2]], axis=1)
    median_radius = np.median(distances)
    
    # Classify inner vs outer
    outer_mask = distances > median_radius
    outer_indices = np.where(outer_mask)[0]
    inner_indices = np.where(~outer_mask)[0]

    # Plot 3D classification
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')
    colors = np.zeros(len(coords))
    colors[outer_indices] = 1  # Outer = 1
    ax.scatter(com[0], com[1], com[2], color='green', s=50, label='COM')
    sc = ax.scatter(headgroups.positions[:, 0], headgroups.positions[:, 1], headgroups.positions[:, 2],
                    c=colors, cmap='bwr', s=10, alpha=0.6)

    # Add median radius circle in YZ plane
    theta = np.linspace(0, 2*np.pi, 100)
    y_circle = com[1] + median_radius * np.cos(theta)
    z_circle = com[2] + median_radius * np.sin(theta)
    x_circle = np.full_like(y_circle, com[0])
    ax.plot(x_circle, y_circle, z_circle, 'g--', linewidth=2, label='Median Radius')

    ax.set_title(f'Frame {frame_number}: Leaflet Classification (Blue=Inner, Red=Outer)')
    fig.colorbar(sc, ax=ax, shrink=0.5, label='Outer (1) / Inner (0)')
    plt.savefig(f"Leaflet_Classification_{frame_number}.png", dpi=300, bbox_inches='tight')
    plt.close(fig)

    return len(inner_indices), len(outer_indices), inner_indices, outer_indices


def analyze_first_n_frames(tpr_file: str, traj_file: str, n_frames: int = 10, segment_width: float = 100):
    universe = mda.Universe(tpr_file, traj_file)
    lipids = universe.select_atoms("resname DLPC")

    # Cylinder bounds
    x_min = np.min(lipids.positions[:, 0])
    x_max = np.max(lipids.positions[:, 0])
    num_segments = int((x_max - x_min) / segment_width)
    print(f"Cylinder length: {x_max - x_min:.2f} Å, Segments: {num_segments}")

    results = []

    for ts in tqdm(universe.trajectory[:n_frames], total=min(n_frames, len(universe.trajectory))):
        frame_data = []
        for segment in range(num_segments):
            seg_min = x_min + segment * segment_width
            seg_max = seg_min + segment_width
            segment_ag = lipids.select_atoms(f"prop x >= {seg_min} and prop x < {seg_max}")

            if len(segment_ag) == 0:
                frame_data.append({"Frame": ts.frame, "Segment": segment, "Inner": 0, "Outer": 0})
                continue

            inner, outer, _, _ = count_lipids_cylindrical_x(segment_ag, ts.frame)
            frame_data.append({"Frame": ts.frame, "Segment": segment, "Inner": inner, "Outer": outer})

        # Compute stats
        inner_vals = [row["Inner"] for row in frame_data]
        outer_vals = [row["Outer"] for row in frame_data]
        results.extend(frame_data)

        # Average row
        results.append({"Frame": ts.frame, "Segment": "Average",
                        "Inner": round(np.mean(inner_vals), 2),
                        "Outer": round(np.mean(outer_vals), 2)})

        # Sum row
        results.append({"Frame": ts.frame, "Segment": "Sum",
                        "Inner": int(np.sum(inner_vals)),
                        "Outer": int(np.sum(outer_vals))})

        # Blank row for spacing
        results.append({"Frame": "", "Segment": "", "Inner": "", "Outer": ""})

    df = pd.DataFrame(results)
    df.to_csv("inner_outer_counts.txt", sep='\t', index=False)
    print("Results saved to inner_outer_counts.txt")


if __name__ == "__main__":
    analyze_first_n_frames("nvt.tpr", "Traj_comp4.xtc", n_frames=10, segment_width=100)

