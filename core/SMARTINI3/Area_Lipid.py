import pandas as pd
import numpy as np
import os
import shutil
import subprocess
import multiprocessing
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import time

start_time = time.time()

outdir = Path('/home/uni08/soleimani/RUNS/POSTDOC/Martini_HD/Rim_Pore_Location/')
indir = Path('/home/uni08/soleimani/RUNS/POSTDOC/Martini_HD/Rim_Pore_Location/APL_Post_Analysis.py/')

NUMBER_FRAMES = 99

class ReadGro:
    line_len: int = 45

    def __init__(self, fname: str):
        self.gro_data = self.read_gro(fname)

    def read_gro(self, fname: str) -> pd.DataFrame:
        processed_line = []
        with open(fname, 'r', encoding='utf8') as f_r:
            for counter, line in enumerate(f_r):
                line = line.rstrip()
                if len(line) != self.line_len:
                    self._process_header_tail(line, counter)
                else:
                    processed_line.append(self._process_line(line))
        return pd.DataFrame(processed_line)

    @staticmethod
    def _process_line(line: str) -> dict:
        return {
            'residue_number': int(line[0:5]),
            'residue_name': line[5:10].strip(),
            'atom_name': line[10:15].strip(),
            'atom_id': int(line[15:20]),
            'x': float(line[20:28]),
            'y': float(line[28:36]),
            'z': float(line[36:44])
        }

    def _process_header_tail(self, line: str, counter: int):
        if counter == 0:
            self.title = line
        elif counter == 1:
            self.number_atoms = int(line)
        elif counter == self.number_atoms + 2:
            self.pbc_box = line

class APL_ANALYSIS:
    def __init__(self, membrane_LX=10, membrane_LY=8, mesh_resolution=30):
        self.membrane_LX = membrane_LX
        self.membrane_LY = membrane_LY
        self.mesh_resolution = mesh_resolution
        self.membrane_area = membrane_LX * membrane_LY

    def _get_xy_grid(self):
        mesh_size_X = self.membrane_LX / self.mesh_resolution
        mesh_size_Y = self.membrane_LY / self.mesh_resolution
        x_mesh, y_mesh = np.meshgrid(
            np.arange(0, self.membrane_LX, mesh_size_X),
            np.arange(0, self.membrane_LY, mesh_size_Y)
        )
        Mesh_NUMBER = x_mesh.size
        return x_mesh, y_mesh, mesh_size_X, mesh_size_Y, Mesh_NUMBER

    @staticmethod
    def write_gromacs_gro(gro_data: pd.DataFrame, output_directory: Path, filename: str, pbc_box=None, title=None):
        output_file_path = output_directory / filename
        with open(output_file_path, 'w', encoding='utf8') as gro_file:
            if title: gro_file.write(f'{title}\n')
            gro_file.write(f'{len(gro_data)}\n')
            for _, row in gro_data.iterrows():
                gro_file.write(
                    f'{row["residue_number"]:>5}{row["residue_name"]:<5}{row["atom_name"]:>5}'
                    f'{row["atom_id"]:>5}{row["x"]:8.3f}{row["y"]:8.3f}{row["z"]:8.3f}\n'
                )
            if pbc_box: gro_file.write(f'{pbc_box}\n')

    @classmethod
    def process_mesh(cls, x_mesh, y_mesh, mesh_size_X, mesh_size_Y, mesh_resolution, xyz_i, max_z, min_z, frame):
        atom_types = ['NC3', 'PO4', 'GL1', 'GL2', 'C1A', 'C2A', 'C3A', 'C1B', 'C2B', 'C3B']
        selected_atoms_info = {}

        x_edges = np.arange(0, x_mesh.max() + mesh_size_X, mesh_size_X)
        y_edges = np.arange(0, y_mesh.max() + mesh_size_Y, mesh_size_Y)

        for i in range(mesh_resolution):
            for j in range(mesh_resolution):
                x_min, x_max = x_edges[i], x_edges[i+1]
                y_min, y_max = y_edges[j], y_edges[j+1]

                mask_type = np.isin(xyz_i[:, 2], atom_types)
                mask_xy = (xyz_i[:, 4] >= x_min) & (xyz_i[:, 4] < x_max) & (xyz_i[:, 5] >= y_min) & (xyz_i[:, 5] < y_max)
                mask_z = (xyz_i[:, 6] >= min_z) & (xyz_i[:, 6] <= max_z)
                indices = np.where(mask_type & mask_xy & mask_z)[0]

                selected_atoms_info[(i, j)] = xyz_i[indices]

        # Create folders
        for idx, atoms in selected_atoms_info.items():
            folder_index = idx[0] + idx[1] * mesh_resolution
            folder_path = outdir / f"Eolder_{folder_index}"
            folder_path.mkdir(parents=True, exist_ok=True)
            if len(atoms) == 0: continue
            df = pd.DataFrame(atoms, columns=['residue_number', 'residue_name', 'atom_name', 'atom_id', 'x', 'y', 'z'])
            APL_ANALYSIS.write_gromacs_gro(df, folder_path, f"grid_{frame}.gro", pbc_box="42.44358 44.57039 18.80638",
                                           title=f"Atoms in grid {idx} frame {frame}")

        return 0

def process_frame(frame):
    gro_data = ReadGro(f"conf{frame}.gro").gro_data
    xyz_i = gro_data[['residue_number','residue_name','atom_name','atom_id','x','y','z']].values
    mesh_gen = APL_ANALYSIS()
    x_mesh, y_mesh, mesh_size_X, mesh_size_Y, Mesh_NUMBER = mesh_gen._get_xy_grid()
    mesh_gen.process_mesh(x_mesh, y_mesh, mesh_size_X, mesh_size_Y, mesh_gen.mesh_resolution, xyz_i, max_z=1000, min_z=-1000, frame=frame)

if __name__ == "__main__":
    frames = range(NUMBER_FRAMES)
    with multiprocessing.Pool(min(multiprocessing.cpu_count(), NUMBER_FRAMES)) as pool:
        pool.map(process_frame, frames)

    mesh_gen = APL_ANALYSIS()
    x_mesh, y_mesh, mesh_size_X, mesh_size_Y, Mesh_NUMBER = mesh_gen._get_xy_grid()

    # Copy scripts to mesh folders
    for i in range(Mesh_NUMBER):
        folder = outdir / f'Eolder_{i}'
        folder.mkdir(exist_ok=True)
        shutil.copy2(indir / "APL_Post_Analysis.py", folder)
        shutil.copy2(indir / "slrun_Julich_APL_Post_Analysis", folder)
        shutil.copy2(indir / "Execute.sh", folder)

    def process_function(NOM):
        folder = outdir / f'Eolder_{NOM}'
        subprocess.call('chmod +x ./Execute.sh', shell=True, cwd=folder)
        subprocess.call(['bash', 'Execute.sh', str(NUMBER_FRAMES), str(mesh_size_X), str(mesh_size_Y), str(mesh_gen.mesh_resolution)], cwd=folder)

    with multiprocessing.Pool(min(multiprocessing.cpu_count(), Mesh_NUMBER)) as pool:
        pool.map(process_function, range(Mesh_NUMBER))

    # Collect results
    My_leaflet_1 = []
    for i in range(Mesh_NUMBER):
        file_path = outdir / f'Eolder_{i}' / 'result_ave.txt'
        raw_data = pd.read_csv(file_path, delim_whitespace=True, header=None, nrows=5)
        My_leaflet_1.append(raw_data.iloc[0,0])

    # Display
    def Display(input_list, layer):
        x, y = np.meshgrid(np.linspace(0, mesh_gen.membrane_LX, mesh_gen.mesh_resolution + 1),
                           np.linspace(0, mesh_gen.membrane_LY, mesh_gen.mesh_resolution + 1))
        plt.figure(figsize=(10,10))
        plt.pcolormesh(x, y, np.zeros_like(x), cmap='YlGnBu')
        plt.scatter(x, y, color='black', s=2)
        segs1 = np.stack((x, y), axis=2)
        segs2 = segs1.transpose(1,0,2)
        plt.gca().add_collection(LineCollection(segs1))
        plt.gca().add_collection(LineCollection(segs2))

        random_numbers_grid = np.array(input_list).reshape(mesh_gen.mesh_resolution, mesh_gen.mesh_resolution)
        cell_size = x[0,1] - x[0,0]
        for j in range(mesh_gen.mesh_resolution):
            for i in range(mesh_gen.mesh_resolution):
                plt.text(x[i,j]+cell_size/2, y[i,j]+cell_size/2, f"{random_numbers_grid[i,j]:.3f}", ha='center', va='center', fontsize=8)
        plt.savefig(f'Plot_leaflet_NUM_{layer}.png')
        plt.close()

    Display(My_leaflet_1, layer=1)

    # Execution time
    execution_time_minutes = (time.time() - start_time) / 60
    with open('execution_time.txt', 'w') as time_file:
        time_file.write(f"Execution time: {execution_time_minutes} minutes\n")

print("END OF THE CODE! GOOD LUCK2 ...")
