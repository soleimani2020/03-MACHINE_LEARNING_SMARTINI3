#!/usr/bin/env python3
import random
import shutil
import subprocess
import logging
from pathlib import Path
from typing import Tuple, List

import numpy as np
from ypstruct import structure

import Table
import GROMACS_FILES

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Constants
DESTINATION_FOLDER0 = Path("/home/uni08/soleimani/RUNS/Thesis/Membrane_Project/PAPER/4000/")
TEMP1 = 200
TEMPN = 420
INTERVAL = 20

# Helper functions
def create_directory(path: Path):
    path.mkdir(parents=True, exist_ok=True)
    logging.debug(f"Directory created: {path}")

def move_files(source: Path, destination: Path):
    try:
        shutil.move(str(source), str(destination))
        logging.info(f"Moved: {source} to {destination}")
    except Exception as e:
        logging.error(f"Failed to move {source}: {e}")

def copy_files(source: Path, destination: Path):
    try:
        shutil.copy(str(source), str(destination))
        logging.info(f"Copied: {source} to {destination}")
    except Exception as e:
        logging.error(f"Failed to copy {source}: {e}")

def run_command(command: List[str], cwd: Path = None, timeout: int = None):
    try:
        subprocess.run(command, cwd=cwd, timeout=timeout, check=True)
        logging.info(f"Command executed: {' '.join(command)}")
    except subprocess.TimeoutExpired:
        logging.warning(f"Command timed out: {' '.join(command)}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed ({e.returncode}): {' '.join(command)}")

# Optimization functions
def optimize(x: np.ndarray, nfe: int, lst: list) -> Tuple[float, int, list]:
    nfe += 1
    epsilon, rmin, r1, rc, bFC, aFC, Friction = x

    # Generate GROMACS files
    GROMACS_FILES.MINIMIZATION(epsilon, rmin, r1, rc, bFC, aFC, Friction, nfe)
    GROMACS_FILES.NVT(epsilon, rmin, r1, rc, bFC, aFC, Friction, nfe)
    GROMACS_FILES.NPT(epsilon, rmin, r1, rc, bFC, aFC, Friction, nfe)
    GROMACS_FILES.FORCEFIELD(epsilon, rmin, r1, rc, bFC, aFC, Friction, nfe)
    GROMACS_FILES.TOPOLOGY(epsilon, rmin, r1, rc, bFC, aFC, Friction, nfe)
    Table.Tables(epsilon, rmin, rc + 1.2, r1, rc, bFC, aFC, Friction, nfe)

    # Create run directory
    destination_folder = DESTINATION_FOLDER0 / "RUNS0"
    create_directory(destination_folder)

    # Run simulation script
    command = [
        "bash", "Calculation002.sh",
        str(epsilon), str(rmin), str(r1), str(rc),
        str(bFC), str(aFC), str(Friction), str(nfe)
    ]
    run_command(command, cwd=destination_folder, timeout=1800)

    # Calculate cost
    cost = calculate_cost(destination_folder)
    return cost, nfe, lst

def calculate_cost(destination_folder: Path) -> float:
    return 0.0  # Placeholder for actual cost calculation

# Genetic algorithm operators : cross-over and mutation 
def single_point_crossover(p1: structure, p2: structure, cut: int) -> Tuple[structure, structure]:
    c1, c2 = p1.deepcopy(), p2.deepcopy()
    c1.position = np.append(p1.position[:cut], p2.position[cut:])
    c2.position = np.append(p2.position[:cut], p1.position[cut:])
    return c1, c2

def mutate(x: structure, mu: float, sigma: np.ndarray) -> structure:
    y = x.deepcopy()
    flag = np.random.rand(*x.position.shape) <= mu
    ind = np.argwhere(flag)
    y.position[ind] += sigma[ind] * np.random.randn(*ind.shape)
    return y

# Main 
if __name__ == "__main__":
    npop = 48  # number of initial population ( optimize according to your resources. I used GWDG HPC so even 280 was pretty fine !)
    pop = [structure() for _ in range(npop)]

    for individual in pop:
        individual.position = np.random.uniform(
            [1.0, 1.0, 1.5, 2.5, 100, 10.0, 2],
            [10.0, 1.5, 2.5, 3.0, 20000, 20000, 2]
        )
        individual.position[2] = individual.position[1] + random.uniform(0.5, 1.5)
        individual.position[3] = individual.position[2] + random.uniform(0.0, 0.5)

    bestsol = structure()
    bestsol.cost = np.inf
    max_iterations = 1

    for it in range(max_iterations):
        for individual in pop:
            individual.cost, _, individual.Objectives = optimize(individual.position, 0, [0]*7)
            if individual.cost < bestsol.cost:
                bestsol = individual.deepcopy()

        logging.info(f"Iteration {it}: Best Cost = {bestsol.cost}")

    with open("results.txt", "w") as f:
        f.write(f"Best Solution: {bestsol.position}\nBest Cost: {bestsol.cost}\n")
