🧠 SMARTINI3: Systematic Parametrization of Realistic Multi-Scale Membrane Models using unsupervised learning 

SMARTINI3 is a minimal yet realistic ultra-coarse-grained (UCG) membrane model, systematically parameterized using unsupervised learning and multi-objective evolutionary algorithms.
It is designed to retain the biophysical fidelity of lipid membranes while enabling large-scale, efficient simulations, especially those involving Martini membrane proteins.



---

## 🌟 Key Features

| # | Feature | Description |
|---|----------|--------------|
| 🧩 **1. Ultra-Coarse-Grained Representation** | Captures essential membrane physics with minimal particles per lipid. |
| 💧 **2. Implicit Solvent Model** | Removes explicit water particles while preserving hydrophobic and hydrophilic balance. |
| ⚙️ **3. Systematic Parametrization via Genetic Algorithms** | Uses machine learning–driven optimization (Open BEAGLE) to fine-tune parameters for accuracy and transferability. |
| 🧬 **4. Realistic Membrane Behavior** | Reproduces essential lipid bilayer properties such as thickness, area per lipid, and bending rigidity. |
| 🧠 **5. Martini Compatibility** | Seamlessly integrates with the Martini coarse-grained model for hybrid or full Martini simulations. |
| 💻 **6. GROMACS Support** | Fully compatible with GROMACS, enabling high-performance molecular dynamics simulations. |
| ⚡ **7. High Computational Efficiency** | Achieves up to 35% improved performance compared to comparable coarse-grained models. |

---
## ⚙️ Installation

Clone the repository:

```bash
git clone https://github.com/soleimani2020/MACHINE_LEARNING_SMARTINI3.git
cd MACHINE_LEARNING_SMARTINI3
```

---

## 💡 Parametrization Framework

SMARTINI3 uses a **genetic algorithms (GAs) evolutionary computation framework** for parameter optimization.

### 🔬 Genetic Algorithms Integration

The parameter optimization leverages genetic algorithms (GAs), enabling systematic and efficient exploration of the parameter space.
