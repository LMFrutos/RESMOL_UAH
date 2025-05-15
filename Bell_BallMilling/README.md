# Bell_BallMilling_Mechanochemistry
Program for determination of activation energy variation in Ball Milling mechanochemistry processes


Compute_W is a lightweight Python 3 program that builds the work surface
W(φ,θ) associated with applying an isotropic external pressure to a molecule.
The code is fully vectorised with NumPy, performs no graphical output, and communicates exclusively through text files.

├── compute_W.py   ← main script
├── R_min.txt      ← starting geometry (local minimum)
├── R_TS.txt       ← target geometry (e.g. transition state)
├── mass.txt       ← external pressure and atomic masses
├── output.txt     ← numerical report          (created at runtime)
└── W.mtx          ← W (φ, θ) in Matrix-Market (created at runtime)



Input Formats
| File                           | Line # | Content                                                    |
| ------------------------------ | ------ | ---------------------------------------------------------- |
| **R\_min.txt** / **R\_TS.txt** | 1      | `N` – number of atoms (integer)                            |
|                                | 2-N+1  | `x  y  z` – Cartesian coordinates in Å (floats)            |
| **mass.txt**                   | 1      | `P_ext_GPa` – external pressure in **GPa** (float)         |
|                                | 2      | `N` – number of atoms (integer, must match the geometries) |
|                                | 3-N+2  | `m_i` – atomic masses in **atomic units** (floats)         |

Outputs
| File           | Description                                                                     |
| -------------- | ------------------------------------------------------------------------------- |
| **output.txt** | Human-readable summary: pressure, sphere radius, forces, ⟨ΔEₐ⟩, min/max W, etc.  |
| **W\.mtx**     | Dense Matrix-Market file containing W (φ, θ) for further plotting or analysis   |


Quick Start
pip install numpy scipy
python compute_W.py




Algorithm Overview
Load and validate input
The script reads the two geometries and the mass/pressure file, ensuring that all three agree on the atom count N.

Remove global translation
Each geometry is shifted so that its mass-weighted centre of mass (COM) sits at the origin.

Remove global rotation
A mass-weighted Kabsch alignment rotates the target geometry onto the reference geometry, minimising the RMSD.

Internal displacements
ΔR = R_TS_aligned – R_min_c contains the true per-atom displacements, free of overall translation and rotation.

Bounding sphere
The radius R_sphere is the largest distance from the origin found in either geometry.

External force from pressure
The isotropic pressure P_ext is converted to a total force
F_ext = P_ext_GPa · 0.01 · π · R_sphere² (in nanonewtons), then divided equally among atoms.

Grid evaluation of W (φ, θ)
A 200 × 100 grid covers the unit sphere. For each direction n(φ, θ) the script finds the largest radial projection of the molecule onto the plane normal to n, scales the external force accordingly, and accumulates the work contribution.

Integration and unit conversion
The program reproduces the established conversion to kcal · mol⁻¹ via the constant
1.88973 × 627.5 / 82.387.

Write results
Key numbers go to output.txt; the complete W matrix is saved as W.mtx for direct import into Matlab, Octave, Python, etc.

Dependencies
Python ≥ 3.9

NumPy ≥ 1.21

SciPy ≥ 1.9 (used only for scipy.io.mmwrite)

A minimal installation file requirements.txt is included:
numpy>=1.21
scipy>=1.9

Set-up example:
python -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate
pip install -r requirements.txt


Contributing
Pull requests and issue reports are welcome.
Author: Luis Manuel Frutos
