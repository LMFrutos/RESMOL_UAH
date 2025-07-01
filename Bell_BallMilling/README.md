Bell-BallMilling-Mechanochemistry
A Python 3 program to compute the change in activation energy for chemical processes under Ball Milling conditions. The core script, compute_W.py, generates the work surface W(phi, theta) for a molecule subjected to external isotropic pressure. NumPy and SciPy are used for vectorized calculations, producing both numerical reports and graphical outputs.

Repository Structure
lua
Copiar
Editar
Bell_BallMilling/
├── Bell_BallMilling.py
├── compute_W.py
├── pyproject.toml
├── requirements.txt
│
├── example_files/
│   ├── input_files/
│   │   ├── R_min.txt
│   │   ├── R_TS.txt
│   │   └── input.txt
│   ├── output_files/
│   │   ├── output.txt
│   │   ├── W.mtx
│   │   └── W_spherical.png
│   └── README_example_info.txt
Installation
It is strongly recommended to use a virtual environment (venv or similar) to avoid conflicts.

Method 1: Standard Installation (Recommended, CLI support)
bash
Copiar
Editar
git clone https://github.com/LMFrutos/RESMOL_UAH.git
cd RESMOL_UAH/Bell_BallMilling
python -m venv .venv
source .venv/bin/activate         # On Windows: .venv\Scripts\activate
pip install -e .
This installs the package in editable mode and adds the CLI command bmm-compute.

Method 2: Development Environment
bash
Copiar
Editar
git clone https://github.com/LMFrutos/RESMOL_UAH.git
cd RESMOL_UAH/Bell_BallMilling
python -m venv .venv
source .venv/bin/activate         # On Windows: .venv\Scripts\activate
pip install -r requirements.txt
Quick Start
Prepare your input files in example_files/input_files/:

R_min.txt, R_TS.txt, input.txt
(Check the example files and README_example_info.txt for guidance.)

Run the program from the main directory (Bell_BallMilling):

If installed with pip (Method 1):

bash
Copiar
Editar
bmm-compute
If running directly (Method 2):

bash
Copiar
Editar
python compute_W.py
The results will appear in example_files/output_files/:

output.txt: numeric report

W.mtx: raw W(phi, theta) data

W_spherical.png: 3D and 2D graphical representations

File Formats
Input files (example_files/input_files/)
File	Description
R_min.txt / R_TS.txt	Line 1: N (number of atoms). Lines 2–N+1: Z x y z (atomic number and Cartesian coordinates, in Å).
input.txt	Line 1: P_ext_GPa (float), external pressure in GPa
Line 2: ntheta (int), angular resolution
Line 3: YES or NO (whether to include van der Waals radii)

See example_files/README_example_info.txt for more details and sample values.

Output files (example_files/output_files/)
output.txt: summary (sphere radius, forces, mean delta Ea, min/max W values, etc.)

W.mtx: plain text; header with nphi, ntheta, followed by phi theta W(phi,theta) (not Matrix-Market format)

W_spherical.png: 2D/3D visual representation of W(phi, theta)

Algorithm Overview
Loads and validates the input geometry and settings.

Centers geometries using the center of mass and Kabsch rotation.

Calculates internal displacements (delta R).

Computes the spherical radius from the center, optionally adding van der Waals radii.

Applies isotropic external force, distributes over atoms.

Generates a spherical grid (phi, theta) and computes the work W in each direction.

Saves the results (output.txt, W.mtx, W_spherical.png).

Dependencies
Python >= 3.9

NumPy >= 1.21

SciPy >= 1.9

Matplotlib >= 3.5

Dependencies are installed automatically with pip install -e . or pip install -r requirements.txt.

Example
A complete example (input and output) is provided in the example_files/ folder.
Refer to README_example_info.txt in that directory for a detailed explanation of file formats and example contents.

Contributing
Pull requests and issues are welcome! Please:

Clearly describe changes or bugs.

Include sample files or data if relevant.

Add tests if you implement new features.

Author
Luis Manuel Frutos

Future improvements
More robust validation of input files

Additional statistical reporting

Export to standard formats (CSV, JSON)

Automated documentation (Sphinx, MkDocs)

Thank you for using Bell-BallMilling-Mechanochemistry!
