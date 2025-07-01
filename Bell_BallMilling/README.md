# Bell-BallMilling-Mechanochemistry Model

A program for the determination of the activation energy variation in Ball Milling mechanochemistry processes.

`compute_W.py` is a Python 3 program that builds the work surface $W(\phi,\theta)$ associated with applying an isotropic external pressure to a molecule. The code uses NumPy for vectorized calculations and generates both a numerical report and a graphical visualization of the results.

## Repository Structure

```
├── Bell_BallMilling.py   ← Main script
├── pyproject.toml        ← Package configuration file (for pip installation)
├── requirements.txt      ← Dependency file for development environment
├── R_min.txt             ← Starting geometry (local minimum)
├── R_TS.txt              ← Target geometry (e.g., transition state)
├── input.txt             ← Input parameters (pressure, resolution, etc.)
│
├── output.txt            ← (Created at runtime) Numerical report
├── W.mtx                 ← (Created at runtime) W(φ,θ) surface in coordinate format
└── W_spherical.png       ← (Created at runtime) Graphical visualization of the results
```

## Installation

Two installation methods are provided. In either case, using a [Python virtual environment](https://docs.python.org/3/library/venv.html) is strongly recommended.

### Method 1: Standard Installation (Recommended)

This method uses the `pyproject.toml` file to install the package and its dependencies, making the script available as a command-line tool.

1.  Clone the repository and navigate to the root directory.
2.  Create and activate a virtual environment:
    ```bash
    python -m venv .venv
    source .venv/bin/activate  # On Windows: .venv\Scripts\activate
    ```
3.  Install the package in editable mode using `pip`:
    ```bash
    pip install -e .
    ```
    The `-e` (editable) flag is optional but allows you to modify the source code without reinstalling.

### Method 2: Development Environment

This method is for users who prefer not to install the package system-wide but to run the script directly from the repository.

1.  Clone the repository and create a virtual environment (steps 1 and 2 from Method 1).
2.  Install the dependencies using `requirements.txt`:
    ```bash
    pip install -r requirements.txt
    ```

## Quick Start

1.  **Prepare input files**: Ensure the `R_min.txt`, `R_TS.txt`, and `input.txt` files are in the directory and have the correct format (see below).
2.  **Run the program**:
    * If you used **Method 1**, simply run the following command in your terminal:
        ```bash
        bmm-compute
        ```
    * If you used **Method 2**, run the Python script directly:
        ```bash
        python compute_W.py
        ```
3.  **Check the results**: The program will generate `output.txt`, `W.mtx`, and `W_spherical.png` in the same directory.

## File Formats

### Input Files

| File | Line(s) | Content |
| :--- | :--- | :--- |
| **`R_min.txt`** / **`R_TS.txt`** | 1 | `N` – number of atoms (integer) |
| | 2 to N+1 | `Z x y z` – **Atomic number** (integer) and Cartesian coordinates in Å (floats) |
| **`input.txt`** | 1 | `P_ext_GPa` – external pressure in **GPa** (float) |
| | 2 | `ntheta` – grid resolution for the θ direction (integer) |
| | 3 | `YES` or `NO` – indicates whether to add van der Waals radii (string) to the computation of molecular sphere radius (it takes into account the atomic radii, otherwise this is set to zero |

### Output Files

| File | Description |
| :--- | :--- |
| **`output.txt`** | Human-readable summary: input parameters, sphere radius, forces, ⟨ΔEₐ⟩, min/max W, etc. |
| **`W.mtx`** | A text coordinate file containing the work function values. **It does not follow the Matrix-Market format**. Its header contains the grid sizes (`nphi`, `ntheta`), followed by the `phi theta W(phi,theta)` data. |
| **`W_spherical.png`** | An image file with the visualization of $W(\phi,\theta)$, including a 3D spherical representation and a 2D rectangular plot. |

## Algorithm Overview

1.  **Load and Validate Input**: The script reads the two geometries and the `input.txt` file, ensuring that the atom count and atomic numbers match.
2.  **Remove Global Translation**: Each geometry is shifted so that its mass-weighted center of mass is at the origin.
3.  **Remove Global Rotation**: A mass-weighted Kabsch alignment rotates the `R_TS` geometry to optimally superimpose it onto `R_min`.
4.  **Internal Displacements**: $\Delta R = R_{TS,aligned} – R_{min,centered}$ is calculated to get the true per-atom displacements.
5.  **Bounding Sphere**: The radius $R_{sphere}$ is determined as the largest distance from the origin found in either geometry. Optionally (based on the flag in `input.txt`), the van der Waals radius of the farthest atom is added.
6.  **External Force**: The isotropic pressure $P_{ext}$ is converted into a total force $F_{ext}$ and distributed among the atoms.
7.  **Grid Evaluation of $W(\phi, \theta)$**: The script iterates over a spherical grid defined by `ntheta`. For each direction, the external force is scaled based on the molecule's cross-section, and its contribution to the work is calculated.
8.  **Write Results**: The numerical report (`output.txt`), the work function data (`W.mtx`), and the plot (`W_spherical.png`) are saved.

## Dependencies

-   Python ≥ 3.9
-   NumPy ≥ 1.21
-   SciPy ≥ 1.9
-   Matplotlib ≥ 3.5

These dependencies are handled automatically by the described installation methods.

## Contributing

Pull requests and issue reports are welcome.

**Author**: Luis Manuel Frutos
