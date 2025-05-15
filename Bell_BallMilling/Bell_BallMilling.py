#!/usr/bin/env python3
"""
Ball-Milling-Bell-model
-----------------------
-----------------------

This program evaluates the work term W(φ,θ) applied to a molecular structure 
under external pressure.

General
----------------------------------------------
1.  The two geometries and the mass / pressure data are *not* hard-coded: they
    are loaded from external text files:
        • R_min.txt
        • R_TS.txt
        • mass.txt
   (formats are described below).
2.  Output files with the results:
        • output.txt  – main numbers and diagnostics
        • W.mtx       – W(φ,θ) in Matrix-Market format 


File formats
------------
R_min.txt and R_TS.txt
    Line 1         : integer N (number of atoms)
    Lines 2 … N+1  : x  y  z   (Cartesian coordinates, Å)

mass.txt
    Line 1         : float     external pressure, P_ext, in GPa
    Line 2         : integer   N (must match geometry files)
    Lines 3 … N+2  : float     atomic masses in a.u. (same ordering as geometry)

The script performs the following major steps
---------------------------------------------
1.  Load input files and verify consistency.
2.  Remove overall translation by shifting both geometries to their
    respective centres of mass (mass-weighted).
3.  Remove overall rotation with a *mass-weighted* Kabsch alignment
    so that R_TS is optimally superimposed onto R_min.
4.  Compute ΔR = R_TS_aligned – R_min_c (internal displacements).
5.  Determine the radius of the smallest sphere that contains either
    geometry; use it to compute the external force for the requested pressure.
6.  Loop over a (φ,θ) grid (nφ × nθ) and evaluate W(φ,θ).
7.  Write summary numbers to *output.txt* and the full W matrix to *W.mtx*.

"""
from pathlib import Path
import numpy as np
from scipy.io import mmwrite


# --------------------------------------------------------------------------
# ------------------------- 1.  Utility functions --------------------------
# --------------------------------------------------------------------------

def read_geometry(filename: str) -> np.ndarray:
    """
    Read a geometry file with format:
        N
        x1 y1 z1
        ...
        xN yN zN
    Returns
    -------
    coords : (N, 3) ndarray of floats (Å)
    """
    with open(filename, "r", encoding="utf-8") as f:
        try:
            N = int(f.readline().strip())
        except ValueError as exc:
            raise ValueError(f"First line of {filename} must be an integer.") from exc
        data = np.loadtxt(f)
    if data.shape != (N, 3):
        raise ValueError(
            f"{filename}: expected {N} lines with 3 columns, found shape {data.shape}"
        )
    return data


def read_mass_and_pressure(filename: str) -> tuple[float, np.ndarray]:
    """
    Read the mass / pressure file with format:
        P_ext_GPa
        N
        m1
        ...
        mN
    Returns
    -------
    P_ext_GPa : float
    masses    : (N,) ndarray
    """
    with open(filename, "r", encoding="utf-8") as f:
        try:
            P_ext_GPa = float(f.readline().strip())
        except ValueError as exc:
            raise ValueError(f"First line of {filename} must be a float.") from exc
        try:
            N = int(f.readline().strip())
        except ValueError as exc:
            raise ValueError(f"Second line of {filename} must be an integer.") from exc
        masses = np.loadtxt(f)
    masses = np.asarray(masses, dtype=float).ravel()
    if masses.size != N:
        raise ValueError(
            f"{filename}: expected {N} masses, but found {masses.size} entries."
        )
    return P_ext_GPa, masses


def kabsch_weighted(A: np.ndarray, B: np.ndarray, w: np.ndarray) -> np.ndarray:
    """
    Mass-weighted Kabsch algorithm.
    Given two coordinate sets A and B (N×3) already translated to their
    respective centres of mass, and a weight vector w (masses, length N),
    return the optimal rotation matrix R such that
        B_rot = B @ R
    minimises RMSD between A and B_rot.
    """
    # Correlation matrix (3 × 3)
    C = A.T @ np.diag(w) @ B
    # Singular-value decomposition
    U, _, Vt = np.linalg.svd(C)
    V = Vt.T
    # Correct for improper rotation if necessary
    d = np.sign(np.linalg.det(V @ U.T))
    D = np.diag([1.0, 1.0, d])
    # Optimal rotation
    R = V @ D @ U.T
    return R


# --------------------------------------------------------------------------
# ------------------------- 2.  Main computation ---------------------------
# --------------------------------------------------------------------------

def main() -> None:
    # ---- 2.1  Read inputs -------------------------------------------------
    geom_min = read_geometry("R_min.txt")
    geom_TS  = read_geometry("R_TS.txt")
    P_ext_GPa, masses = read_mass_and_pressure("mass.txt")

    if geom_min.shape != geom_TS.shape:
        raise RuntimeError("R_min and R_TS contain different numbers of atoms.")
    if geom_min.shape[0] != masses.size:
        raise RuntimeError("Mass vector is inconsistent with geometry files.")

    N = geom_min.shape[0]                   # number of atoms
    masses = masses.astype(float)
    mass_total = masses.sum()

    # ---- 2.2  Remove overall translation (shift to centre of mass) -------
    COM_min = (masses @ geom_min) / mass_total     # 1 × 3
    COM_TS  = (masses @ geom_TS)  / mass_total
    R_min_c = geom_min - COM_min                   # centred geometries
    R_TS_c  = geom_TS  - COM_TS

    # ---- 2.3  Remove overall rotation (weighted Kabsch) ------------------
    R_opt = kabsch_weighted(R_min_c, R_TS_c, masses)
    R_TS_aligned = R_TS_c @ R_opt                # now aligned with R_min_c

    # ---- 2.4  Internal displacement --------------------------------------
    delta_R = R_TS_aligned - R_min_c             # N × 3

    # ---- 3.   Sphere radius enclosing both geometries --------------------
    max_dist_Rmin = np.linalg.norm(geom_min, axis=1).max()
    max_dist_RTS  = np.linalg.norm(geom_TS,  axis=1).max()
    R_sphere = max(max_dist_Rmin, max_dist_RTS)  # Å

    # ---- 4.   External force corresponding to the requested pressure -----
    # Force (nN) according to the original formula:
    #   F_ext = P_ext_GPa * 0.01 * π * R_sphere**2
    F_ext = P_ext_GPa * 0.01 * np.pi * R_sphere**2    # nN
    Fext_i = F_ext / N                                # per-atom force

    # ---- 5.   Loop over (φ,θ) grid and evaluate W ------------------------
    nphi   = 200
    ntheta = 100
    dphi   = 2.0 * np.pi / nphi
    dtheta = np.pi / ntheta
    factor_conv = 1.88973 * 627.5 / 82.387            # conversion to a.u. -> kcal/mol

    W      = np.zeros((nphi, ntheta))
    R_mat  = np.zeros_like(W)
    sum_Ea = 0.0

    # Pre-compute norms and unit vectors of R_min 
    norms_Rmin = np.linalg.norm(geom_min, axis=1)
    # Avoid division by zero for any possible origin-coincident atom
    with np.errstate(divide="ignore", invalid="ignore"):
        dR_min_unit = np.where(
            norms_Rmin[:, None] > 0,
            geom_min / norms_Rmin[:, None],
            0.0                                              # (not expected here)
        )

    for iphi in range(nphi):
        phi = iphi * dphi
        cos_phi, sin_phi = np.cos(phi), np.sin(phi)

        for itheta in range(ntheta):
            theta = itheta * dtheta
            sin_theta, cos_theta = np.sin(theta), np.cos(theta)

            # Outward normal vector n(φ,θ)
            n_vec = np.array([sin_theta * cos_phi,
                              sin_theta * sin_phi,
                              cos_theta])

            # ---- 5.1  Largest distance of projection of R_min_c onto
            #           plane perpendicular to n_vec  -----------------------
            #   For each atom i:
            #       Rproj = R_min_c[i] - (R_min_c[i]·n) * n
            #       r_i   = |Rproj|
            #   Keep the maximum r_i
            proj_along_n = R_min_c @ n_vec             # N
            Rproj        = R_min_c - proj_along_n[:, None] * n_vec
            r_i          = np.linalg.norm(Rproj, axis=1)
            r_max        = r_i.max()

            R_mat[iphi, itheta] = (r_max / R_sphere)**2

            # ---- 5.2  Work term W(φ,θ) -----------------------------------
            dot_dR_n = dR_min_unit @ n_vec             # N
            # Feff (N × 3) : broadcast over atoms
            Feff = (-dot_dR_n[:, None] * n_vec) * Fext_i * (r_max / R_sphere)**2
            dotFeff_deltaR = np.einsum("ij,ij->i", Feff, delta_R)  # length N

            # Subtracts each contribution
            suma = -factor_conv * dotFeff_deltaR.sum()
            W[iphi, itheta] = suma

            # Increment integral of Ea variation over sphere
            sum_Ea += (-factor_conv
                       * (dphi * dtheta * sin_theta)
                       * dotFeff_deltaR.sum())

    sum_Ea /= (4.0 * np.pi)
    min_W, max_W = W.min(), W.max()

    # ---- 6.   Write outputs ---------------------------------------------
    with open("output.txt", "w", encoding="utf-8") as f:
        f.write("# =====  Results generated by compute_W.py  =====\n")
        f.write(f"Number of atoms, N                 : {N:d}\n")
        f.write(f"External pressure, P_ext (GPa)     : {P_ext_GPa:.6f}\n")
        f.write(f"Sphere radius, R_sphere (Å)        : {R_sphere:.6f}\n")
        f.write(f"Total external force, F_ext (nN)   : {F_ext:.6f}\n")
        f.write(f"Per-atom external force, F_i (nN)  : {Fext_i:.6f}\n")
        f.write(f"Average ΔEa over sphere (kcal/mol) : {sum_Ea:.6f}\n")
        f.write(f"Minimum W value (kcal/mol)         : {min_W:.6f}\n")
        f.write(f"Maximum W value (kcal/mol)         : {max_W:.6f}\n")
        f.write("# Matrix W(φ,θ) is stored separately in W.mtx\n")

    # Dense Matrix-Market output (coordinate pairs are implicit)
    mmwrite("W.mtx", W, comment="W(phi,theta) generated by compute_W.py")

    print("Computation Ball-Milling-Bell-model finished.")
    print("  → Numerical summary : output.txt")
    print("  → Matrix-Market file: W.mtx")


# --------------------------------------------------------------------------
# --------------------------- 3.  Script entry -----------------------------
# --------------------------------------------------------------------------

if __name__ == "__main__":
    main()
