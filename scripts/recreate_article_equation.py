"""
Recreation of the DMD limitation example from the article using the equation:

f(x,t) = f1(x,t) + f2(x,t)
f1(x,t) = sech(x + 6 - t) exp(i 2.3 t)
f2(x,t) = 2 sech(x) tanh(x) exp(i 2.8 t)

This script computes the data and the DMD results from scratch. No pre-generated
results are used.
"""

from __future__ import annotations

import gc
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["figure.max_open_warning"] = 0
import numpy as np
import pandas as pd


# -----------------------------
# User-editable experiment setup
# -----------------------------
NX = 400
NT = 200
X_MIN = -10.0
X_MAX = 10.0
T_MIN = 0.0
T_MAX = 4.0 * np.pi
RANKS = [2, 5, 10]
OUTPUT_DIR = Path("outputs")


# -----------------------------
# Signal definition
# -----------------------------
def sech(z: np.ndarray) -> np.ndarray:
    """Hyperbolic secant: sech(z) = 1/cosh(z)."""
    return 1.0 / np.cosh(z)


def make_signal(x: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns f, f1, f2 as matrices of shape (space points, time points).

    Rows correspond to spatial locations x.
    Columns correspond to time snapshots t.
    """
    X, T = np.meshgrid(x, t, indexing="ij")

    f1 = sech(X + 6.0 - T) * np.exp(1j * 2.3 * T)
    f2 = 2.0 * sech(X) * np.tanh(X) * np.exp(1j * 2.8 * T)
    f = f1 + f2
    return f, f1, f2


# -----------------------------
# DMD implementation
# -----------------------------
def exact_dmd(X: np.ndarray, t: np.ndarray, rank: int, svd_cache: tuple[np.ndarray, np.ndarray, np.ndarray] | None = None) -> dict:
    """
    Exact DMD implementation.

    X is a snapshot matrix with columns x_k.
    X1 = [x_1, x_2, ..., x_{m-1}]
    X2 = [x_2, x_3, ..., x_m]

    DMD finds a linear operator A such that X2 ≈ A X1.
    Instead of forming huge A directly, it forms a reduced operator A_tilde
    in rank-r SVD coordinates.

    svd_cache lets run_all compute the expensive SVD once and reuse it for
    rank 2, 5, and 10. This makes the script faster and reduces memory use.
    """
    dt = float(t[1] - t[0])
    X1 = X[:, :-1]
    X2 = X[:, 1:]

    if svd_cache is None:
        U, s, Vh = np.linalg.svd(X1, full_matrices=False)
    else:
        U, s, Vh = svd_cache
    r = min(rank, len(s))

    Ur = U[:, :r]
    sr = s[:r]
    Vr = Vh.conj().T[:, :r]

    # Reduced-order DMD operator: A_tilde = U_r^* X2 V_r Sigma_r^{-1}
    A_tilde = Ur.conj().T @ X2 @ Vr @ np.diag(1.0 / sr)

    eigenvalues, W = np.linalg.eig(A_tilde)

    # Exact DMD modes
    Phi = X2 @ Vr @ np.diag(1.0 / sr) @ W

    # Continuous-time DMD eigenvalues: lambda = exp(omega dt)
    omega = np.log(eigenvalues) / dt

    # Initial amplitudes
    b = np.linalg.pinv(Phi) @ X[:, 0]

    # Reconstruct all time snapshots
    time_dynamics = np.zeros((r, len(t)), dtype=complex)
    for k, tk in enumerate(t):
        time_dynamics[:, k] = b * np.exp(omega * tk)
    X_dmd = Phi @ time_dynamics

    rel_error = np.linalg.norm(X - X_dmd, ord="fro") / np.linalg.norm(X, ord="fro")

    total_energy = np.sum(s**2)
    cumulative_energy_rank = np.sum(s[:r] ** 2) / total_energy

    return {
        "rank": r,
        "dt": dt,
        "X1_shape": X1.shape,
        "X2_shape": X2.shape,
        "svd_singular_values": s,
        "A_tilde": A_tilde,
        "eigenvalues": eigenvalues,
        "omega": omega,
        "modes": Phi,
        "amplitudes": b,
        "reconstruction": X_dmd,
        "relative_error": rel_error,
        "cumulative_energy_rank": cumulative_energy_rank,
    }


# -----------------------------
# Plotting helpers
# -----------------------------
def save_complex_field_plots(x: np.ndarray, t: np.ndarray, f: np.ndarray, f1: np.ndarray, f2: np.ndarray) -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    extent = [t[0], t[-1], x[0], x[-1]]

    for name, data in [("f_total", f), ("f1_translating", f1), ("f2_stationary", f2)]:
        fig, ax = plt.subplots(figsize=(8, 5))
        im = ax.imshow(np.asarray(np.real(data), dtype=np.float32), aspect="auto", origin="lower", extent=extent, interpolation="nearest")
        ax.set_title(f"Real part of {name}")
        ax.set_xlabel("t")
        ax.set_ylabel("x")
        fig.colorbar(im, ax=ax, label="Real amplitude")
        fig.tight_layout()
        fig.savefig(OUTPUT_DIR / f"01_{name}_real.png", dpi=120)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 5))
        im = ax.imshow(np.asarray(np.abs(data), dtype=np.float32), aspect="auto", origin="lower", extent=extent, interpolation="nearest")
        ax.set_title(f"Magnitude of {name}")
        ax.set_xlabel("t")
        ax.set_ylabel("x")
        fig.colorbar(im, ax=ax, label="Magnitude")
        fig.tight_layout()
        fig.savefig(OUTPUT_DIR / f"02_{name}_magnitude.png", dpi=120)
        plt.close(fig)
        gc.collect()


def save_svd_plot(s: np.ndarray) -> None:
    energy = s**2 / np.sum(s**2)
    cumulative = np.cumsum(s**2) / np.sum(s**2)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(np.arange(1, len(s) + 1), s, marker="o", markersize=3, linewidth=1)
    ax.set_title("Singular values of X1")
    ax.set_xlabel("Mode number")
    ax.set_ylabel("Singular value")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "03_svd_singular_values.png", dpi=120)
    plt.close(fig)
    gc.collect()

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(np.arange(1, len(cumulative) + 1), cumulative, marker="o", markersize=3, linewidth=1)
    ax.set_title("Cumulative SVD energy")
    ax.set_xlabel("Rank")
    ax.set_ylabel("Cumulative energy")
    ax.set_ylim(0, 1.01)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "04_svd_cumulative_energy.png", dpi=120)
    plt.close(fig)
    gc.collect()

    pd.DataFrame({
        "mode_number": np.arange(1, len(s) + 1),
        "singular_value": s,
        "energy_fraction": energy,
        "cumulative_energy": cumulative,
    }).to_csv(OUTPUT_DIR / "svd_energy_table.csv", index=False)


def save_rank_outputs(x: np.ndarray, t: np.ndarray, X: np.ndarray, result: dict) -> dict:
    r = result["rank"]
    rank_dir = OUTPUT_DIR / f"rank_{r}"
    rank_dir.mkdir(parents=True, exist_ok=True)
    X_dmd = result["reconstruction"]
    eig = result["eigenvalues"]
    omega = result["omega"]
    Phi = result["modes"]

    # Reconstruction comparison at selected times
    selected_indices = [0, len(t) // 3, 2 * len(t) // 3, len(t) - 1]
    fig, ax = plt.subplots(figsize=(8, 5))
    for idx in selected_indices:
        ax.plot(x, np.real(X[:, idx]), linewidth=2, label=f"True t={t[idx]:.2f}")
        ax.plot(x, np.real(X_dmd[:, idx]), linestyle="--", linewidth=1.5, label=f"DMD t={t[idx]:.2f}")
    ax.set_title(f"True vs DMD reconstruction, rank {r}")
    ax.set_xlabel("x")
    ax.set_ylabel("Real amplitude")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(rank_dir / "reconstruction_line_comparison.png", dpi=120)
    plt.close(fig)
    gc.collect()

    # Error field
    extent = [t[0], t[-1], x[0], x[-1]]
    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(np.asarray(np.real(X - X_dmd), dtype=np.float32), aspect="auto", origin="lower", extent=extent, interpolation="nearest")
    ax.set_title(f"Real reconstruction error field, rank {r}")
    ax.set_xlabel("t")
    ax.set_ylabel("x")
    fig.colorbar(im, ax=ax, label="Real error")
    fig.tight_layout()
    fig.savefig(rank_dir / "reconstruction_error_field.png", dpi=120)
    plt.close(fig)
    gc.collect()

    # Eigenvalues in complex plane
    theta = np.linspace(0, 2 * np.pi, 400)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(np.cos(theta), np.sin(theta), linestyle="--", linewidth=1, label="Unit circle")
    ax.scatter(np.real(eig), np.imag(eig), s=45, label="DMD eigenvalues")
    ax.axhline(0, linewidth=0.8)
    ax.axvline(0, linewidth=0.8)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"DMD eigenvalues, rank {r}")
    ax.set_xlabel("Real(lambda)")
    ax.set_ylabel("Imag(lambda)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(rank_dir / "dmd_eigenvalues.png", dpi=120)
    plt.close(fig)
    gc.collect()

    # DMD mode magnitudes
    fig, ax = plt.subplots(figsize=(8, 5))
    for j in range(min(r, 10)):
        mode = Phi[:, j]
        scale = np.max(np.abs(mode))
        if scale > 0:
            mode = mode / scale
        ax.plot(x, np.abs(mode), label=f"Mode {j+1}")
    ax.set_title(f"DMD mode magnitudes, rank {r}")
    ax.set_xlabel("x")
    ax.set_ylabel("Normalized |mode|")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(rank_dir / "dmd_mode_magnitudes.png", dpi=120)
    plt.close(fig)
    gc.collect()

    eig_table = pd.DataFrame({
        "mode": np.arange(1, r + 1),
        "lambda_real": np.real(eig),
        "lambda_imag": np.imag(eig),
        "lambda_abs": np.abs(eig),
        "lambda_angle_rad": np.angle(eig),
        "omega_real_growth_decay": np.real(omega),
        "omega_imag_angular_frequency": np.imag(omega),
        "frequency_cycles_per_time": np.imag(omega) / (2 * np.pi),
        "amplitude_abs": np.abs(result["amplitudes"]),
    })
    eig_table.to_csv(rank_dir / "dmd_eigenvalue_frequency_table.csv", index=False)

    summary = {
        "rank": r,
        "relative_reconstruction_error": float(result["relative_error"]),
        "cumulative_svd_energy_at_rank": float(result["cumulative_energy_rank"]),
        "dt": float(result["dt"]),
        "X1_shape": list(result["X1_shape"]),
        "X2_shape": list(result["X2_shape"]),
    }
    with open(rank_dir / "summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    return summary


def main() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)

    x = np.linspace(X_MIN, X_MAX, NX)
    t = np.linspace(T_MIN, T_MAX, NT)
    X, f1, f2 = make_signal(x, t)

    # Save generated dataset so you can inspect what DMD is receiving.
    np.savez(OUTPUT_DIR / "generated_dataset_exact_equation.npz", x=x, t=t, f=X, f1=f1, f2=f2)

    save_complex_field_plots(x, t, X, f1, f2)

    # Compute this expensive SVD one time and reuse it for all requested ranks.
    # This avoids repeated memory-heavy decompositions on slower laptops.
    X1 = X[:, :-1]
    svd_cache = np.linalg.svd(X1, full_matrices=False)
    save_svd_plot(svd_cache[1])

    all_summaries = []
    for r in RANKS:
        result = exact_dmd(X, t, r, svd_cache=svd_cache)
        all_summaries.append(save_rank_outputs(x, t, X, result))
        gc.collect()

    pd.DataFrame(all_summaries).to_csv(OUTPUT_DIR / "rank_error_summary.csv", index=False)

    print("DMD limitation recreation completed.")
    print(f"Equation: sech(x+6-t)*exp(i*2.3*t) + 2*sech(x)*tanh(x)*exp(i*2.8*t)")
    print(f"Grid: NX={NX}, NT={NT}, x=[{X_MIN}, {X_MAX}], t=[{T_MIN}, 4*pi]")
    print("\nSummary:")
    for item in all_summaries:
        print(
            f"  rank={item['rank']:>2} | error={item['relative_reconstruction_error']:.6e} "
            f"| SVD energy={item['cumulative_svd_energy_at_rank']:.6f}"
        )
    print(f"\nOutputs saved in: {OUTPUT_DIR.resolve()}")
    print("Important files:")
    print("  outputs/rank_error_summary.csv")
    print("  outputs/svd_energy_table.csv")
    print("  outputs/rank_2/dmd_eigenvalue_frequency_table.csv")
    print("  outputs/rank_5/dmd_eigenvalue_frequency_table.csv")
    print("  outputs/rank_10/dmd_eigenvalue_frequency_table.csv")


if __name__ == "__main__":
    main()
