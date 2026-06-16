"""
half_life_extraction.py
=======================
Weighted linear least-squares half-life measurement from linearized decay data.

Extracts the decay constant λ and half-life t½ for each spectroscopic range
of the ²⁷Al neutron activation experiment. The approach uses the linearization:

    ln(N(t)/N₀) = −λt

which converts a non-linear exponential fit into a linear regression problem,
making the result less sensitive to initial parameter guesses and allowing the
use of Poisson-derived weights in log space.

Key engineering decisions:
    1. FITTING WINDOW: Only the early portion of each decay curve (where the
       activation product dominates over background) is used for the fit.
       Late-time points are excluded because the background plateau invalidates
       the single-component exponential model.

    2. WEIGHTS: In log space, the Poisson uncertainty on ln(N) ≈ 1/√N, so
       the optimal weight for each point is w_i = N_i (the raw count). This
       upweights the early high-statistics points that best constrain λ.

    3. FORCED-THROUGH-ORIGIN FIT: The linearized form y = -λt passes through
       the origin by construction (ln(N₀/N₀) = 0 at t = 0). A standard two-
       parameter regression (y = mx + b) is used here for robustness, but
       the intercept should be near zero for a clean measurement.

    4. UNCERTAINTY PROPAGATION: σ_t½ = (ln2 / λ²) · σ_λ, derived from the
       standard propagation of independent errors.

Author : Tadeck Jones
Course : PHY3004W — Neutron Activation Analysis, UCT (2024)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# Reference nuclear data (from NuDat3 / Basunia 2013)
# ---------------------------------------------------------------------------

NUCLEAR_REFERENCE = {
    "28Al": {"t_half_s": 134.7,  "sigma_t": 0.1,  "gamma_keV": [1778.99]},
    "27Mg": {"t_half_s": 567.6,  "sigma_t": 0.6,  "gamma_keV": [843.76, 1014.44]},
}

# Decay window definitions: (filename, label, t_fit_end [s])
RANGE_CATALOGUE = [
    ("Peak1_count_range.txt", "Range 1", 390),
    ("Peak2_count_range.txt", "Range 2", 570),
    ("Peak3_count_range.txt", "Range 3", 450),
    ("Peak4_count_range.txt", "Range 4", 390),
]


# ---------------------------------------------------------------------------
# Core fitting function
# ---------------------------------------------------------------------------

def extract_half_life(time, counts, fit_window_end, label=""):
    """
    Extract decay constant and half-life via linearized weighted regression.

    The decay law N(t) = N₀ exp(-λt) is linearized by taking:
        y_i = ln(N_i / N_0)   →   y_i = -λ · t_i

    Poisson weighting: var(ln N) ≈ 1/N, so weights w_i = N_i.

    Parameters
    ----------
    time          : ndarray  Time stamps [s] (first point is reference t₀)
    counts        : ndarray  Integrated counts at each time step
    fit_window_end: float    Upper time bound for the fit window [s]
    label         : str      Range label for printed output

    Returns
    -------
    lam       : float  Decay constant λ [s⁻¹]
    sigma_lam : float  1-σ uncertainty on λ [s⁻¹]
    t_half    : float  Half-life t½ = ln2/λ [s]
    sigma_t   : float  1-σ uncertainty on t½ [s]
    intercept : float  Fit intercept (should be ≈ 0)
    """
    # --- Select fitting window ---
    mask   = time <= fit_window_end
    t      = time[mask]
    N      = counts[mask].astype(float)

    # --- Linearize: y = ln(N / N[0]) ---
    # Use the first point as the reference N₀ (t = 30 s measurement)
    N0     = N[0]
    y      = np.log(N / N0)

    # --- Poisson weights in log space: w_i = N_i ---
    weights = N.copy()                          # w_i = N_i  (not 1/σ², but N)

    # --- Weighted linear fit: y = slope*t + intercept ---
    def linear(t_, slope, intercept):
        return slope * t_ + intercept

    # Sigma for curve_fit: σ_i = 1/√N_i → absolute_sigma=True
    sigma_y = 1.0 / np.sqrt(np.maximum(N, 1))  # Poisson: σ(ln N) ≈ 1/√N

    popt, pcov = curve_fit(
        linear, t, y,
        p0=[-0.005, 0.0],
        sigma=sigma_y,
        absolute_sigma=True
    )

    slope     = popt[0]          # = -λ
    intercept = popt[1]          # should be ≈ 0

    lam       = -slope
    sigma_lam = np.sqrt(pcov[0, 0])

    t_half    = np.log(2) / lam
    sigma_t   = np.log(2) * sigma_lam / lam**2

    # Goodness-of-fit
    residuals = y - linear(t, *popt)
    chi2      = np.sum((residuals / sigma_y)**2)
    dof       = len(t) - 2
    chi2_red  = chi2 / dof

    if label:
        print(f"  {label}:")
        print(f"    Fit window       : t ≤ {fit_window_end} s  ({len(t)} points)")
        print(f"    Decay constant λ : {lam:.5f} ± {sigma_lam:.5f}  s⁻¹")
        print(f"    Half-life t½     : {t_half:.1f} ± {sigma_t:.1f}  s")
        print(f"    Intercept        : {intercept:.3f}  (ideal: 0.00)")
        print(f"    χ²/dof           : {chi2_red:.2f}")
        print()

    return lam, sigma_lam, t_half, sigma_t, intercept


# ---------------------------------------------------------------------------
# Multi-range analysis
# ---------------------------------------------------------------------------

def run_all_ranges(data_dir=None, figures_dir=None):
    """
    Run half-life extraction for all four spectroscopic ranges and produce
    a 4-panel plot of the linearized fits with residuals.

    Parameters
    ----------
    data_dir    : str or None  Path to data/ folder (auto-detected if None)
    figures_dir : str or None  Path to figures/ folder for saving
    """
    if data_dir is None:
        data_dir = os.path.join(os.path.dirname(__file__), "../data/")
    if figures_dir is None:
        figures_dir = os.path.join(os.path.dirname(__file__), "../figures/")

    os.makedirs(figures_dir, exist_ok=True)

    print("\n" + "=" * 65)
    print("  ²⁷Al NEUTRON ACTIVATION — HALF-LIFE EXTRACTION")
    print("=" * 65)
    print("\n  Reference values from NuDat3:")
    for iso, v in NUCLEAR_REFERENCE.items():
        print(f"    {iso}: t½ = {v['t_half_s']} ± {v['sigma_t']} s,  "
              f"γ at {v['gamma_keV']} keV")
    print()

    colors    = ["steelblue", "darkorange", "forestgreen", "crimson"]
    results   = []

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    axes_flat = axes.flatten()

    for idx, (fname, label, t_end) in enumerate(RANGE_CATALOGUE):
        fpath = os.path.join(data_dir, fname)
        if not os.path.exists(fpath):
            print(f"  [WARNING] File not found: {fpath}")
            continue

        data    = np.loadtxt(fpath, skiprows=1, usecols=(0, 1))
        t, N    = data[:, 0], data[:, 1]

        lam, sig_lam, t_half, sig_t, intercept = extract_half_life(
            t, N, t_end, label=label
        )
        results.append((label, lam, sig_lam, t_half, sig_t))

        # --- Plot linearized data and fit ---
        ax  = axes_flat[idx]
        clr = colors[idx]

        mask = t <= t_end
        t_fit = t[mask]
        N_fit = N[mask].astype(float)
        y_fit = np.log(N_fit / N_fit[0])
        sig_y = 1.0 / np.sqrt(np.maximum(N_fit, 1))

        t_model = np.linspace(t_fit[0], t_fit[-1], 200)
        y_model = -lam * (t_model - t_model[0])

        ax.errorbar(t_fit, y_fit, yerr=sig_y,
                    fmt="+", color=clr, markersize=6, capsize=3,
                    label="Measured ln(N/N₀)", linewidth=0.8)
        ax.plot(t_model, -lam * (t_model - t_model[0]) + intercept, "--k",
                linewidth=1.5,
                label=f"WLS fit:  λ = {lam:.4f} s⁻¹\n"
                      f"         t½ = {t_half:.0f} ± {sig_t:.0f} s")
        ax.set_xlabel("Time [s]", fontsize=10)
        ax.set_ylabel("ln(N / N₀)", fontsize=10)
        ax.set_title(label, fontsize=11, fontweight="bold")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.suptitle(
        "Linearized Exponential Decay Fits — ²⁷Al Neutron Activation\n"
        "[ln(N/N₀) = −λt; Poisson-weighted; truncated at background onset]",
        fontsize=11
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    save_path = os.path.join(figures_dir, "linearized_fits_panel.png")
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"  Saved: {save_path}")
    plt.show()

    # --- Summary table ---
    print("\n" + "=" * 65)
    print("  HALF-LIFE SUMMARY TABLE")
    print("=" * 65)
    print(f"  {'Range':<10} {'λ (s⁻¹)':>14} {'σ_λ':>10} {'t½ (s)':>10} {'σ_t½':>8}")
    print("  " + "-" * 55)
    for label, lam, sig_lam, t_half, sig_t in results:
        print(f"  {label:<10} {lam:>14.5f} {sig_lam:>10.5f} {t_half:>10.1f} {sig_t:>8.1f}")

    print()
    print("  Comparison with nuclear database (NuDat3):")
    print(f"    ²⁸Al: t½ = 134.7 s  →  Range 4 = {results[3][3]:.0f} ± {results[3][4]:.0f} s  "
          f"({abs(results[3][3]-134.7)/results[3][4]:.1f}σ discrepancy)")
    print(f"    ²⁷Mg: t½ = 567.6 s  →  Range 1 = {results[0][3]:.0f} ± {results[0][4]:.0f} s  "
          f"({abs(results[0][3]-567.6)/results[0][4]:.1f}σ discrepancy)")
    print()
    print("  → Systematic offset attributed to manual transfer delay (~60–120 s).")
    print("    This lost the highest-activity early portion of the ²⁸Al decay,")
    print("    causing all measured t½ values to be systematically larger than true.")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_all_ranges()
