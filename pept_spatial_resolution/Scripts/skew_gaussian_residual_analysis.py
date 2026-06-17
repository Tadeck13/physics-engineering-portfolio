"""
skew_gaussian_residual_analysis.py
====================================
Spatial resolution analysis for PEPT circular motion tracking.

Implements the full residual-to-resolution-metric pipeline:
    1. Load optimised tracking parameters and fit results from CSV
    2. Reconstruct the combined residuals δ(R) from fit parameters (or load raw)
    3. Fit skew Gaussian distribution to each radius's residuals
    4. Extract μ_δ(R) and σ_δ(R) as spatial resolution metrics
    5. Plot μ and σ vs radius A to identify the tracking breakdown threshold

The skew Gaussian is required (vs symmetric Gaussian) because δ(R) is defined as:
    δ(R) = √(d(X)² + d(Y)² + d(Z)²)
which is strictly non-negative. This produces a right-skewed distribution
bounded at zero, with the skew parameter α capturing the asymmetry.

Key finding: for A ≥ 0.50 mm, both μ_δ(R) ≈ 0.21 mm and σ_δ(R) ≈ 0.30 mm
are statistically constant, establishing these as intrinsic system parameters
of the HR++ + Birmingham algorithm combination. Below A = 0.50 mm, both
fit parameter uncertainties increase dramatically, indicating tracking breakdown.

Attribution:
    - GATE simulation framework: R. Perin, UCT / iThemba LABS (Perin et al. 2023)
    - Birmingham algorithm optimisation: R. Perin et al. (2023)
    - Residual analysis implementation: T. Jones (2024)

References:
    [1] Parker et al., NIM A, 326(3):592–607, 1993. (Birmingham algorithm)
    [2] Perin et al., Applied Sciences, 13(11):6690, 2023. (HR++ GATE simulation)
    [3] Leadbeater et al., SAIP 2018. (0.5 mm experimental uncertainty)

Author : Tadeck Jones
Date   : September 2024
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# Skew Gaussian model — the core statistical model for δ(R) distributions
# ---------------------------------------------------------------------------

def skew_gaussian_pdf(x, mu, sigma, alpha, amplitude):
    """
    Skew Gaussian probability density function.

    Constructed as the product of a normal PDF and its CDF (scaled):
        f(x) = 2·A · g(x, µ, σ) · G(α·x, µ, σ)

    where g is the Gaussian PDF and G is the Gaussian CDF. The skew
    parameter α controls the asymmetry of the distribution:
        α > 0  → right-skewed (long right tail) — typical for δ(R)
        α = 0  → reduces to symmetric Gaussian
        α < 0  → left-skewed

    Parameters
    ----------
    x         : array  Values at which to evaluate the PDF
    mu        : float  Location parameter (≈ peak position)
    sigma     : float  Scale parameter (≈ spread)
    alpha     : float  Skewness parameter
    amplitude : float  Overall scaling factor A

    Returns
    -------
    array  Evaluated PDF values [mm⁻¹]
    """
    norm_pdf = (1.0 / (sigma * np.sqrt(2 * np.pi))) * np.exp(
        -0.5 * ((x - mu) / sigma) ** 2
    )
    norm_cdf = 0.5 * (1 + erf((alpha * (x - mu)) / (sigma * np.sqrt(2))))
    return 2.0 * amplitude * norm_pdf * norm_cdf


# ---------------------------------------------------------------------------
# Residual fitting
# ---------------------------------------------------------------------------

def fit_residuals(delta_R, n_bins=30, label="", plot=False):
    """
    Fit a skew Gaussian to a combined PEPT residual distribution.

    Parameters
    ----------
    delta_R : array  Combined 3D residuals δ(R) = √(dX² + dY² + dZ²) [mm]
    n_bins  : int    Number of histogram bins for PMF construction
    label   : str    Label for plot title
    plot    : bool   Show diagnostic plot of fit

    Returns
    -------
    mu        : float  Fitted location parameter [mm]
    sigma     : float  Fitted scale parameter [mm]
    sigma_mu  : float  1-σ uncertainty on mu [mm]
    sigma_sig : float  1-σ uncertainty on sigma [mm]
    alpha     : float  Fitted skewness parameter
    """
    # Construct empirical PMF from histogram
    counts, edges = np.histogram(delta_R, bins=n_bins, density=True)
    centres        = 0.5 * (edges[:-1] + edges[1:])

    # Initial parameter estimates from moments
    mu0    = np.mean(delta_R)
    sigma0 = np.std(delta_R)
    alpha0 = 2.0           # expected right skew for δ(R)
    A0     = counts.max()

    p0 = [mu0, sigma0, alpha0, A0]

    try:
        popt, pcov = curve_fit(
            skew_gaussian_pdf, centres, counts,
            p0=p0, maxfev=10_000,
            bounds=([0, 0, -np.inf, 0], [np.inf, np.inf, np.inf, np.inf])
        )
    except RuntimeError as e:
        print(f"  [WARNING] Fit did not converge for {label}: {e}")
        return mu0, sigma0, np.nan, np.nan, alpha0

    mu, sigma, alpha, amplitude = popt
    sigma_mu  = np.sqrt(pcov[0, 0])
    sigma_sig = np.sqrt(pcov[1, 1])

    if plot:
        x_model = np.linspace(0, max(centres), 300)
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(centres, counts, width=(edges[1]-edges[0]),
               color="white", edgecolor="black", label=r"$\delta(R)$", linewidth=0.8)
        ax.plot(x_model, skew_gaussian_pdf(x_model, *popt),
                color="red", linewidth=1.5, label="Skew Gaussian fit")
        ax.set_xlabel(r"$\delta(R)$ (mm)", fontsize=11)
        ax.set_ylabel(r"PMF (mm$^{-1}$)", fontsize=11)
        ax.set_title(f"Combined Residuals — {label}", fontsize=11)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    return mu, sigma, sigma_mu, sigma_sig, alpha


# ---------------------------------------------------------------------------
# Load optimised parameters from CSV
# ---------------------------------------------------------------------------

def load_tracking_params(csv_path):
    """
    Load radius, tracking parameters, and skew Gaussian fit results from CSV.

    CSV format (from optTrackingParams.csv):
        radius (mm), SS, f_opt, mu (mm), u(mu), sigma (mm), u(sigma), L (Hz)

    Returns
    -------
    dict with arrays: radii, SS, f_opt, mu, u_mu, sigma, u_sigma, L
    """
    data = np.loadtxt(csv_path, delimiter=",", comments="#")
    return {
        "radii"  : data[:, 0],
        "SS"     : data[:, 1].astype(int),
        "f_opt"  : data[:, 2],
        "mu"     : data[:, 3],
        "u_mu"   : data[:, 4],
        "sigma"  : data[:, 5],
        "u_sigma": data[:, 6],
        "L"      : data[:, 7],
    }


# ---------------------------------------------------------------------------
# Spatial resolution summary plot
# ---------------------------------------------------------------------------

def plot_resolution_vs_radius(params, save_path=None):
    """
    Plot μ_δ(R) and σ_δ(R) vs radius A — the primary result figure.

    Replicates Figure 9 from Jones 2024. The key visual signature is that
    μ and σ are statistically constant for A ≥ 0.50 mm, but the uncertainties
    on both parameters increase sharply at A = 0.10 mm.

    Parameters
    ----------
    params    : dict  Output from load_tracking_params()
    save_path : str   Save path for figure (optional)
    """
    radii   = params["radii"]
    mu      = params["mu"]
    u_mu    = params["u_mu"]
    sigma   = params["sigma"]
    u_sigma = params["u_sigma"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    # Panel (a): mean μ_δ(R) vs A
    ax = axes[0]
    ax.errorbar(radii, mu, yerr=u_mu,
                fmt="o", color="black", markersize=6, capsize=4, linewidth=1.2)
    # Plateau line (A ≥ 0.5 mm)
    plateau_mu = np.mean(mu[radii >= 0.5])
    ax.axhline(plateau_mu, color="gray", linestyle="--", linewidth=0.8,
               label=f"Plateau  μ = {plateau_mu:.3f} mm")
    # Mark breakdown region
    ax.axvspan(0, 0.5, alpha=0.10, color="red", label="Breakdown interval")
    ax.set_xlabel("A (mm)", fontsize=12)
    ax.set_ylabel(r"$\mu_{\delta(R)}$ (mm)", fontsize=12)
    ax.set_title(r"(a) Fit Mean $\mu_{\delta(R)}$ vs. Radius $A$", fontsize=11)
    ax.set_xlim(-0.1, 5.3)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel (b): standard deviation σ_δ(R) vs A
    ax2 = axes[1]
    ax2.errorbar(radii, sigma, yerr=u_sigma,
                 fmt="o", color="black", markersize=6, capsize=4, linewidth=1.2)
    plateau_sigma = np.mean(sigma[radii >= 0.5])
    ax2.axhline(plateau_sigma, color="gray", linestyle="--", linewidth=0.8,
                label=f"Plateau  σ = {plateau_sigma:.3f} mm")
    ax2.axvspan(0, 0.5, alpha=0.10, color="red", label="Breakdown interval")
    ax2.set_xlabel("A (mm)", fontsize=12)
    ax2.set_ylabel(r"$\sigma_{\delta(R)}$ (mm)", fontsize=12)
    ax2.set_title(r"(b) Fit Std Dev $\sigma_{\delta(R)}$ vs. Radius $A$", fontsize=11)
    ax2.set_xlim(-0.1, 5.3)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    plt.suptitle(
        "PEPT Spatial Resolution — Birmingham Algorithm + Siemens HR++\n"
        r"Tracking breaks down for $A < 0.5$ mm  "
        r"[$\mu \approx 0.21$ mm, $\sigma \approx 0.30$ mm for $A \geq 0.5$ mm]",
        fontsize=11
    )
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {save_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------

def print_summary_table(params):
    """Print a formatted summary of tracking parameters and resolution metrics."""
    print("\n" + "=" * 80)
    print("  PEPT SPATIAL RESOLUTION SUMMARY — HR++ + Birmingham Algorithm")
    print("=" * 80)
    print(f"\n  {'A (mm)':>8} {'SS':>6} {'f_opt':>7} {'μ (mm)':>10} {'u(μ)':>8} "
          f"{'σ (mm)':>10} {'u(σ)':>8} {'L (Hz)':>10}")
    print("  " + "-" * 75)

    radii   = params["radii"]
    plateau_mu    = np.mean(params["mu"][radii >= 0.5])
    plateau_sigma = np.mean(params["sigma"][radii >= 0.5])

    for i in range(len(radii)):
        flag = "← breakdown" if radii[i] < 0.5 else ""
        print(f"  {radii[i]:>8.2f} {params['SS'][i]:>6d} {params['f_opt'][i]:>7.3f} "
              f"{params['mu'][i]:>10.4f} {params['u_mu'][i]:>8.4f} "
              f"{params['sigma'][i]:>10.4f} {params['u_sigma'][i]:>8.4f} "
              f"{params['L'][i]:>10.1f}  {flag}")

    print()
    print(f"  Plateau (A ≥ 0.5 mm):")
    print(f"    μ_δ(R) = {plateau_mu:.3f} mm   (tracking noise floor)")
    print(f"    σ_δ(R) = {plateau_sigma:.3f} mm   (spread of tracking uncertainty)")
    print()
    print(f"  Tracking breakdown interval: 0.10 mm < A < 0.50 mm")
    print(f"  Consistent with Leadbeater et al. 2018 experimental minimum: ~0.5 mm")
    print()

    # Significance of breakdown
    u_mu_breakdown = params["u_mu"][radii < 0.5][0]
    u_mu_plateau   = np.mean(params["u_mu"][radii >= 0.5])
    print(f"  Uncertainty increase at A=0.10 mm vs plateau:")
    print(f"    u(μ):  {u_mu_breakdown:.4f} mm vs {u_mu_plateau:.4f} mm  "
          f"→ {u_mu_breakdown/u_mu_plateau:.1f}× larger  (tracking becoming unreliable)")
    print("=" * 80)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    data_dir    = os.path.join(os.path.dirname(__file__), "../data/")
    figures_dir = os.path.join(os.path.dirname(__file__), "../figures/")
    os.makedirs(figures_dir, exist_ok=True)

    print("\n=== PEPT Spatial Resolution Analysis — Circular Motion Benchmark ===\n")

    # Load pre-computed fit parameters (from Gaussian process-optimised tracking)
    csv_path = os.path.join(data_dir, "optTrackingParams.csv")
    params   = load_tracking_params(csv_path)

    # Print summary table
    print_summary_table(params)

    # Generate resolution vs radius plot
    print("Generating spatial resolution plot...")
    plot_resolution_vs_radius(
        params,
        save_path=os.path.join(figures_dir, "resolution_vs_radius_reproduced.png")
    )

    # Demonstrate skew Gaussian with a synthetic residual distribution
    # (matching the A=5mm measured parameters for illustration)
    print("\nDemonstrating skew Gaussian fit on synthetic δ(R) distribution...")
    rng        = np.random.default_rng(42)
    # Generate synthetic δ(R) ≈ |Normal(0.207, 0.308)| (rough approximation)
    synthetic_residuals = np.abs(rng.normal(loc=0.0, scale=0.36, size=5000))
    mu_fit, sig_fit, u_mu, u_sig, alpha_fit = fit_residuals(
        synthetic_residuals, n_bins=30, label="Synthetic A=5mm (demo)", plot=True
    )
    print(f"  Demo fit: μ = {mu_fit:.3f} ± {u_mu:.3f} mm  |  σ = {sig_fit:.3f} ± {u_sig:.3f} mm")
    print(f"  Skewness α = {alpha_fit:.2f}  (positive = right-skewed, as expected for |residual|)")
