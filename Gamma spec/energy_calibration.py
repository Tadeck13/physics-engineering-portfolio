"""
energy_calibration.py
=====================
Weighted linear least-squares energy calibration for a NaI(Tl) detector.

Establishes the linear relationship between MCA channel number and gamma-ray
energy using a suite of known calibration isotopes. Propagates uncertainties
on peak centroid positions through to the final calibrated energy output,
including off-diagonal covariance terms.

Calibration model (forward):   C = m * Egamma + b
Identification (inverse):      Egamma = (C - b) / m

Key engineering note:
    The non-zero intercept b = +23.21 channels is a known systematic arising
    from the ORTEC 572A amplifier's Baseline Restorer (BLR) dc-level shift
    and the MCA lower-level discriminator threshold. It is NOT a hardware
    fault — it is absorbed analytically by the calibration fit and corresponds
    to a ~72 keV equivalent zero-energy offset. Ignoring it would introduce a
    systematic error of (b/m) = 23.21/0.32 ≈ 72.5 keV on every measurement.

Author : Tadeck Jones
Course : PHY3004W — Gamma-Ray Spectroscopy, UCT (2024)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# Calibration dataset (centroids from Gaussian fits in gaussian_peak_fitting.py)
# ---------------------------------------------------------------------------

CALIBRATION_SOURCES = {
    # Key: source label
    # Val: known energy [keV], fitted centroid [ch], centroid 1-sigma [ch]
    "137Cs"  : {"energy": 661.66,  "centroid": 236.00, "sigma_c": 0.17},
    "60Co_a" : {"energy": 1173.23, "centroid": 377.01, "sigma_c": 0.42},
    "60Co_b" : {"energy": 1332.49, "centroid": 395.03, "sigma_c": 0.42},
    "22Na"   : {"energy": 1274.54, "centroid": 426.00, "sigma_c": 0.50},
    "133Ba"  : {"energy": 356.01,  "centroid": 156.02, "sigma_c": 0.30},
    "54Mn"   : {"energy": 834.85,  "centroid": 292.00, "sigma_c": 0.35},
    "57Co"   : {"energy": 122.06,  "centroid":  83.00, "sigma_c": 0.25},
}


# ---------------------------------------------------------------------------
# Linear calibration model
# ---------------------------------------------------------------------------

def linear_model(energy, slope, intercept):
    """Forward calibration model:  Channel = slope * Energy[keV] + intercept."""
    return slope * energy + intercept


# ---------------------------------------------------------------------------
# Calibration fit
# ---------------------------------------------------------------------------

def calibrate(sources=None, plot=True):
    """
    Perform weighted linear least-squares energy calibration.

    Weights are derived from centroid uncertainties: w_i = 1 / sigma_c_i**2.
    The full 2x2 covariance matrix is returned for downstream uncertainty
    propagation in channel_to_energy().

    Parameters
    ----------
    sources : dict or None
        Calibration source dictionary. Uses CALIBRATION_SOURCES if None.
    plot : bool
        If True, display calibration curve + residuals panel.

    Returns
    -------
    slope     : float    m  [channels / keV]
    intercept : float    b  [channels]
    cov       : ndarray  2x2 covariance matrix [[s_m^2, s_mb], [s_mb, s_b^2]]
    """
    if sources is None:
        sources = CALIBRATION_SOURCES

    energies  = np.array([v["energy"]   for v in sources.values()])
    centroids = np.array([v["centroid"] for v in sources.values()])
    sigma_c   = np.array([v["sigma_c"]  for v in sources.values()])
    labels    = list(sources.keys())

    # Weighted curve_fit — absolute_sigma=True keeps covariance in original units
    popt, pcov = curve_fit(
        linear_model, energies, centroids,
        p0=[0.32, 23.0],
        sigma=sigma_c,
        absolute_sigma=True
    )
    slope, intercept = popt
    sigma_m = np.sqrt(pcov[0, 0])
    sigma_b = np.sqrt(pcov[1, 1])

    # --- Goodness-of-fit: reduced chi-squared ---
    residuals = centroids - linear_model(energies, *popt)
    chi2      = np.sum((residuals / sigma_c) ** 2)
    dof       = len(energies) - 2          # N data points - 2 fit parameters
    chi2_red  = chi2 / dof

    # --- Console summary ---
    print("\n" + "=" * 55)
    print("  NaI(Tl) ENERGY CALIBRATION — RESULTS")
    print("=" * 55)
    print(f"  Slope (m)        : {slope:.5f} ± {sigma_m:.5f}  ch/keV")
    print(f"  Intercept (b)    : {intercept:.3f} ± {sigma_b:.3f}  ch")
    print(f"  Calibration eqn  : C = {slope:.4f}·Eγ + {intercept:.2f}")
    print(f"  Zero-energy equiv : b/m = {intercept/slope:.1f} keV  (BLR + MCA offset)")
    print(f"  Chi-squared/dof  : {chi2_red:.3f}  (ideal ≈ 1.0)")
    print("=" * 55)

    if plot:
        e_model = np.linspace(0, 1450, 400)
        fig, axes = plt.subplots(
            2, 1, figsize=(9, 8), gridspec_kw={"height_ratios": [3, 1]}
        )

        # --- Calibration curve (top panel) ---
        ax = axes[0]
        ax.errorbar(
            energies, centroids, yerr=sigma_c,
            fmt="+k", markersize=9, capsize=4, label="Calibration sources"
        )
        ax.plot(
            e_model, linear_model(e_model, *popt), "--g",
            label=f"Weighted fit: C = {slope:.4f}·Eγ + {intercept:.2f}"
        )
        # Annotate each point with its isotope label
        for e, c, lbl in zip(energies, centroids, labels):
            ax.annotate(lbl, xy=(e, c), xytext=(8, -12),
                        textcoords="offset points", fontsize=8, color="navy")
        ax.set_xlabel("Known Energy [keV]")
        ax.set_ylabel("MCA Channel")
        ax.set_title("NaI(Tl) Energy Calibration Curve — Weighted Linear Fit")
        ax.legend()
        ax.grid(True, alpha=0.3)

        # --- Residuals (bottom panel) ---
        ax2 = axes[1]
        ax2.axhline(0, color="gray", linewidth=0.8, linestyle="--")
        ax2.errorbar(
            energies, residuals, yerr=sigma_c,
            fmt="ok", markersize=6, capsize=4
        )
        ax2.set_xlabel("Known Energy [keV]")
        ax2.set_ylabel("Residual [ch]")
        ax2.set_title(f"Calibration Residuals  (χ²/dof = {chi2_red:.3f})")
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    return slope, intercept, pcov


# ---------------------------------------------------------------------------
# Uncertainty-propagated energy conversion
# ---------------------------------------------------------------------------

def channel_to_energy(channel, slope, intercept, pcov, sigma_channel=0.0):
    """
    Convert a measured MCA channel to energy with full uncertainty propagation.

    The inverse calibration is:
        Egamma = (C - b) / m

    Uncertainty is propagated using the first-order Taylor expansion,
    including the off-diagonal covariance term Cov(m, b):

        sigma_E^2 = (dE/dC * sigma_C)^2
                  + (dE/db * sigma_b)^2
                  + (dE/dm * sigma_m)^2
                  + 2 * (dE/dm)(dE/db) * Cov(m, b)

    Parameters
    ----------
    channel       : float   Measured peak centroid [channels]
    slope         : float   Calibration slope m
    intercept     : float   Calibration intercept b
    pcov          : ndarray 2x2 covariance matrix from calibrate()
    sigma_channel : float   1-sigma uncertainty on centroid [channels]

    Returns
    -------
    energy  : float  Calibrated energy [keV]
    sigma_E : float  1-sigma uncertainty on energy [keV]
    """
    energy  = (channel - intercept) / slope

    sigma_m = np.sqrt(pcov[0, 0])
    sigma_b = np.sqrt(pcov[1, 1])
    cov_mb  = pcov[0, 1]

    # Partial derivatives of E = (C - b) / m
    dE_dC = 1.0 / slope
    dE_db = -1.0 / slope
    dE_dm = -(channel - intercept) / slope**2

    sigma_E = np.sqrt(
        (dE_dC * sigma_channel) ** 2 +
        (dE_db * sigma_b)       ** 2 +
        (dE_dm * sigma_m)       ** 2 +
        2 * dE_dm * dE_db * cov_mb    # off-diagonal covariance term
    )

    return energy, sigma_E


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Build calibration from the 7-source dataset
    slope, intercept, pcov = calibrate(plot=True)

    # Apply to unknown source (centroid from Gaussian fit)
    unknown_centroid = 251.18
    unknown_sigma_c  = 0.17

    E_unknown, sigma_E = channel_to_energy(
        unknown_centroid, slope, intercept, pcov,
        sigma_channel=unknown_sigma_c
    )

    print(f"\n  Unknown source photopeak: channel {unknown_centroid:.2f} ± {unknown_sigma_c:.2f}")
    print(f"  Calibrated energy       : {E_unknown:.2f} ± {sigma_E:.2f} keV")
    print(f"\n  Note: Nearest known isotope is ¹³⁷Cs at 661.66 keV")
    print(f"  Separation: {abs(E_unknown - 661.66):.1f} keV "
          f"({abs(E_unknown - 661.66)/sigma_E:.0f}σ — statistically distinct)")
