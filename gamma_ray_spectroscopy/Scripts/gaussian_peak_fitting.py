"""
gaussian_peak_fitting.py
========================
Gaussian photopeak extraction for NaI(Tl) gamma-ray MCA spectra.

Fits a single Gaussian (with optional local linear background subtraction)
to a user-defined region of interest (ROI) in a 1024-channel MCA spectrum.
Returns peak centroid, FWHM, and their 1-sigma uncertainties for use in
energy calibration and detector characterisation.

Key physics:
  - Photoelectric peaks in NaI(Tl) are well-described by a symmetric Gaussian
    because the dominant broadening mechanism is Poisson fluctuation in
    scintillation photon yield (σ/N ∝ 1/√N).
  - Background subtraction removes the Compton continuum contribution to the
    peak ROI, which otherwise inflates the fitted σ and biases the centroid.

Author : Tadeck Jones
Course : PHY3004W — Gamma-Ray Spectroscopy, UCT (2024)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# Model functions
# ---------------------------------------------------------------------------

def gaussian(channel, amplitude, centroid, sigma):
    """
    Symmetric Gaussian photopeak model.

    Parameters
    ----------
    channel   : array  MCA channel numbers
    amplitude : float  Peak height (counts)
    centroid  : float  Peak centre (channels)
    sigma     : float  Gaussian width parameter (channels); FWHM = 2.355 * sigma
    """
    return amplitude * np.exp(-0.5 * ((channel - centroid) / sigma) ** 2)


# ---------------------------------------------------------------------------
# Core fitting function
# ---------------------------------------------------------------------------

def fit_photopeak(channels, counts, roi_start, roi_end,
                  subtract_background=True, plot=False, label=""):
    """
    Fit a Gaussian to a photopeak in a gamma-ray spectrum.

    Parameters
    ----------
    channels           : ndarray  Channel array (integers, typically 0–1023)
    counts             : ndarray  Counts-per-channel (raw MCA output)
    roi_start, roi_end : int      ROI window boundaries [channels]
    subtract_background: bool     If True, subtract a local linear baseline
                                  fitted to the ROI edges before Gaussian fit.
                                  Recommended for peaks sitting on a Compton
                                  continuum (e.g., 60Co 1332 keV peak).
    plot               : bool     Display the fitted peak with residuals
    label              : str      Source label for plot title

    Returns
    -------
    centroid   : float  Peak centroid [channels]
    fwhm       : float  Full Width at Half Maximum [channels]
    sigma_c    : float  1-sigma uncertainty on centroid [channels]
    resolution : float  Energy resolution R = FWHM / centroid (dimensionless)

    Notes
    -----
    Weights for curve_fit are derived from Poisson statistics:
        sigma_y[i] = sqrt(max(counts[i], 1))
    The floor at 1 prevents division-by-zero in empty channels while keeping
    the weight matrix non-singular. absolute_sigma=True ensures pcov is the
    true covariance matrix, not scaled by chi-squared.
    """
    # --- Extract ROI ---
    mask = (channels >= roi_start) & (channels <= roi_end)
    x    = channels[mask].astype(float)
    y    = counts[mask].astype(float)

    # --- Poisson-derived uncertainties: sigma_N = sqrt(N), floored at 1 ---
    sigma_y = np.sqrt(np.maximum(y, 1))

    # --- Optional: subtract local linear background (Compton continuum) ---
    bg_values = np.zeros_like(y)
    if subtract_background:
        # Fit a straight line through the two edge values of the ROI
        bg_coeffs = np.polyfit([x[0], x[-1]], [y[0], y[-1]], deg=1)
        bg_values = np.polyval(bg_coeffs, x)
        y_fit = y - bg_values
    else:
        y_fit = y

    # --- Initial parameter guesses ---
    peak_idx  = np.argmax(y_fit)
    p0 = [y_fit[peak_idx],  # amplitude
          x[peak_idx],       # centroid estimate
          5.0]               # sigma estimate (channels)

    # --- Weighted least-squares Gaussian fit ---
    popt, pcov = curve_fit(
        gaussian, x, y_fit,
        p0=p0,
        sigma=sigma_y,
        absolute_sigma=True,   # pcov reflects true sigma, not reduced chi-sq
        maxfev=10_000
    )

    centroid   = popt[1]
    fwhm       = 2.3548 * abs(popt[2])  # FWHM = 2*sqrt(2*ln2) * sigma_gauss
    sigma_c    = np.sqrt(pcov[1, 1])    # 1-sigma uncertainty on centroid
    resolution = fwhm / centroid        # dimensionless; multiply by 100 for %

    # --- Optional: diagnostic plot ---
    if plot:
        x_model = np.linspace(x[0], x[-1], 500)
        fig, axes = plt.subplots(2, 1, figsize=(9, 7),
                                 gridspec_kw={"height_ratios": [3, 1]})

        ax = axes[0]
        ax.step(x, y, where="mid", color="black", linewidth=0.9, label="Raw data")
        if subtract_background:
            ax.plot(x, bg_values, ":", color="gray", linewidth=1, label="Linear BG")
        ax.plot(x_model,
                gaussian(x_model, *popt) + np.polyval(bg_coeffs if subtract_background else [0, 0], x_model),
                "--r", linewidth=1.5,
                label=f"Gaussian fit\n"
                      f"μ = {centroid:.2f} ± {sigma_c:.2f} ch\n"
                      f"FWHM = {fwhm:.2f} ch\n"
                      f"R = {resolution*100:.1f}%")
        ax.fill_between(x_model, gaussian(x_model, *popt), alpha=0.2, color="cyan")
        ax.set_ylabel("Counts per Channel")
        ax.set_title(f"Photopeak Gaussian Fit — {label}")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

        # Residuals panel
        ax2 = axes[1]
        y_model_at_x = gaussian(x, *popt) + (bg_values if subtract_background else 0)
        residuals = y - y_model_at_x
        ax2.bar(x, residuals / sigma_y, width=1.0, color="steelblue", alpha=0.7)
        ax2.axhline(0, color="black", linewidth=0.8)
        ax2.set_xlabel("MCA Channel")
        ax2.set_ylabel("Residual (σ)")
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    return centroid, fwhm, sigma_c, resolution


# ---------------------------------------------------------------------------
# Spectrum I/O
# ---------------------------------------------------------------------------

def load_spectrum(filepath):
    """
    Load a two-column MCA spectrum file (channel <whitespace> counts).

    Parameters
    ----------
    filepath : str  Path to .txt file (first row is header, skipped)

    Returns
    -------
    channels : ndarray  int, shape (N,)
    counts   : ndarray  float, shape (N,)
    """
    data     = np.loadtxt(filepath, skiprows=1)
    channels = data[:, 0].astype(int)
    counts   = data[:, 1]
    return channels, counts


# ---------------------------------------------------------------------------
# Quick demo — run as script
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import os

    data_dir = os.path.join(os.path.dirname(__file__), "../data/")

    # Demonstrate on 137Cs — 661.66 keV photopeak near channel 250
    print("=" * 55)
    print("  Demo: ¹³⁷Cs photopeak fit (661.66 keV)")
    print("=" * 55)

    channels, counts = load_spectrum(data_dir + "137Cs.txt")
    centroid, fwhm, sigma_c, R = fit_photopeak(
        channels, counts,
        roi_start=220, roi_end=280,
        subtract_background=True,
        plot=True,
        label="¹³⁷Cs  —  661.66 keV"
    )

    print(f"  Centroid   : {centroid:.2f} ± {sigma_c:.2f} channels")
    print(f"  FWHM       : {fwhm:.2f} channels")
    print(f"  Resolution : {R * 100:.1f}%  (spec ≤ 8.5%)")
