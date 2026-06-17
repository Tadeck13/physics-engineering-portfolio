"""
spectral_analysis_pipeline.py
==============================
End-to-end gamma-ray spectroscopy analysis pipeline for NaI(Tl) detector.

Orchestrates the complete workflow from raw MCA spectra to calibrated
energies and isotope identification, generating all figures for the project
portfolio. Designed to be reproducible: running this script re-creates every
quantitative result from the raw data files alone.

Workflow:
    1. Load all calibration source spectra
    2. Fit Gaussian photopeaks with Poisson-weighted uncertainty propagation
    3. Build weighted linear energy calibration (slope + intercept + covariance)
    4. Verify secondary spectral features (Compton edges, backscatter peaks)
    5. Apply calibration to identify unknown source with full sigma_E
    6. Plot energy resolution vs. gamma-ray energy (R ∝ 1/√E verification)
    7. Export summary table to console

Author : Tadeck Jones
Course : PHY3004W — Gamma-Ray Spectroscopy, UCT (2024)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add code/ directory to path for sibling imports
sys.path.insert(0, os.path.dirname(__file__))
from gaussian_peak_fitting import load_spectrum, fit_photopeak
from energy_calibration import calibrate, channel_to_energy


# ---------------------------------------------------------------------------
# Source catalogue — ROI windows tuned for THIS gain configuration
# (137Cs photopeak set to channel ~250 during amplifier calibration)
# ---------------------------------------------------------------------------

SOURCE_CATALOGUE = [
    # (data filename,          known energy [keV],  ROI start, ROI end, label)
    ("137Cs.txt",              661.66,              220, 280, "¹³⁷Cs"),
    ("Co60.txt",              1332.49,              370, 420, "⁶⁰Co (1332 keV)"),
    ("Co60.txt",              1173.23,              340, 390, "⁶⁰Co (1173 keV)"),
    ("22Na.txt",              1274.54,              185, 215, "²²Na (1274 keV)"),
    ("22Na.txt",               511.00,              186, 216, "²²Na (511 keV annihil.)"),
    ("54Mn.txt",               834.85,              285, 325, "⁵⁴Mn"),
    ("133Ba.txt",              356.01,               25,  45, "¹³³Ba"),
    ("57Co.txt",               122.06,               43,  65, "⁵⁷Co"),
    ("unknown_source.txt",       None,              220, 280, "Unknown Source"),
]

# Compton edge and backscatter energies for verification
# Compton edge: E_C  = E_γ * (2α / (1 + 2α))  where α = E_γ / 511 keV
# Backscatter:  E_bs = E_γ / (1 + 2α)
COMPTON_SOURCES = [
    ("¹³⁷Cs",  661.66),
    ("⁶⁰Co",  1332.49),
    ("²²Na",   511.00),
]


# ---------------------------------------------------------------------------
# Helper: compute Compton edge and backscatter energies
# ---------------------------------------------------------------------------

def compton_features(e_gamma_kev):
    """
    Calculate secondary spectral feature energies from first principles.

    Parameters
    ----------
    e_gamma_kev : float  Incident gamma-ray energy [keV]

    Returns
    -------
    e_compton   : float  Compton edge energy [keV]  (max recoil electron KE)
    e_backscatter : float  Backscatter peak energy [keV]  (180° scattered photon)
    """
    m_e_c2 = 511.0                           # electron rest mass energy [keV]
    alpha   = e_gamma_kev / m_e_c2           # reduced photon energy (dimensionless)
    e_compton     = e_gamma_kev * (2 * alpha / (1 + 2 * alpha))
    e_backscatter = e_gamma_kev / (1 + 2 * alpha)
    return e_compton, e_backscatter


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_pipeline(data_dir=None, figures_dir=None):
    """
    Execute the complete spectroscopic analysis pipeline.

    Parameters
    ----------
    data_dir    : str or None  Path to data/ folder (auto-detected if None)
    figures_dir : str or None  Path to figures/ folder for saving plots
    """
    if data_dir is None:
        data_dir = os.path.join(os.path.dirname(__file__), "../data/")
    if figures_dir is None:
        figures_dir = os.path.join(os.path.dirname(__file__), "../figures/")

    print("\n" + "=" * 70)
    print("  NaI(Tl) GAMMA-RAY SPECTROSCOPY — FULL ANALYSIS PIPELINE")
    print("=" * 70)

    # -----------------------------------------------------------------------
    # Step 1: Gaussian photopeak fitting for all sources
    # -----------------------------------------------------------------------
    print("\n[1/5]  Fitting photopeaks...")
    print()
    hdr = f"{'Source':<28} {'E_known':>10} {'Centroid':>10} {'σ_c':>6} {'FWHM':>8} {'R [%]':>7}"
    print(hdr)
    print("-" * len(hdr))

    fit_results = []
    calibration_points = {}

    for fname, e_known, roi0, roi1, label in SOURCE_CATALOGUE:
        fpath = os.path.join(data_dir, fname)
        if not os.path.exists(fpath):
            print(f"  {label:<28} FILE NOT FOUND: {fpath}")
            continue
        try:
            ch, cts = load_spectrum(fpath)
            centroid, fwhm, sigma_c, R = fit_photopeak(
                ch, cts, roi0, roi1,
                subtract_background=True, plot=False
            )
            result = {
                "label": label, "e_known": e_known,
                "centroid": centroid, "fwhm": fwhm,
                "sigma_c": sigma_c, "R": R
            }
            fit_results.append(result)

            # Collect unique calibration points (known energy only)
            if e_known is not None and label not in calibration_points:
                calibration_points[label] = {
                    "energy": e_known, "centroid": centroid, "sigma_c": sigma_c
                }

            e_str = f"{e_known:>10.2f}" if e_known else "   (unknown)"
            print(f"  {label:<28} {e_str} {centroid:>10.2f} {sigma_c:>6.2f} "
                  f"{fwhm:>8.2f} {R*100:>7.1f}")
        except Exception as err:
            print(f"  {label:<28} ERROR: {err}")

    # -----------------------------------------------------------------------
    # Step 2: Energy calibration
    # -----------------------------------------------------------------------
    print("\n[2/5]  Building energy calibration from Gaussian centroids...")
    slope, intercept, pcov = calibrate(sources=calibration_points, plot=False)

    # -----------------------------------------------------------------------
    # Step 3: Verify Compton features
    # -----------------------------------------------------------------------
    print("\n[3/5]  Verifying secondary spectral features (Compton edges):")
    print()
    print(f"  {'Isotope':<10} {'Eγ (keV)':>10} {'Compton Edge':>14} {'Ch (pred)':>12} {'Backscatter':>14} {'Ch (pred)':>12}")
    print("  " + "-" * 76)
    for name, e_g in COMPTON_SOURCES:
        ec, ebs = compton_features(e_g)
        ch_ec  = slope * ec  + intercept
        ch_ebs = slope * ebs + intercept
        print(f"  {name:<10} {e_g:>10.2f} {ec:>14.1f} {ch_ec:>12.1f} "
              f"{ebs:>14.1f} {ch_ebs:>12.1f}")
    print("  (Compare predicted channels to features visible in raw spectra)")

    # -----------------------------------------------------------------------
    # Step 4: Identify unknown source
    # -----------------------------------------------------------------------
    print("\n[4/5]  Identifying unknown source...")
    unknown = next((r for r in fit_results if r["e_known"] is None), None)
    if unknown:
        E_unk, sigma_E = channel_to_energy(
            unknown["centroid"], slope, intercept, pcov,
            sigma_channel=unknown["sigma_c"]
        )
        print(f"\n  Unknown centroid  : {unknown['centroid']:.2f} ± {unknown['sigma_c']:.2f} channels")
        print(f"  Calibrated energy : {E_unk:.2f} ± {sigma_E:.2f} keV")

        # Compare against isotope library
        library = {
            "¹³⁷Cs"        : 661.66,
            "⁵⁴Mn"         : 834.85,
            "⁶⁰Co (1173)"  : 1173.23,
            "²²Na (511)"   : 511.00,
            "⁶⁰Co (1332)"  : 1332.49,
            "²²⁶Ra (609)"  : 609.31,
            "²²⁶Ra (768)"  : 768.36,
        }
        print()
        print(f"  {'Isotope':<20} {'E_lib (keV)':>12} {'ΔE (keV)':>10} {'|Δ/σ|':>8}  {'Match?':>8}")
        print("  " + "-" * 65)
        for iso, e_lib in sorted(library.items(), key=lambda kv: abs(kv[1] - E_unk)):
            delta   = abs(e_lib - E_unk)
            n_sigma = delta / sigma_E
            match   = "✓ MATCH" if n_sigma < 3.0 else ""
            print(f"  {iso:<20} {e_lib:>12.2f} {delta:>10.2f} {n_sigma:>8.1f}  {match}")
    else:
        print("  No unknown source file found in catalogue.")

    # -----------------------------------------------------------------------
    # Step 5: Energy resolution plot (R ∝ 1/√E verification)
    # -----------------------------------------------------------------------
    print("\n[5/5]  Generating energy resolution plot...")
    cal_results = [r for r in fit_results if r["e_known"] is not None]
    energies    = np.array([r["e_known"] for r in cal_results])
    resolutions = np.array([r["R"] * 100 for r in cal_results])   # as percent

    sort_idx = np.argsort(energies)
    E_sorted = energies[sort_idx]
    R_sorted = resolutions[sort_idx]

    # Fit R = k / sqrt(E)
    k_fit = np.mean(R_sorted * np.sqrt(E_sorted))
    E_model = np.linspace(80, 1400, 400)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.scatter(E_sorted, R_sorted, color="black", zorder=5,
               s=60, label="Measured resolution")
    ax.plot(E_model, k_fit / np.sqrt(E_model), "--b", linewidth=1.5,
            label=f"1/√E fit  (k = {k_fit:.1f}  keV^{{0.5}})")
    ax.set_xlabel("γ-Ray Energy [keV]", fontsize=12)
    ax.set_ylabel("Energy Resolution R [%]", fontsize=12)
    ax.set_title("NaI(Tl) Energy Resolution vs. γ-Ray Energy\n"
                 "(R ∝ 1/√E — Poisson statistical origin)", fontsize=12)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    out_path = os.path.join(figures_dir, "resolution_vs_energy.png")
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"  Saved: {out_path}")
    plt.show()

    print("\n[✓]  Pipeline complete.\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_pipeline()
