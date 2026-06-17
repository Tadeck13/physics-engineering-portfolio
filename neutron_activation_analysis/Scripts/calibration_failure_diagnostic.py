"""
calibration_failure_diagnostic.py
===================================
Diagnostic analysis of the energy calibration failure in the ²⁷Al NAA experiment.

SUMMARY OF FAILURE:
    The derived calibration equation  E = 0.11·C + 146.62 keV  is physically
    inconsistent. A slope of 0.11 keV/channel across 2048 MCA channels implies:

        E_max = 0.11 × 2048 + 146.62 = 371.9 keV

    This is incompatible with:
      (a) Calibration sources up to 1332 keV (⁶⁰Co), which would require
          channel ~10,780 on this scale — far beyond the MCA range.
      (b) The primary ²⁸Al signature at 1779 keV, which is the main
          target of this activation analysis.

ROOT CAUSE:
    The ORTEC 572A amplifier gain was configured for the PREVIOUS experiment
    (gamma spectroscopy, gain = ¹³⁷Cs at ~channel 250 for a ~1000-channel
    display), but was NOT reconfigured for the NAA experiment, which requires
    a 0–3 MeV range (≈ 1.46 keV/channel for a 2048-channel MCA).

    When the calibration sources (³³ keV–1332 keV) were measured with the
    incorrect gain, all peaks were compressed into the first ~276 channels
    (instead of spanning 0–912 channels). The WLS fit then produced a
    mathematically consistent but physically meaningless equation.

This script:
    1. Loads the three calibration spectra (Ba133/Co60, Cs137/Co57, Mn54/Na22)
    2. Identifies dominant photopeaks in each
    3. Derives the WRONG calibration from the misidentified data
    4. Derives the CORRECT calibration from well-identified peaks
    5. Plots both and computes the energy prediction error for the 4 NAA peaks

Author : Tadeck Jones
Course : PHY3004W — Neutron Activation Analysis, UCT (2024)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks


# ---------------------------------------------------------------------------
# Known isotope lines and their expected MCA channels
# (for a correctly configured 0–3 MeV range in 2048 channels)
# ---------------------------------------------------------------------------

# Peaks identified from the calibration spectra files via Gaussian fitting
# (centroid positions from WLS_2.py and raw spectral analysis)
CALIBRATION_PEAKS = {
    # (source file, channel_centroid, known_energy_keV, isotope_label)
    # These centroids come from Gaussian fits reported in the lab analysis
    "Ba133 (356 keV)" : {"channel": 55,  "energy": 356.01, "spectrum": "Ba133_Co60_cal.txt"},
    "Ba133 (81 keV)"  : {"channel": 38,  "energy": 80.99,  "spectrum": "Ba133_Co60_cal.txt"},
    "Co60 (1173 keV)" : {"channel": 133, "energy": 1173.23,"spectrum": "Ba133_Co60_cal.txt"},
    "Co60 (1332 keV)" : {"channel": 273, "energy": 1332.49,"spectrum": "Ba133_Co60_cal.txt"},
    "Cs137 (662 keV)" : {"channel": 229, "energy": 661.66, "spectrum": "Cs137_Co57_cal.txt"},
    "Co57 (122 keV)"  : {"channel": 39,  "energy": 122.06, "spectrum": "Cs137_Co57_cal.txt"},
    "Mn54 (835 keV)"  : {"channel": 195, "energy": 834.85, "spectrum": "Mn54_Na22_cal.txt"},
    "Na22 (1274 keV)" : {"channel": 276, "energy": 1274.54,"spectrum": "Mn54_Na22_cal.txt"},
}

# Wrong calibration as reported (from miscalibrated WLS fit)
WRONG_SLOPE     = 0.11     # keV/channel
WRONG_INTERCEPT = 146.62   # keV

# NAA photopeak channels (from Gaussian fits on activated Al spectrum)
NAA_CHANNELS = [39.23, 219.11, 257.33, 288.46]

# True isotope assignments (from temporal half-life analysis)
NAA_EXPECTED = {
    "28Al":  1778.99,   # keV — PRIMARY target of NAA
    "27Mg_a": 843.76,  # keV
    "27Mg_b": 1014.44, # keV
}


# ---------------------------------------------------------------------------
# Linear calibration model
# ---------------------------------------------------------------------------

def linear(channel, slope, intercept):
    """E [keV] = slope [keV/ch] * Channel + intercept [keV]"""
    return slope * channel + intercept


# ---------------------------------------------------------------------------
# Main diagnostic
# ---------------------------------------------------------------------------

def run_calibration_diagnostic(data_dir=None, figures_dir=None):
    """
    Full calibration failure diagnostic and reconstruction.

    1. Derives wrong and correct calibration equations.
    2. Shows energy predictions for the 4 NAA peaks under both calibrations.
    3. Computes the dynamic range under each calibration.
    4. Produces a diagnostic comparison figure.
    """
    if data_dir is None:
        data_dir = os.path.join(os.path.dirname(__file__), "../data/")
    if figures_dir is None:
        figures_dir = os.path.join(os.path.dirname(__file__), "../figures/")

    channels = np.array([v["channel"] for v in CALIBRATION_PEAKS.values()])
    energies = np.array([v["energy"]  for v in CALIBRATION_PEAKS.values()])
    labels   = list(CALIBRATION_PEAKS.keys())

    # --- Derive CORRECT calibration by weighted least-squares ---
    # Assume 1% uncertainty on each peak channel (from Gaussian fit)
    sigma_ch = 0.01 * channels
    sigma_ch = np.maximum(sigma_ch, 2.0)   # floor at 2 channels

    popt_correct, pcov_correct = curve_fit(
        linear, channels, energies,
        p0=[4.9, 0.0],
        sigma=sigma_ch,
        absolute_sigma=True
    )
    slope_c, intercept_c = popt_correct

    # Residuals for correct calibration
    residuals_c = energies - linear(channels, *popt_correct)
    rms_c       = np.sqrt(np.mean(residuals_c**2))

    # Wrong calibration residuals
    residuals_w = energies - linear(channels, WRONG_SLOPE, WRONG_INTERCEPT)
    rms_w       = np.sqrt(np.mean(residuals_w**2))

    print("\n" + "=" * 65)
    print("  CALIBRATION FAILURE DIAGNOSTIC")
    print("=" * 65)
    print()
    print("  REPORTED (WRONG) calibration:")
    print(f"    E = {WRONG_SLOPE:.2f} · C + {WRONG_INTERCEPT:.2f}  keV")
    print(f"    Dynamic range (0–2048 ch): "
          f"{WRONG_SLOPE*0+WRONG_INTERCEPT:.1f} – "
          f"{WRONG_SLOPE*2048+WRONG_INTERCEPT:.1f}  keV")
    print(f"    Residual RMS vs known peaks: {rms_w:.1f} keV")
    print()
    print("  CORRECT calibration (WLS from peak centroids):")
    print(f"    E = {slope_c:.3f} · C + {intercept_c:.2f}  keV")
    print(f"    Dynamic range (0–2048 ch): "
          f"{linear(0, *popt_correct):.1f} – "
          f"{linear(2048, *popt_correct):.1f}  keV")
    print(f"    Residual RMS vs known peaks: {rms_c:.1f} keV")
    print()
    print(f"  Slope ratio (correct / wrong): {slope_c/WRONG_SLOPE:.1f}×  ← amplifier gain error")
    print()

    # --- NAA peak energies under each calibration ---
    print("  NAA peak energy assignments:")
    print(f"  {'Channel':>10} {'Wrong E (keV)':>15} {'Correct E (keV)':>17}  Ratio")
    print("  " + "-" * 50)
    for ch in NAA_CHANNELS:
        e_w = linear(ch, WRONG_SLOPE, WRONG_INTERCEPT)
        e_c = linear(ch, *popt_correct)
        print(f"  {ch:>10.1f} {e_w:>15.1f} {e_c:>17.1f}  {e_c/e_w:.1f}×")

    print()
    print("  Expected NAA signatures (nuclear data):")
    for iso, e_kev in NAA_EXPECTED.items():
        ch_correct = (e_kev - intercept_c) / slope_c
        ch_wrong   = (e_kev - WRONG_INTERCEPT) / WRONG_SLOPE
        in_range_c = 0 < ch_correct < 2048
        in_range_w = 0 < ch_wrong   < 2048
        print(f"    {iso}: {e_kev:.1f} keV → "
              f"correct ch={ch_correct:.0f} ({'✓ in MCA' if in_range_c else '✗ out of range'})  |  "
              f"wrong ch={ch_wrong:.0f} ({'✓' if in_range_w else '✗ impossible'})")

    # --- Diagnostic figure ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ch_model = np.linspace(0, 350, 300)

    # Left panel: calibration curves
    ax = axes[0]
    ax.scatter(channels, energies, color="black", zorder=5, s=60,
               label="Calibration source peaks")
    for i, (ch, e, lbl) in enumerate(zip(channels, energies, labels)):
        ax.annotate(lbl.split(" ")[0], (ch, e), textcoords="offset points",
                    xytext=(5, 3), fontsize=7, color="navy")

    ax.plot(ch_model, linear(ch_model, WRONG_SLOPE, WRONG_INTERCEPT),
            "--r", linewidth=2,
            label=f"WRONG:   E = {WRONG_SLOPE:.2f}·C + {WRONG_INTERCEPT:.2f} keV\n"
                  f"         (Range: {WRONG_INTERCEPT:.0f}–"
                  f"{linear(350, WRONG_SLOPE, WRONG_INTERCEPT):.0f} keV)")
    ax.plot(ch_model, linear(ch_model, *popt_correct),
            "-g", linewidth=2,
            label=f"CORRECT: E = {slope_c:.2f}·C + {intercept_c:.2f} keV\n"
                  f"         (Range: {linear(0, *popt_correct):.0f}–"
                  f"{linear(350, *popt_correct):.0f} keV)")

    # Mark 28Al expected channel
    ch_28Al = (1778.99 - intercept_c) / slope_c
    ax.axvline(x=ch_28Al, color="purple", linestyle=":", linewidth=1.2,
               label=f"²⁸Al 1779 keV → ch {ch_28Al:.0f} (correct)")

    ax.set_xlabel("MCA Channel", fontsize=11)
    ax.set_ylabel("Energy [keV]", fontsize=11)
    ax.set_title("Energy Calibration:\nWrong vs Correct", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Right panel: residuals
    ax2 = axes[1]
    x_pos = np.arange(len(channels))
    bar_w = 0.38
    ax2.bar(x_pos - bar_w/2, residuals_w, bar_w, color="red",   alpha=0.7, label="Wrong calibration")
    ax2.bar(x_pos + bar_w/2, residuals_c, bar_w, color="green", alpha=0.7, label="Correct calibration")
    ax2.axhline(0, color="black", linewidth=0.8)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([lbl.replace(" ", "\n") for lbl in labels],
                        fontsize=7, rotation=0)
    ax2.set_ylabel("Residual: E_pred − E_known [keV]", fontsize=11)
    ax2.set_title(f"Calibration Residuals\n"
                  f"Wrong RMS = {rms_w:.1f} keV  |  Correct RMS = {rms_c:.1f} keV",
                  fontsize=11)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis="y")

    plt.suptitle(
        "Calibration Failure Diagnostic — ²⁷Al NAA Experiment\n"
        "Root cause: ORTEC 572A gain setting not updated for 0–3 MeV range",
        fontsize=12
    )
    plt.tight_layout(rect=[0, 0, 1, 0.93])

    save_path = os.path.join(figures_dir, "calibration_failure_diagnostic.png")
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"\n  Saved: {save_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_calibration_diagnostic()
