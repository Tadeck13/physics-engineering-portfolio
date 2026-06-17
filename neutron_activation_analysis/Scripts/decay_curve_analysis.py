"""
decay_curve_analysis.py
=======================
Time-resolved decay curve visualization for neutron-activated ²⁷Al.

Loads the four integrated-count time-series (one per spectroscopic range),
normalizes each to its initial measurement, and plots all four on a shared
logarithmic y-scale. The log-linear appearance of the early portion of each
curve confirms exponential decay, while the late-time flattening onto a
non-zero background plateau motivates the truncated fitting window used in
half_life_extraction.py.

Physical context:
    Expected decay products from ²⁷Al irradiation with AmBe (0.025 eV–11 MeV):
      - ²⁷Al(n,γ)²⁸Al : t½ = 134.7 s,  γ at 1779 keV        (thermal/epithermal)
      - ²⁷Al(n,p)²⁷Mg  : t½ = 567.6 s,  γ at 844, 1014 keV  (fast neutrons > 3 MeV)

    Background contribution from the AmBe source itself (4.44 MeV from ⁹Be(α,n)¹²C*)
    and NORM sources causes the count rate to plateau at late times rather than
    reaching zero. This motivates fitting only the early, decay-dominated window.

Author : Tadeck Jones
Course : PHY3004W — Neutron Activation Analysis, UCT (2024)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# ---------------------------------------------------------------------------
# Range catalogue — one entry per spectroscopic range
# ---------------------------------------------------------------------------

# Note on energies: the originally reported energies (151.3, 172.9, 177.46,
# 181.18 keV) are artifacts of a miscalibrated amplifier (E = 0.11C + 146.62 keV).
# The calibration equation is physically inconsistent — it implies a maximum
# detectable energy of ~372 keV across 2048 channels, incompatible with
# calibration sources at 1332 keV. True peak energies are indeterminate without
# a valid calibration, but the temporal analysis below is calibration-independent.

RANGES = [
    {
        "file":  "Range_1_IntegratedCounts.txt",
        "label": "Range 1",
        "color": "steelblue",
        "fit_end_s": 390,          # last time point where decay dominates BG
        "t_half_expected_s": None,  # unclear match; slow component
    },
    {
        "file":  "Range_2_IntegratedCounts.txt",
        "label": "Range 2",
        "color": "darkorange",
        "fit_end_s": 570,
        "t_half_expected_s": 567.6,  # ²⁷Mg
    },
    {
        "file":  "Range_3_IntegratedCounts.txt",
        "label": "Range 3",
        "color": "forestgreen",
        "fit_end_s": 450,
        "t_half_expected_s": 567.6,  # ²⁷Mg
    },
    {
        "file":  "Range_4_IntegratedCounts.txt",
        "label": "Range 4",
        "color": "crimson",
        "fit_end_s": 390,
        "t_half_expected_s": 134.7,  # ²⁸Al
    },
]

# Reference isotope lines for visual comparison
REFERENCE_ISOTOPES = [
    {"name": "²⁸Al",  "t_half_s": 134.7,  "color": "red",    "linestyle": ":"},
    {"name": "²⁷Mg",  "t_half_s": 567.6,  "color": "purple", "linestyle": "--"},
]


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_decay_data(filepath):
    """
    Load a two-column time-series file (time [s] \\t integrated counts).

    Returns
    -------
    time   : ndarray  Time stamps [s]
    counts : ndarray  Integrated counts per time interval
    """
    data   = np.loadtxt(filepath, skiprows=1)
    time   = data[:, 0]
    counts = data[:, 1]
    return time, counts


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_all_decay_curves(ranges=None, data_dir=None, save_path=None):
    """
    Plot all decay curves (normalized) on a shared log scale.

    Parameters
    ----------
    ranges   : list of dict   Range catalogue (uses module-level RANGES if None)
    data_dir : str or None    Path to data/ directory
    save_path : str or None   Save figure if provided
    """
    if ranges is None:
        ranges = RANGES
    if data_dir is None:
        data_dir = os.path.join(os.path.dirname(__file__), "../data/")

    fig, ax = plt.subplots(figsize=(11, 6))

    for rng in ranges:
        fpath = os.path.join(data_dir, rng["file"])
        t, N  = load_decay_data(fpath)

        # Normalize to first measured point (t = 30 s)
        N_norm = N / N[0]

        ax.semilogy(t, N_norm,
                    color=rng["color"], linewidth=1.2, alpha=0.85,
                    label=rng["label"])

        # Mark the fit window cutoff
        ax.axvline(x=rng["fit_end_s"], color=rng["color"],
                   linewidth=0.7, linestyle="--", alpha=0.5)

    # Overlay reference exponential decay curves
    t_ref = np.linspace(0, 3600, 500)
    for iso in REFERENCE_ISOTOPES:
        lam   = np.log(2) / iso["t_half_s"]
        N_ref = np.exp(-lam * t_ref)
        ax.semilogy(t_ref, N_ref,
                    color=iso["color"], linewidth=1.5,
                    linestyle=iso["linestyle"],
                    label=f'{iso["name"]} reference  (t½ = {iso["t_half_s"]} s)',
                    alpha=0.7)

    ax.set_xlabel("Time [s]", fontsize=12)
    ax.set_ylabel("Normalized Count Rate  N(t) / N(30s)", fontsize=12)
    ax.set_title(
        "Time-Resolved Decay Curves — ²⁷Al Neutron Activation\n"
        "(log scale: linear portion = pure exponential decay; "
        "late plateau = background floor)",
        fontsize=11
    )
    ax.set_xlim(0, 3650)
    ax.set_ylim(0.1, 1.1)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, which="both", alpha=0.3)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(
        lambda y, _: f"{y:.1f}" if y >= 0.1 else f"{y:.2f}"
    ))

    # Annotate background plateau
    ax.annotate(
        "Background plateau\n(AmBe + NORM)",
        xy=(2500, 0.5), fontsize=9, color="gray",
        ha="center", style="italic"
    )

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved: {save_path}")
    plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    figures_dir = os.path.join(os.path.dirname(__file__), "../figures/")
    os.makedirs(figures_dir, exist_ok=True)

    print("\n=== Decay Curve Analysis — ²⁷Al Neutron Activation ===\n")
    print("Expected isotopes from ²⁷Al + AmBe irradiation:")
    for iso in REFERENCE_ISOTOPES:
        print(f"  {iso['name']}: t½ = {iso['t_half_s']} s")
    print()
    print("Plotting all 4 spectroscopic ranges...")

    plot_all_decay_curves(
        save_path=os.path.join(figures_dir, "decay_curves_all_ranges.png")
    )
