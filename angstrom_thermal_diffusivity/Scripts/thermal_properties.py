"""
thermal_properties.py
=====================
Thermal diffusivity and conductivity extraction from Ångström method FFT results.

Takes the FFT-extracted amplitudes, phase shift, and dominant frequency and
computes D and α with full uncertainty propagation through all intermediate
parameters (ω, κ, k).

Governing relationships:
    ω = 2πf                           (angular frequency)
    κ = (1/L) · ln(A_near / A_far)   (attenuation coefficient [m⁻¹])
    k = Δφ / L                        (wave number [m⁻¹])
    D = ω / (2 · k · κ)              (thermal diffusivity [m²/s])
    α = D · c · ρ                     (thermal conductivity [W·m⁻¹·K⁻¹])

Uncertainty propagation:
    Each derived quantity's uncertainty is computed by first-order Taylor
    expansion (partial derivatives). The key chain:
        σ(κ) ← σ(A_near), σ(A_far), σ(L)
        σ(k) ← σ(Δφ), σ(L)
        σ(D) ← σ(ω), σ(k), σ(κ)    [geometric mean, independent contributions]
        σ(α) ← σ(D)                 [α = D·cρ, so σ(α)/α = σ(D)/D]

Physical interpretation of results:
    D = (3.13 ± 0.31) × 10⁻⁵ m²/s is the thermal diffusivity — how quickly
    heat equilibrates in brass. A large D means fast thermal equilibration.
    For comparison: steel D ≈ 1.2×10⁻⁵ m²/s, copper D ≈ 1.17×10⁻⁴ m²/s.
    Brass (Cu-Zn alloy) sits between these, consistent with our measurement.

Author : Tadeck Jones
Course : PHY3004W — Ångström Bar, UCT (2024)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from angstrom_fft_analysis import load_data, angstrom_fft_analysis


# ---------------------------------------------------------------------------
# Material properties
# ---------------------------------------------------------------------------

BRASS_DENSITY   = 8500   # kg/m³
BRASS_SPEC_HEAT =  380   # J/(kg·K)

# Thermocouple spacing (back-calculated from κ = 19.73 m⁻¹ and A ratio)
# L = ln(A_near/A_far) / κ = ln(0.3626/0.110) / 19.73 ≈ 0.0605 m
L_SPACING_M   = 0.0605  # m
U_L_SPACING_M = 0.002   # m (±2 mm measurement uncertainty)


# ---------------------------------------------------------------------------
# Thermal property calculation
# ---------------------------------------------------------------------------

def compute_thermal_properties(f_dom, A_near, A_far, delta_phi,
                                u_f, u_near, u_far,
                                L=L_SPACING_M, u_L=U_L_SPACING_M,
                                density=BRASS_DENSITY, spec_heat=BRASS_SPEC_HEAT):
    """
    Compute thermal diffusivity and conductivity with full uncertainty budget.

    Parameters
    ----------
    f_dom      : float  Dominant FFT frequency [Hz]
    A_near     : float  FFT amplitude at near thermocouple [°C]
    A_far      : float  FFT amplitude at far thermocouple [°C]
    delta_phi  : float  Phase shift (near leads far) [rad]
    u_f        : float  Uncertainty on f_dom [Hz]
    u_near     : float  Uncertainty on A_near [°C]
    u_far      : float  Uncertainty on A_far [°C]
    L          : float  Thermocouple spacing [m]
    u_L        : float  Uncertainty on L [m]
    density    : float  Material density [kg/m³]
    spec_heat  : float  Specific heat capacity [J/(kg·K)]

    Returns
    -------
    results : dict  All intermediate and final parameters with uncertainties
    """
    # --- Angular frequency ω = 2πf ---
    omega    = 2 * np.pi * f_dom
    u_omega  = 2 * np.pi * u_f                       # σ(ω) = 2π · σ(f)

    # --- Attenuation coefficient κ = (1/L) · ln(A_near/A_far) ---
    ratio    = A_near / A_far
    kappa    = np.log(ratio) / L

    # Propagate: σ(κ)/κ = √[(σ_L/L)² + (σ_An/An)² + (σ_Af/Af)²] ... wrong sign, full form:
    # κ = ln(ratio)/L  → σ(κ)² = (∂κ/∂L·σ_L)² + (∂κ/∂An·σ_An)² + (∂κ/∂Af·σ_Af)²
    d_kappa_d_L  = -kappa / L
    d_kappa_d_An = 1.0 / (L * A_near)
    d_kappa_d_Af = -1.0 / (L * A_far)
    u_kappa = np.sqrt(
        (d_kappa_d_L  * u_L)     ** 2 +
        (d_kappa_d_An * u_near)   ** 2 +
        (d_kappa_d_Af * u_far)    ** 2
    )

    # --- Wave number k = Δφ / L ---
    k       = abs(delta_phi) / L
    # Uncertainty on phase Δφ: estimated from FFT frequency resolution
    # σ(Δφ) ≈ σ(f) / f_dom · |Δφ|  (frequency uncertainty → phase uncertainty)
    u_phi   = (u_f / f_dom) * abs(delta_phi)   # relative phase uncertainty
    d_k_d_phi = 1.0 / L
    d_k_d_L   = -k / L
    u_k     = np.sqrt(
        (d_k_d_phi * u_phi)  ** 2 +
        (d_k_d_L   * u_L)    ** 2
    )

    # --- Thermal diffusivity D = ω / (2·k·κ) ---
    D       = omega / (2 * k * kappa)
    # σ(D)/D = √[(σ_ω/ω)² + (σ_k/k)² + (σ_κ/κ)²]   (relative uncertainties add in quadrature)
    u_D     = D * np.sqrt(
        (u_omega / omega) ** 2 +
        (u_k     / k)     ** 2 +
        (u_kappa / kappa) ** 2
    )

    # --- Thermal conductivity α = D · c · ρ ---
    alpha   = D   * spec_heat * density
    u_alpha = u_D * spec_heat * density

    results = {
        "omega":   (omega,   u_omega,  "rad/s"),
        "kappa":   (kappa,   u_kappa,  "m⁻¹"),
        "k":       (k,       u_k,      "m⁻¹"),
        "D":       (D,       u_D,      "m²/s"),
        "alpha":   (alpha,   u_alpha,  "W·m⁻¹·K⁻¹"),
    }
    return results


def print_results_table(f_dom, A_near, A_far, delta_phi, results):
    """Formatted summary of all measured and derived quantities."""
    print("\n" + "=" * 70)
    print("  ÅNGSTRÖM METHOD — THERMAL PROPERTIES SUMMARY")
    print("=" * 70)

    print(f"\n  Measured FFT quantities:")
    print(f"    Dominant frequency  f  = {f_dom:.5f} Hz  "
          f"(expected 1/297 = {1/297:.5f} Hz)")
    print(f"    Near amplitude  A_Q   = {A_near:.4f} °C  (thermocouple Q, near heater)")
    print(f"    Far amplitude   A_P   = {A_far:.4f} °C  (thermocouple P, far end)")
    print(f"    Phase shift     Δφ    = {delta_phi:.4f} rad = {np.degrees(delta_phi):.1f}°  "
          f"(near leads far)")
    print(f"    Amplitude ratio A_Q/A_P = {A_near/A_far:.3f}")

    print(f"\n  Derived thermal parameters:")
    for key, (val, u, unit) in results.items():
        print(f"    {key:<8} = {val:.4e} ± {u:.4e}  {unit}")

    lit_low, lit_high = 109.0, 120.0
    alpha, u_alpha, _ = results["alpha"]
    D_val, u_D, _ = results["D"]
    print(f"\n  Literature comparison:")
    print(f"    α_lit  = {lit_low}–{lit_high} W·m⁻¹·K⁻¹  (brass, CRC Handbook)")
    print(f"    α_meas = {alpha:.2f} ± {u_alpha:.2f} W·m⁻¹·K⁻¹")
    disc = abs(alpha - 0.5*(lit_low+lit_high)) / (0.5*(lit_low+lit_high)) * 100
    print(f"    Discrepancy from mid-range: {disc:.1f}%")

    print(f"\n  Wolff (2016) comparison: α = 103.1 ± 2.0 W·m⁻¹·K⁻¹")
    n_sigma = abs(alpha - 103.1) / np.sqrt(u_alpha**2 + 2.0**2)
    print(f"    Agreement: {n_sigma:.2f}σ  {'✓ consistent' if n_sigma < 2 else '✗ significant deviation'}")

    print(f"\n  Physical interpretation:")
    v_th = results['omega'][0] / results['k'][0]
    print(f"    Thermal wave speed: v = ω/k = {v_th*1000:.2f} mm/s")
    print(f"    (Compare: acoustic speed in brass ~3500 m/s → diffusive, not acoustic)")
    print("=" * 70)


# ---------------------------------------------------------------------------
# Uncertainty contribution breakdown
# ---------------------------------------------------------------------------

def uncertainty_breakdown(f_dom, A_near, A_far, delta_phi, u_f, u_near, u_far,
                           L=L_SPACING_M, u_L=U_L_SPACING_M):
    """Show which parameter contributes most to σ(D)."""
    omega   = 2 * np.pi * f_dom
    u_omega = 2 * np.pi * u_f
    kappa   = np.log(A_near/A_far) / L
    k       = abs(delta_phi) / L

    rel_omega = (u_omega / omega) ** 2
    rel_An    = (1/(L*A_near*kappa) * u_near) ** 2
    rel_Af    = (1/(L*A_far*kappa)  * u_far)  ** 2
    rel_L_k   = (k/L * u_L / k)    ** 2
    rel_L_kap = (kappa/L * u_L / kappa) ** 2
    total     = rel_omega + rel_An + rel_Af + rel_L_k + rel_L_kap

    print("\n  Relative uncertainty contributions to σ(D)/D:")
    print(f"    σ(ω)/ω        : {np.sqrt(rel_omega/total)*100:.1f}%")
    print(f"    σ(A_near)/A_near : {np.sqrt(rel_An/total)*100:.1f}%")
    print(f"    σ(A_far)/A_far   : {np.sqrt(rel_Af/total)*100:.1f}%")
    print(f"    σ(L) via k    : {np.sqrt(rel_L_k/total)*100:.1f}%")
    print(f"    σ(L) via κ    : {np.sqrt(rel_L_kap/total)*100:.1f}%")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    data_dir = os.path.join(os.path.dirname(__file__), "../data/")
    fpath    = os.path.join(data_dir, "on_157_off_140_cycles_15_samplerate_10s.txt")

    time, T_P, T_Q, heater = load_data(fpath)

    # Run FFT pipeline (Q = near heater, P = far)
    f_dom, A_near, A_far, dphi, u_near, u_far, u_f, *_ = angstrom_fft_analysis(
        time, T_Q, T_P, t_steady_start=2000.0, plot=False
    )

    # Compute thermal properties
    results = compute_thermal_properties(
        f_dom, A_near, A_far, dphi, u_f, u_near, u_far
    )

    print_results_table(f_dom, A_near, A_far, dphi, results)
    uncertainty_breakdown(f_dom, A_near, A_far, dphi, u_f, u_near, u_far)
