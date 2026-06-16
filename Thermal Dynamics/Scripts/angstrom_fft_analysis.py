"""
angstrom_fft_analysis.py
=========================
FFT-based thermal wave analysis for the Ångström bar experiment.

Extracts the dominant frequency, amplitude attenuation, and phase shift of a
periodic heat wave propagating through a brass rod from two thermocouple records.

Physical model:
    The heat equation in a rod subject to sinusoidal heating at x=0 has the
    travelling-wave solution:
        T(x, t) = A₀ · exp(-κx) · cos(ωt - kx + φ₀)

    where κ is the attenuation coefficient and k = ω/v_phase is the wave number.
    Thermal diffusivity: D = ω/(2kκ).

FFT strategy:
    1. Select steady-state window (after transient warm-up)
    2. Detrend (subtract mean) to remove DC offset
    3. Interpolate to uniform time grid (accommodates irregular sampling)
    4. Apply FFT and identify dominant frequency bin
    5. Fit Gaussian to dominant peak → sub-bin amplitude precision and σ
    6. Extract phase from complex FFT coefficients at dominant frequency

Key outputs:
    f_dom       — dominant frequency [Hz]
    A_near, A_far — FFT amplitudes at near (Q) and far (P) thermocouples [°C]
    delta_phi   — phase shift Δφ = φ_Q - φ_P [rad] (negative = far lags near)
    sigma_near, sigma_far — uncertainties from Gaussian peak fit [°C]

Author : Tadeck Jones
Course : PHY3004W — Ångström Bar, UCT (2024)
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_data(filepath):
    """
    Load Ångström bar temperature data.

    Expected columns: Time [s], Heater status [0/1], Temp_P [°C], Temp_Q [°C]

    Returns
    -------
    time    : ndarray  Time stamps [s]
    T_P     : ndarray  Temperature at far thermocouple P [°C]
    T_Q     : ndarray  Temperature at near thermocouple Q [°C]
    heater  : ndarray  Heater ON/OFF flag [0 or 1]
    """
    data   = np.loadtxt(filepath, skiprows=1)
    time   = data[:, 0]
    heater = data[:, 1]
    T_P    = data[:, 2]
    T_Q    = data[:, 3]
    return time, T_P, T_Q, heater


# ---------------------------------------------------------------------------
# Gaussian model for sub-bin amplitude fitting
# ---------------------------------------------------------------------------

def _gaussian(x, amplitude, centroid, sigma):
    return amplitude * np.exp(-0.5 * ((x - centroid) / sigma) ** 2)


def _gaussian_peak_fit(frequencies, amplitudes):
    """
    Fit a Gaussian to a peak in the FFT spectrum.

    Returns the peak amplitude and its 1-σ uncertainty from the covariance.
    """
    try:
        dom_idx = np.argmax(amplitudes)
        p0 = [amplitudes[dom_idx], frequencies[dom_idx],
              (frequencies[-1] - frequencies[0]) / 4]
        popt, pcov = curve_fit(_gaussian, frequencies, amplitudes, p0=p0, maxfev=10_000)
        return abs(popt[0]), np.sqrt(pcov[0, 0])
    except RuntimeError:
        return amplitudes[dom_idx], 0.0


# ---------------------------------------------------------------------------
# Core FFT analysis
# ---------------------------------------------------------------------------

def angstrom_fft_analysis(time, T_near, T_far, t_steady_start=2000.0,
                           plot=False, label=""):
    """
    Full FFT-based Ångström thermal wave analysis.

    Parameters
    ----------
    time            : ndarray  Time [s]
    T_near          : ndarray  Temperature at NEAR thermocouple (Q) [°C]
    T_far           : ndarray  Temperature at FAR thermocouple (P) [°C]
    t_steady_start  : float    Start of steady-state window [s]
    plot            : bool     Show diagnostic plots
    label           : str      Dataset label

    Returns
    -------
    f_dom       : float  Dominant frequency [Hz]
    A_near      : float  FFT amplitude at near thermocouple [°C]
    A_far       : float  FFT amplitude at far thermocouple [°C]
    delta_phi   : float  Phase shift Δφ = φ_near - φ_far [rad]  (positive → near leads)
    u_near      : float  Uncertainty on A_near from Gaussian fit [°C]
    u_far       : float  Uncertainty on A_far from Gaussian fit [°C]
    u_f         : float  Frequency resolution [Hz]  = 1 / T_total
    freqs_pos   : ndarray  Positive-frequency array [Hz]
    amp_near    : ndarray  Amplitude spectrum for near thermocouple
    amp_far     : ndarray  Amplitude spectrum for far thermocouple
    """
    # --- Steady-state selection and detrend ---
    ss      = time >= t_steady_start
    t_ss    = time[ss]
    Tn_ss   = T_near[ss] - np.mean(T_near[ss])   # near (Q)
    Tf_ss   = T_far[ss]  - np.mean(T_far[ss])    # far  (P)

    T_total = t_ss[-1] - t_ss[0]   # total steady-state duration [s]
    u_f     = 1.0 / T_total         # frequency resolution [Hz]

    # --- Interpolate to uniform grid (handles occasional missing samples) ---
    dt_med  = np.median(np.diff(t_ss))
    t_uni   = np.arange(t_ss[0], t_ss[-1] + dt_med, dt_med)
    Tn_uni  = interp1d(t_ss, Tn_ss, kind="linear", fill_value="extrapolate")(t_uni)
    Tf_uni  = interp1d(t_ss, Tf_ss, kind="linear", fill_value="extrapolate")(t_uni)

    # --- FFT ---
    N         = len(t_uni)
    freqs     = np.fft.fftfreq(N, d=dt_med)
    FFT_near  = np.fft.fft(Tn_uni)
    FFT_far   = np.fft.fft(Tf_uni)

    # Positive-frequency one-sided spectrum
    pos       = freqs > 0
    freqs_pos = freqs[pos]
    amp_near  = np.abs(FFT_near[pos]) / N * 2
    amp_far   = np.abs(FFT_far[pos])  / N * 2

    # --- Dominant frequency and phase shift ---
    dom       = np.argmax(amp_near + amp_far)    # peak in sum of both spectra
    f_dom     = freqs_pos[dom]
    phi_near  = np.angle(FFT_near[pos][dom])
    phi_far   = np.angle(FFT_far[pos][dom])
    delta_phi = phi_near - phi_far               # near leads far by this amount

    # --- Sub-bin amplitude fitting (Gaussian on ±5 bins around peak) ---
    roi       = slice(max(0, dom - 5), min(len(freqs_pos), dom + 6))
    A_near, u_near = _gaussian_peak_fit(freqs_pos[roi], amp_near[roi])
    A_far,  u_far  = _gaussian_peak_fit(freqs_pos[roi], amp_far[roi])

    if plot:
        _plot_fft(freqs_pos, amp_near, amp_far, A_near, A_far,
                  f_dom, u_f, label)

    return f_dom, A_near, A_far, delta_phi, u_near, u_far, u_f, freqs_pos, amp_near, amp_far


def _plot_fft(freqs_pos, amp_near, amp_far, A_near, A_far,
              f_dom, u_f, label):
    """Diagnostic two-sided FFT amplitude spectrum plot."""
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(freqs_pos,  amp_near, "r-", linewidth=0.9, label="Amplitude Q (near)")
    ax.plot(freqs_pos,  amp_far,  "b-", linewidth=0.9, label="Amplitude P (far)")
    ax.plot(-freqs_pos, amp_near, "r-", linewidth=0.9, alpha=0.5)
    ax.plot(-freqs_pos, amp_far,  "b-", linewidth=0.9, alpha=0.5)
    ax.axvline(f_dom, color="gray", linestyle="--", linewidth=0.8,
               label=f"f_dom = {f_dom:.5f} Hz")
    ax.set_xlabel("Frequency (Hz)", fontsize=11)
    ax.set_ylabel("Amplitude (°C)", fontsize=11)
    ax.set_title(f"FFT Amplitude Spectrum — Ångström Bar  {label}", fontsize=11)
    ax.legend(fontsize=9)
    ax.set_xlim(-0.05, 0.05)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    data_dir = os.path.join(os.path.dirname(__file__), "../data/")
    fpath    = os.path.join(data_dir, "on_157_off_140_cycles_15_samplerate_10s.txt")

    print("=" * 65)
    print("  ÅNGSTRÖM BAR — FFT ANALYSIS")
    print("=" * 65)

    time, T_P, T_Q, heater = load_data(fpath)

    print(f"\n  Dataset: {len(time)} samples, t = [{time[0]:.0f}, {time[-1]:.0f}] s")
    print(f"  Heating period: T_on=157s + T_off=140s = 297s → f={1/297:.6f} Hz")

    f_dom, A_near, A_far, dphi, u_near, u_far, u_f, freqs, amp_n, amp_f = angstrom_fft_analysis(
        time, T_Q, T_P, t_steady_start=2000.0, plot=True, label="(10s sample rate)"
    )

    print(f"\n  FFT Results:")
    print(f"    Dominant frequency : {f_dom:.6f} ± {u_f:.6f} Hz")
    print(f"    Amplitude Q (near) : {A_near:.4f} ± {u_near:.4f} °C")
    print(f"    Amplitude P (far)  : {A_far:.4f}  ± {u_far:.4f}  °C")
    print(f"    Phase shift Δφ     : {dphi:.4f} rad  = {np.degrees(dphi):.1f}°")
    print(f"    Amplitude ratio    : A_near/A_far = {A_near/A_far:.3f}")
