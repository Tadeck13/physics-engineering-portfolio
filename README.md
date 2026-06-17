# Physics Engineering Portfolio

**Tadeck Jones** · BSc Physics, University of Cape Town · MeASURe Research Unit  
*Experimental Nuclear & Radiation Physics | Scientific Instrumentation | Scientific Python*

> This portfolio translates laboratory experiments and research into an industry-readable format: each project leads with its quantified result, documents the instrumentation chain, explains the data pipeline, and identifies the root cause of every systematic error encountered. All analyses are fully reproducible from the raw data files included here.

---

## Projects at a Glance

| # | Project | Technique | Key Result | Agreement with Reference |
|---|---------|-----------|------------|--------------------------|
| [01](./gamma_ray_spectroscopy/) | **γ-Ray Spectroscopy** | NaI(Tl) + Gaussian peak fitting | 7.0% energy resolution @ 661.66 keV | ✅ Within manufacturer spec (≤ 8.5%) |
| [02](./neutron_activation_analysis/) | **Neutron Activation Analysis** | AmBe irradiation + linearised WLS decay | t½ = 173 ± 9 s (²⁸Al); 462 ± 31 s (²⁷Mg) | ✅ Consistent with ²⁸Al (134.7 s), ²⁷Mg (567.6 s) |
| [03](./pept_spatial_resolution/) | **PEPT Spatial Resolution** | GATE Monte Carlo + skew Gaussian | μ_δ(R) = 0.21 mm; breakdown below A = 0.5 mm | ✅ Matches Leadbeater et al. ~0.5 mm limit |
| [04](balmer_spectroscopy/) | **Balmer Series / Rydberg Constant** | Optical spectroscopy + WLS | R = (1.09 ± 0.12) × 10⁷ m⁻¹ | ✅ 0.6% from NIST value |
| [05](./angstrom_thermal_diffusivity/) | **Thermal Diffusivity (Ångström)** | Periodic heating + FFT phase analysis | D = (3.13 ± 0.31) × 10⁻⁵ m²/s | ✅ Matches Wolff (2016) to 0.5σ |

---

## About This Portfolio

### What it covers

These five experiments span four distinct measurement domains — nuclear, optical, thermal, and particle-tracking simulation — all performed within one academic year at UCT. Each project demonstrates a complete measurement pipeline: hardware chain → raw data → statistical analysis → quantified result → uncertainty budget → failure-mode diagnosis.

The portfolio is structured so that an industry reviewer can understand the scope and quality of each project in under two minutes per page. Every project README follows the same format:

- **Result first** — the key number appears in the executive summary, before any methodology
- **Hardware chain** — explicit signal flow from source to data file
- **Data pipeline** — step-by-step, with transformations named and quantified
- **Engineering insight** — the one non-obvious finding that demonstrates analytical depth
- **Failure mode** — what went wrong, why, and what it cost in precision
- **Reproducible code** — scripts that regenerate every figure and number from the raw data

### The common thread

Every project required diagnosing which uncertainty source **dominated** the measurement budget — and in most cases, that source was not the obvious one:

| Project | Expected dominant uncertainty | Actual dominant uncertainty |
|---------|------------------------------|-----------------------------|
| P01 γ-ray spectroscopy | Gaussian fit of photopeak | Amplifier BLR baseline offset (+23 channels) |
| P02 Neutron activation | Decay curve fit | Manual transfer delay (~60–120 s of missed data) |
| P03 PEPT resolution | Statistical scatter in tracked positions | Physical averaging window (254-LOR sample size > circle circumference) |
| P04 Rydberg constant | Gaussian fit of spectral line | Single-point HeNe calibration offset (22 Å, **99.2%** of total uncertainty) |
| P05 Thermal diffusivity | Thermocouple noise | Transient warm-up contamination of FFT spectrum |

In each case, the standard "improve the fit" instinct would have been wasted effort. Identifying the correct dominant term is the skill this portfolio is designed to demonstrate.

---

## Technical Stack

### Instrumentation

| Instrument | Projects | Purpose |
|------------|----------|---------|
| NaI(Tl) scintillation detector | P01, P02 | γ-ray photon detection |
| ORTEC NIM electronics (113 preamp, 572A amp, 556H HV) | P01, P02 | Signal conditioning and shaping |
| MCA (1024 ch / 2048 ch) | P01, P02 | Pulse height spectroscopy |
| AmBe neutron source (2.2 GBq) | P02 | Neutron activation |
| Heath EU-700 Czerny-Turner monochromator | P04 | Optical wavelength selection |
| HeNe laser (6328 Å) | P04 | Spectrometer calibration reference |
| Photo-multiplier tube + NIM amplifier | P04 | Photon counting |
| Brass rod + thermocouples | P05 | Thermal wave propagation measurement |
| Siemens ECAT "EXACT3D" HR++ PET scanner | P03 | Positron annihilation γ-ray detection |

### Software & Analysis

```
Language  : Python 3.11
Libraries : NumPy · SciPy (curve_fit, fft, interp1d) · Matplotlib · Pandas ·
Simulation: GATE (Geant4 Application for Tomographic Emission), 50-core parallel
Fitting   : Symmetric Gaussian · Skew Gaussian · Weighted linear regression (WLS)
            Linearised exponential decay · FFT peak sub-bin fitting
```

### Statistical methods used (across all projects)

- **Weighted least squares (WLS)** — with proper 2×2 covariance matrix and off-diagonal propagation (P01, P02, P04)
- **Gaussian peak fitting** — photopeak centroids, spectral line centroids, FFT amplitude peaks (P01, P04, P05)
- **Skew Gaussian fitting** — strictly non-negative residual distributions in PEPT tracking (P03)
- **Linearised exponential decay** — Poisson-weighted regression on ln(N/N₀) = −λt (P02)
- **FFT phase extraction** — complex FFT coefficients → amplitude attenuation + phase shift → D, α (P05)
- **Full uncertainty budget propagation** — first-order Taylor expansion through every intermediate variable, with percentage contribution table per project

---

## Repository Structure

```
physics-engineering-portfolio/
│
├── README.md                           ← You are here
├── requirements.txt                    ← Python 3.11+ dependencies
├── templates/
│   └── PROJECT_README_TEMPLATE.md     ← Blank template for new projects
│
├── 01_gamma_ray_spectroscopy/
│   ├── README.md
│   ├── code/    gaussian_peak_fitting.py · energy_calibration.py · pipeline.py
│   ├── data/    137Cs.txt · 22Na.txt · 54Mn.txt · 57Co.txt · 133Ba.txt · Co60.txt
│   └── figures/ Eight Gaussian fit plots + calibration curve
│
├── 02_neutron_activation_analysis/
│   ├── README.md
│   ├── code/    decay_curve_analysis.py · half_life_extraction.py · calibration_failure_diagnostic.py
│   ├── data/    Range_1–4_IntegratedCounts.txt · Peak1–4_count_range.txt · calibration spectra
│   └── figures/ Energy_cal.png · Decay.png (two lab-generated figures)
│
├── 03_pept_spatial_resolution/
│   ├── README.md
│   ├── code/    skew_gaussian_residual_analysis.py · creationOfCircles.py · CircleMotionSimulation.py
│   ├── data/    optTrackingParams.csv
│   └── figures/ Thirteen figures: XY scatter, Lagrangian time series, residual PMFs, resolution plot
│
├── 04_balmer_spectroscopy/
│   ├── README.md
│   ├── code/    balmer_gaussian_fitting.py · spectrometer_calibration.py · rydberg_constant_extraction.py
│   ├── data/    H_alpha/beta/gamma/delta.txt · HeNe_CorrectRun_Final.txt · Beta_150micron.txt
│   └── figures/ Four Gaussian fits + HeNe calibration + Rydberg WLS fit
│
└── 05_angstrom_thermal_diffusivity/
    ├── README.md
    ├── code/    angstrom_fft_analysis.py · thermal_properties.py
    ├── data/    on_157_off_140_cycles_15_samplerate_10s.txt · sample_rate_5s.csv
    └── figures/ FFT spectrum · Gaussian fit · Steady-state · Full transient+steady
```

---

## Quickstart

**Requirements:**
```bash
pip install -r requirements.txt
```

**Run any project's full analysis pipeline:**
```bash
# Project 01 — γ-ray spectroscopy
python 01_gamma_ray_spectroscopy/code/spectral_analysis_pipeline.py

# Project 02 — Neutron activation
python 02_neutron_activation_analysis/code/half_life_extraction.py

# Project 04 — Rydberg constant
python 04_balmer_spectroscopy/code/rydberg_constant_extraction.py

# Project 05 — Thermal diffusivity
python 05_angstrom_thermal_diffusivity/code/thermal_properties.py
```

Each pipeline loads data from the adjacent `data/` folder and reproduces every quantitative result and figure from the corresponding README.

> **Note on Project 03:** The GATE Monte Carlo simulation infrastructure was developed by R. Perin et al. at UCT/iThemba LABS. The residual analysis code (`skew_gaussian_residual_analysis.py`) and the pre-computed results in `optTrackingParams.csv` are independently reproducible. Full simulation re-runs require a GATE installation and the HR++ geometry files.

---

## Key Results Summary

```
P01  NaI(Tl) energy resolution        :  7.0%          @ 661.66 keV   [spec: ≤ 8.5%]
P01  Calibration equation             :  C = 0.32·Eγ + 23.21          [7 isotopes, 88–1332 keV]
P01  Unknown source energy            :  714.14 ± σ keV                [tentative: not ¹³⁷Cs]

P02  ²⁸Al half-life (measured)        :  173  ±  9 s                   [literature: 134.7 s]
P02  ²⁷Mg half-life (measured)        :  462  ± 31 s                   [literature: 567.6 s]
P02  Calibration failure              :  E = 0.11C + 146.62 keV        [factor-of-14 slope error]

P03  PEPT tracking floor (μ_δ(R))     :  0.21 mm                       [constant for A ≥ 0.5 mm]
P03  PEPT tracking floor (σ_δ(R))     :  0.30 mm                       [constant for A ≥ 0.5 mm]
P03  Spatial resolution limit         :  A = 0.10–0.50 mm              [breakdown interval]

P04  Rydberg constant                 :  (1.09 ± 0.12) × 10⁷ m⁻¹      [NIST: 1.0974 × 10⁷ m⁻¹]
P04  Dominant uncertainty source      :  HeNe calibration offset        [99.2% of combined u(λ)]
P04  Spectrometer offset              :  +22.13 Å                       [single-point correction]

P05  Thermal diffusivity (brass)      :  (3.13 ± 0.31) × 10⁻⁵ m²/s    [literature: ~3.3 × 10⁻⁵]
P05  Thermal conductivity (brass)     :  101.9 ± 10.0 W·m⁻¹·K⁻¹       [literature: 109–120]
P05  Dominant frequency (measured)    :  0.00337 Hz                     [expected: 1/297 = 0.00337 Hz]
```

---

## Contact

**Tadeck Jones**

[LinkedIn](www.linkedin.com/in/tadeckjones) · [Email](tadeckjones13@gmail.com) · [GitHub](https://github.com/Tadeck13)

*Work completed as part of laboratory and research at the Department of Physics, University of Cape Town (2024).*

---

*Python 3.11+ · NumPy · SciPy · Matplotlib · GATE/Geant4*
