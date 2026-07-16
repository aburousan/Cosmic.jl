# Cosmic

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://aburousan.github.io/Cosmic.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://aburousan.github.io/Cosmic.jl/dev/)
[![Build Status](https://github.com/aburousan/Cosmic.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aburousan/Cosmic.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/aburousan/Cosmic.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/aburousan/Cosmic.jl)

[Cosmic](https://github.com/aburousan/Cosmic.jl) is a cosmology library in pure Julia,
following Baumann's *Cosmology* and Dodelson & Schmidt's *Modern Cosmology*.
It grew out of a cosmology course at the Indian Association for the Cultivation of Science.

## The idea

A cosmology is a **list of species**, each of which knows only its own ρ(a):

```julia
E(a)² = Σ_s ρ_s(a) / ρ_c0
```

Everything above that — the expansion rate, every distance, the thermal history —
is generic in the composition. Adding a new component costs one method, not a new
struct in every file. That is what lets massive neutrinos, arbitrary dark energy,
and extra relativistic relics coexist without special-casing.

## Installation

Requires Julia 1.9+.

```
pkg> add https://github.com/aburousan/Cosmic.jl
```

## Usage

```julia
using Cosmic

c = cosmology()                      # Planck 2018 ΛCDM
age(c)                               # 13.787 Gyr
luminosity_distance(c, 1.0)          # Mpc
z_equality(c)                        # 3388

# Things a fixed Ω_m/Ω_r/Ω_Λ struct cannot express:
cosmology(m_ν = [0.1, 0.1, 0.1])     # degenerate massive neutrinos
cosmology(w0 = -0.9, wa = -0.3)      # CPL dark energy
cosmology(w_func = a -> -1 + 0.1*(1-a)^2)   # arbitrary w(a)
cosmology(Ω_k = 0.01)                # open universe

# Thermal history
r = recombination(c)
x_e(r, 1100.0)                       # free-electron fraction
z_star(r)                            # 1089.8  — last scattering
r_drag(r)                            # 147.2 Mpc — the BAO standard ruler
```

## Physics

| Layer | Contents |
|---|---|
| `constants.jl` | CODATA constants; internal Mpc / c=1 unit system |
| `species.jl` | γ, massless ν, **massive ν** (full Fermi-Dirac phase-space integral), CDM, baryons, curvature, dark energy with arbitrary w(a) |
| `background.jl` | E(a), conformal time, all distances with correct curvature, ages |
| `bbn.jl` | Big Bang nucleosynthesis: 63-reaction / 12-nuclide network, n↔p weak rates derived from the Born integral with Coulomb + radiative + finite-mass + thermal-QED corrections, Monte-Carlo error propagation (`bbn_mc`). Y_p and D/H are **solved for**, not read from a table |
| `hyrec.jl` | Hydrogen recombination at HyRec-2020 accuracy: effective four-level atom from tabulated multilevel rates, SWIFT radiative-transfer correction |
| `thermodynamics.jl` | Recombination (RECFAST 1.5 default, `hydrogen = :hyrec` for HyRec-2020), exotic energy injection (DM annihilation/decay), x_e, τ, visibility, sound horizon |
| `nuclear_deposition.jl` | Post-BBN ³H beta-decay and ⁷Be recombination/electron-capture histories sourced by the solved nuclear abundances; electromagnetic energy only, with neutrino energy excluded |
| `perturbations.jl` | Einstein-Boltzmann hierarchy in conformal Newtonian gauge: photon polarization, free-streaming massless neutrinos, and the momentum-resolved massive-ν hierarchy |
| `gauges.jl` | Transformations to synchronous / comoving / uniform-density / spatially-flat gauges; gauge invariants ℛ, ζ, Bardeen potentials |
| `observables.jl` | Transfer functions, P(k), σ₈, growth D(z) and f(z) |
| `cmb.jl` | C_ℓ^TT/EE/TE by line-of-sight integration (Seljak & Zaldarriaga) |
| `distortions.jl` | CMB spectral distortions μ, y, and the residual (r-type) signal, by the exact Green's-function / PCA method (Chluba & Jeong 2014) |
| `thermalization.jl` | μ and y **derived from scratch**: solves the photon Boltzmann equation with the exact Klein–Nishina ⊗ Maxwell–Jüttner redistribution kernel, exact double Compton, and exact free–free Gaunt factors — no Green's-function table |

The CMB source function is decomposed into its four physical terms — Sachs-Wolfe,
integrated Sachs-Wolfe, Doppler, and the polarization quadrupole — each with a
weight you can set to zero to see what it contributes:

```julia
cmb_spectra(c, r; dop = 0.0)   # what the spectrum looks like with no Doppler term
```

The perturbation equations as implemented — including a derivation of the adiabatic
initial conditions, so they can be checked rather than taken on faith — are written
out in [`docs/src/perturbation_equations.md`](docs/src/perturbation_equations.md).

Massive neutrinos are integrated over the real Fermi-Dirac distribution rather than
interpolated between a⁻⁴ and a⁻³ power laws. As a check, the code reproduces the
standard Ω_ν h² = Σm_ν / 93.14 eV relation to 0.1% without that number appearing
anywhere in it.

Recombination follows RECFAST 1.5, including the corrections calibrated against the
multi-level-atom codes: a retuned hydrogen fudge factor and a double-Gaussian
correction to the Lyman-α escape rate (matched to HyRec and CosmoRec), and for
helium a Sobolev escape probability, H I continuum opacity (Switzer & Hirata 2008),
and the 2³P₁ → 1¹S₀ intercombination channel. Passing
`fudge_H = 1.14, gaussian_H = false, sobolev_He = false, triplets = false`
recovers the 1999 RECFAST, which makes the size of each correction easy to see.

`recombination(c; hydrogen = :hyrec)` replaces the fudged three-level hydrogen atom
with HyRec-2020's model: an effective four-level atom whose recombination
coefficients α_2s(T_M, T_R), α_2p(T_M, T_R) come from the full l-resolved multilevel
calculation (tabulated in `data/hyrec/`), photoionization by detailed balance, and
the SWIFT fitted correction for the Lyman-α radiative transfer. It reproduces the
HyRec−RECFAST shift CLASS sees (Δz_star = −0.03 here vs −0.025 in CLASS).

`recombination(c; injection = Injection(p_ann = ...))` adds exotic energy injection —
dark-matter annihilation and/or decay heating and ionizing the gas, in CLASS's
parameterisation with Chen–Kamionkowski channel splitting. The extra ionization
floor after freeze-out and the corresponding τ shift are how the CMB constrains
⟨σv⟩/m_DM; the deposition enters the x_e, T_b equations and, through them,
everything downstream.

## Validation

**Against CAMB 1.6.0**, run with an identical cosmology. Not a table of round
numbers — this is what CAMB actually outputs:

| Quantity | Cosmic | CAMB | diff |
|---|---|---|---|
| Age | 13.804049 Gyr | 13.803485 | 0.004% |
| z_eq | 3388.3696 | 3388.3695 | 0.000% |
| k_eq | 0.010340352 | 0.010340352 | 0.000% |
| z_* | 1089.790 | 1089.796 | 0.001% |
| z_drag | 1059.953 | 1059.991 | 0.004% |
| r_s(z_drag) | 147.233 Mpc | 147.230 | 0.003% |
| 100 θ_* | 1.039733 | 1.039727 | 0.001% |
| σ₈ (massless ν) | 0.822268 | 0.822326 | 0.007% |
| P(k), 0.1 < k < 2 | — | — | **≤ 0.04%** |
| P(k), k = 1e-3 | — | — | 0.9% |
| **C_ℓ^TT**, ℓ = 2…1400 | — | — | **mean 0.999, scatter 0.4%** |
| **C_ℓ^TT**, ℓ = 540…1400 | — | — | **0.999 flat (0.1%)** |
| **C_ℓ^EE**, ℓ = 2…1400 | — | — | **1.000, scatter 0.6%** |
| Acoustic peaks | ℓ = 220, 540, 820 | 221, 537, 816 | — |

The perturbation layer is separately checked against limits it cannot fake:

- The superhorizon potential in radiation domination reproduces the analytic
  Φ = (1 + 2f_ν/5)/(3/2 + 2f_ν/5) to four digits.
- Φ drops by **9/10** across matter-radiation equality — the canonical test of a
  Boltzmann hierarchy, which fails if the photon, neutrino or metric equations
  have a sign error anywhere.
- The Poisson equation k²Φ = -(3/2)H₀²Ω_m δ_m/a holds sub-horizon.
- P(k) turns over at k_eq = a_eq H(a_eq), predicted independently by the background.

### Curved geometry and tensors

A non-flat model never falls back to flat Bessel functions, at any wavenumber.
The radial functions are the true hyperspherical ones (Frobenius series near the
origin, Miller/forward recurrence away from it), scalars live at k² = q² − K and
tensors at k² = q² − 3K, and a closed universe gets its actual discrete eigenmode
sum with the exact weight q√K/k² — not a quadrature pretending the spectrum is
continuous. CLASS switches to rescaled flat Bessels above ν = 4000 as a precision
option; there is deliberately no such branch here.

At matched, converged k_max:

| | reference | agreement |
|---|---|---|
| open scalar TT/TE/EE | CLASS 3.3.4 | ≤0.6% (TE/EE mostly ≤0.1%) |
| closed scalar TT/TE/EE | CAMB 1.6.6 | TT ≤0.4%, TE/EE ≤0.7% |
| flat + open tensor TT/TE/EE/BB | CLASS 3.3.4 | ≤1% |
| closed tensor TT/TE/EE/BB | CAMB 1.6.6 | ≤3% (≤0.3% at ℓ ≤ 10) |

CAMB is the closed reference because CLASS 3.3.4 itself crashes on closed models
(a harmonic-module spline index, then heap corruption).

### CMB lensing

`lensed_cmb_spectra(c, rec)` returns the lensed TT/TE/EE/BB and, alongside them,
the lensing-potential spectrum C_ℓ^φφ and the ISW–lensing cross C_ℓ^Tφ — all
from one scalar mode sweep. The potential is the line-of-sight integral of the
Weyl potential φ+ψ against the geometric lensing kernel; the lensed spectra come
from the full-sky correlation-function method of Challinor & Lewis (2005), the
same algorithm as CLASS's `lensing.c` — the deflection variance σ²(β) is kept
inside the exponential (not expanded in C_ℓ^φφ), which is what makes the lensed
B-modes and the smoothed damping tail come out right.

Validated against a live CLASS 3.3.4 lensed run (flat, massless ν, HyRec,
matched reionization) at converged kmax and grid:

| | agreement vs CLASS |
|---|---|
| lensed TT, ℓ = 30…2000 | 0.3% |
| lensed EE, ℓ = 30…2000 | 0.4% |
| lensed TE, ℓ = 30…2000 | \|ΔTE\|/√(TT·EE) ≤ 0.0023 |
| C_ℓ^φφ, ℓ ≤ 1000 | ≤ 0.9% |

Two things need their own resolution and are easy to get wrong. C_ℓ^φφ at high ℓ
draws on much larger k than the temperature transfer does — the kernel reaches
small comoving distance, hence high k — so it needs `kmax ≈ 0.7`, and a run
tuned only for TT undershoots φφ by ~20% at ℓ = 2000. And the lensed BB is a
pure off-diagonal EE→BB transfer, so its top few hundred multipoles are
under-converged unless the unlensed spectra are carried a few hundred ℓ beyond
`lmax` (the `dl_buffer`). Both are documented at the function.

Two comparison traps are worth recording, because each cost real time and each
produces numbers that *look* like physics bugs. Low-ℓ C_ℓ takes large
contributions from k ≫ ℓ/η₀ — 35% of TT at ℓ = 20, and a factor of ~9 for tensor
EE/BB, whose radial functions tend to −j_ℓ/2 and j_ℓ′/2 at large argument instead
of decaying — so reference runs must be matched in *absolute* k_max, not in ℓ
range; CLASS and CAMB tie k_max to ℓ_max and an isolated low-ℓ_max reference run
is silently under-converged. And the tensor transfer oscillates in k with period
2π/η₀, so a grid coarser than ~40 points per cycle aliases EE and BB against each
other by several percent. The `tensor_cmb_spectra` defaults encode both lessons.

Run the suite with `] test Cosmic` — offline, no plotting.

### Errors worth knowing about

Each was caught only by comparing against CAMB or CLASS, and each is the kind that
produces a *plausible* answer rather than an obviously broken one:

- **Gauge.** The matter power spectrum is conventionally the comoving-gauge δ_m, not
  the Newtonian-gauge one. They agree deep inside the horizon and diverge badly
  outside it: using the wrong one is **72% wrong at k = 5e-4 /Mpc** while looking
  fine at k ≳ 0.01. See [`δ_matter_comoving`](src/perturbations.jl).
- **The Doppler sign** in the line-of-sight source. With it flipped, the acoustic
  *peaks* still land in the right places — they are set by geometry — but the
  troughs, which the Doppler term fills, come out up to 4× too deep. The spectrum
  looks like a CMB spectrum and is wrong.
- **The baryon sound speed.** c_s² = k_B T_b/(μ m_H c²)·(1 − ⅓ dlnT_b/dlna). Two things
  are easy to miss: the mean molecular weight μ = (1+4f_He)/(1+f_He+x_e) ≈ 0.59 before
  recombination (not 1), and the bracket equals 4/3 only while T_b ∝ 1/a — after thermal
  decoupling the gas cools as 1/a² and it becomes 5/3. Getting either wrong misplaces
  the baryon Jeans scale and costs ~1% in P(k) at k ≳ 1.
- **Hierarchy truncation.** ℓ_max = 12 for the photons is a common textbook choice and it
  leaves a **2% error in P(k) at k = 1/Mpc** — power reflected off the top of the truncated
  hierarchy leaks back down into the low multipoles that source the metric. ℓ_max_γ = 25,
  ℓ_max_ν = 32 brings P(k) to 0.04% of CAMB.
- **k-sampling.** The BAO wiggles are periodic in *linear* k with period 2π/r_s ≈ 0.04 /Mpc,
  so a log grid that resolves the turnover at k_eq gives 2–3 points per wiggle by k ~ 0.2 and
  the interpolated P(k) comes out visibly kinked. `bao_k_grid` follows CLASS: 10 points per
  decade in general, 70 inside the BAO window.

## Status

Complete and CAMB/CLASS-validated: species composition, background, recombination,
reionization, the Boltzmann hierarchy (including the momentum-resolved massive-ν
hierarchy), P(k), scalar CMB C_ℓ (TT/EE/TE), tensor C_ℓ (TT/TE/EE/BB), CMB lensing
(C_ℓ^φφ, C_ℓ^Tφ, and the lensed TT/TE/EE/BB), and gauge transformations — in flat,
open, and closed geometry, with exact hyperspherical radial functions everywhere
K ≠ 0.

Spectral distortions μ, y (and the residual r-type signal) use the **exact**
Green's-function / PCA method of Chluba & Jeong (2014): the branching between
distortion types comes from Chluba's thermalization Green's function G_th(x,z)
(shipped in `data/`, computed by CosmoTherm from the full Compton + double-Compton
+ bremsstrahlung problem), and the Silk-damping heating rate reproduces CLASS's
integrand exactly. Validated against a live CLASS 3.3.4 run
(`sd_branching_approx = exact`, FIRAS detector) on the same cosmology:

| | Cosmic | CLASS | |
|---|---|---|---|
| μ | 2.0916e-8 | 2.0921e-8 | 0.02% |
| y | 3.3753e-9 | 3.3760e-9 | 0.02% |
| k_D(z) | — | matches r_d to <0.5% | |

The near-unobservable temperature shift g is UV-sensitive and reproduced only
roughly (it needs CLASS's separate high-z WKB module); it is a monopole absorbed
into T_CMB. `distortion_spectrum` returns the full ΔI(ν); `firas_chi2` scores it
against COBE/FIRAS.

### μ and y derived from first principles

`thermalization.jl` does not read the Green's function at all — it *solves* the photon
Boltzmann equation (exact Klein–Nishina redistribution, double Compton, and
bremsstrahlung, linearised around Planck and stiff-integrated in ln(1+z)) and
projects the evolved spectrum. The Compton trace is generated with FORM; its
thermal Maxwell–Jüttner integral is evaluated in Double64 and cached only after
detailed-balance, conservation, H-theorem, node-doubling, and Kompaneets-limit
gates:

```julia
c = cosmology(m_ν = Float64[])
thermalization_distortions(c, recombination(c); sources = (nuclear = false,))
# (μ = 1.63e-8, y = 3.73e-9, Δρ_ρ = 2.77e-8)   visibility-split convention
```

(`nuclear = false` drops the ³H/⁷Be trace-decay heating so the source list matches
what CLASS integrates; it is on by default.)

The evolved spectrum is one physical Δn(x); μ and y depend on which estimator
reads it. Both, at nx = 512, against their proper references:

| same spectrum, read as | μ | y | reference | agreement |
|---|---|---|---|---|
| visibility split (Chluba 2013) | 1.630e-8 | 3.731e-9 | 1.647e-8, 3.657e-9 | 1.0%, 2.0% |
| FIRAS-band PCA (CLASS's) | 2.082e-8 | 3.332e-9 | 2.0921e-8, 3.3760e-9 | 0.5%, 1.3% |

The PCA readoff converges on CLASS as the frequency grid refines — μ/CLASS =
0.907, 0.974, 0.991, 0.995 at nx = 64, 128, 256, 512 — and swapping the exact
kernel in for the Kompaneets operator moved the converged closure by ~0.1%: the
exact treatment is physics-neutral where it should be and removes the
approximation everywhere else.

Energy is conserved to machine precision and photon number to 1e-16. The
double-Compton matrix element is independently FORM-derived and soft-factorization
gated (I₄ = ∫x⁴n_Pl(1+n_Pl)dx = 4π⁴/15 is the Bose stimulated-emission factor).

> **On the μ vs CLASS discrepancy (resolved — convention, not physics).** An earlier
> version of this file argued CLASS's `FIRAS_branching_ratios.dat` was unphysical because
> its J_μ exceeds the 1.401 photon-number ceiling. That argument was wrong: the ceiling
> only binds when the decomposition has no temperature-shift component, and CLASS's has
> one (negative), with which the file conserves energy exactly. The ~17% gap itself is
> real but is bookkeeping: μ depends on how a distortion is split between a chemical
> potential and a temperature shift. Much of ΛCDM's dissipation lands at z ≈ 5×10⁴–3×10⁵
> where the spectrum is genuinely between the μ and y shapes; the visibility split books
> that signal as y and residuals, the FIRAS-band PCA fit books more of it as μ against a
> negative ΔT. The ΔI(ν) is identical either way. This is closed by experiment, not
> assertion — reading the same evolved spectrum with CLASS's own estimator reproduces
> CLASS to 0.5%, and the check is pinned in the test suite. μ alone is not an observable:
> quote it only together with its convention.

See [`docs/src/spectral_distortions.md`](docs/src/spectral_distortions.md) for the
four bugs that each, on their own, silently zero out the answer.

## Big Bang nucleosynthesis

Y_He is not an input. `cosmology(Yp = :bbn)` solves the nuclear network for it:

```julia
b = bbn(ω_b = 0.02242, N_eff = 3.044)
b.Y_p     # 0.24560   helium mass fraction
b.DH      # 2.44e-5   deuterium (PRIMAT chain; rates = :parthenope gives 2.50e-5)
```

The n↔p weak rates are *derived* — a single Born integral over the electron momentum
(Pitrou et al., Phys. Rep. 754 (2018), eqs 77–78, 81), with the zero-temperature
Coulomb (Fermi function) and resummed Sirlin radiative corrections. G_F, V_ud and g_A
never appear: the overall constant is fixed by the measured neutron lifetime, since
that is exactly what τ_n measures.

The network is 63 reactions among 12 nuclides (through ⁶Li, ⁸Li, ⁷Be, ⁸B, ⁶He;
`network = :small` runs the 12-reaction core alone). Nuclear cross-sections *are*
measured — nobody computes ⁴He photodisintegration from QCD — so the forward rates
are a tabulated compilation: **PRIMAT throughout, by default**. One chain, one
provenance; `rates = :parthenope` swaps the 12 core reactions for PArthENoPE 3.0's
post-LUNA fits. The tables are verified as far as theory permits: the three
deuterium-burning rates are recomputed in the test suite by Gamow integration of the
published LUNA-era S-factors (they match to 0.1%, cross-checked independently in the
Wolfram engine), and all 63 reverse-rate coefficient sets are re-derived from AME
masses and ground-state spins by `detailed_balance` (they match to 0.5%, and one
transcription error in an upstream reaction label was caught exactly this way). The
tabulated 1σ envelopes are not decoration either: `bbn_mc` propagates them, redrawing
every rate lognormally plus τ_n within its error, and returns the abundance
distributions.

Against both sBBN tables CLASS ships (τ_n = 878.4 s in all three codes):

| ω_b | ΔN_eff | Cosmic | PRIMAT | diff | PArthENoPE | diff |
|---|---|---|---|---|---|---|
| 0.0224 | −1 | 0.23104 | 0.23127 | −0.10% | 0.23087 | +0.07% |
| 0.0200 | 0 | 0.24449 | 0.24457 | −0.03% | 0.24422 | +0.11% |
| 0.0224 | 0 | 0.24559 | 0.24567 | −0.04% | 0.24531 | +0.11% |
| 0.0240 | 0 | 0.24625 | 0.24633 | −0.04% | 0.24597 | +0.11% |
| 0.0224 | +1 | 0.25818 | 0.25816 | +0.01% | 0.25782 | +0.14% |

Within 0.1% of the PRIMAT table over the whole (ω_b, ΔN_eff) plane — the correct
comparison, since the default chain *is* PRIMAT's — while the two reference codes
disagree with each other by 0.15% at identical inputs.

The sharpest check is against PRyMordial run here on the same inputs with the same
compilation, which isolates the network solver itself from every data choice:

| | Cosmic (`rates=:primat`) | PRyMordial | diff |
|---|---|---|---|
| Y_p | 0.245597 | 0.245606 | −0.004% |
| D/H | 2.4375e-5 | 2.4368e-5 | +0.03% |
| ⁷Li/H | 5.49e-10 | 5.55e-10 | −1.1% |

With `rates = :parthenope` the same network gives D/H = 2.5016e-5 against PArthENoPE's
own published 2.51e-5 at these inputs (Pisanti et al 2021, eq. 16): −0.3%. The 2.6% D/H
and 13% ⁷Li/H moves between the two configurations are the genuine spread between the
compilations — mostly d(p,γ)³He and the d+d channels for deuterium, ³He(α,γ)⁷Be and the
⁷Be destruction rates for lithium (swapping just those four moves ⁷Li/H from 4.75e-10 to
5.37e-10 with Y_p and D/H untouched). Quoting any abundance without naming the rate
compilation is meaningless at this precision.

Fiducial abundances (ω_b = 0.02242, default chain) with Monte-Carlo error propagation
over all 63 rate envelopes plus τ_n (`bbn_mc`, 500 draws):

| | value | 1σ (rates + τ_n) |
|---|---|---|
| Y_p | 0.24560 | 0.00010 |
| D/H | 2.437e-5 | 0.027e-5 |
| ³He/H | 1.031e-5 | 0.013e-5 |
| ⁷Li/H | 5.4e-10 | 0.59e-10 |

(⁷Li/H is still 3× the observed Spite plateau, as it is for every standard BBN code —
that is the lithium problem, and it survives here too.)

Getting there took four corrections beyond the Born rate and one bug in the background,
each of which was worth more than the one before it looked like it should be:

| | Y_p vs CLASS |
|---|---|
| Born rate only | −1.66% |
| + Coulomb (Fermi function) + resummed Sirlin radiative | −0.34% |
| + finite nucleon mass (Fokker-Planck, Pitrou §III.G) | +0.18% |
| + non-instantaneous ν decoupling *(then a ramp; now the solved history)* | +0.013% |
| + finite-temperature plasma QED (Brown & Sawyer) | +0.020% |
| + solved decoupling history everywhere: T_νe in the weak rates, a(T_γ) as the clock, n_b/m̄ | −0.08% |
| + one consistent rate compilation (PRIMAT chain) | **−0.035%** |

The finite-temperature correction is computed, not read: `bbn_thermal.jl` does the four
2D integrals (thermal-photon exchange, differenced bremsstrahlung, thermal mass shift,
e–e scattering) with tensor-product Gauss–Legendre. PRIMAT and PRyMordial ship this as a
precomputed table and PRyMordial evaluates it with `vegas` — but a table is a *cache*, not
a limit: those files were produced by doing these integrals once. Here the cache
(`data/bbn/thermal_qed.dat`) is generated by `generate_thermal_cache()`, takes two seconds,
and deleting it just recomputes it. The integrands are smooth on a rectangular domain,
which is exactly where a deterministic rule beats Monte Carlo.

Three things about this calculation fail *silently*, and all three did:

- **N_eff is not constant during nucleosynthesis.** 3.044 is a *late-time* number: the
  excess over 3 is energy handed to the neutrinos *by* e± annihilation, because decoupling
  is not quite complete when the positrons go. Before the annihilation that transfer has
  not happened and there are exactly three species at T_ν = T_γ. Using 3.044 throughout
  overstates ρ_ν by 1.5%, hence H by 0.36%, hence — through T_freeze ∝ H^{1/3} and
  n/p = e^{−Q/T} — **Y_p by 0.18%**. It survives every other correction and looks like
  nothing.
- **Starting the network at T9 = 10.** That is where the rate tables begin, but Γ/H ≈ 1.6
  there — the weak interactions are *already freezing out*, so imposing weak equilibrium
  as an initial condition forces n/p low, and every missing neutron is two missing helium
  nucleons. It costs 7% of Y_p. The fix is the same split PRIMAT makes: run n↔p alone from
  T9 = 100 (where Γ/H ≈ 900 and equilibrium is genuinely exact) and hand over at T9 = 10.
- **Normalising a corrected rate with an uncorrected phase-space factor.** The naive argument
  says the QED corrections cancel — K = 1/(τ_n·C·λ₀) against Γ = K·C·Λ — and they would, if C
  were constant. It is not: the Coulomb factor F(b) → 2πα/b diverges for *slow* electrons,
  which is the free-decay regime, while the thermal rates at T ~ 1 MeV sit at b → 1 where
  F → 1.02. So the normalisation is enhanced more than the rates it normalises. Getting this
  wrong is worth 1.7% of Y_p, in the direction that looks like "close enough".

## N_eff, derived

`N_eff = 3.044` was the last hardcoded number in the code, and it is not a free
parameter — it is a prediction. `neff()` computes it.

It exceeds 3 for two reasons, and both are calculated rather than assumed. Neutrinos do
not decouple *quite* instantaneously, so when e± annihilate a little of their entropy
leaks into the neutrinos instead of all of it going to the photons; and the e±γ plasma is
not an ideal gas, so the amount of entropy that has to be dumped is not the free-field one.

Neither reference code computes this. PRIMAT reads the neutrino history from precomputed
NEVO tables; PRyMordial ships NUDEC_BSM's `nue_ann.txt`/`nue_scatt.txt` and interpolates.
Here every ingredient is derived, and each is separately checkable:

| ingredient | how | check |
|---|---|---|
| QED plasma pressure | closed form: Frenkel–Galitskii–Migdal at O(e²), Blaizot ring at O(e³) | **0.03%** vs PRIMAT's table, 0.005–40 MeV |
| \|M\|² for e⁺e⁻ → νν̄ | FeynCalc | — |
| \|M\|² for νe → νe | FeynCalc | mₑ² term flips sign vs annihilation (crossing) |
| 2-body massless phase space | analytic tensor reduction | reproduces its own trace |
| F_ann(T) | phase-space integral of the matrix element | **<1%** vs NUDEC_BSM for T ≥ 1 MeV |
| F_scat(T) | full 2→2 collision integral with CM boost | detailed balance = 3.8e-17; **0.1–0.4%** vs NUDEC_BSM |

The O(e³) term deserves a note: it is *half-integer* order in α, not O(e⁴), because the
soft-photon loop is infrared divergent order by order and has to be Debye-resummed. A
straight perturbative expansion in α never produces it.

The check worth trusting is the flavour split. g_L has opposite sign for ν_e and ν_μ, so
the g_L g_R mₑ² term enters with opposite sign, and the finite-mass corrections come out
genuinely *different* (F_ann = 0.988 vs 0.953 at 1 MeV). Nothing tells the code to do that
— it falls out of the matrix element. A sign error there would still give plausible
numbers, just identical ones.

Every load-bearing formula was re-derived independently (Mathematica, FORM, or
symbolic algebra in Julia) rather than trusted from a single source.

## Not done yet

Being precise about this, because "does what CLASS and CAMB do combined" is not yet true:

- **Vector modes** — absent.
- **Isocurvature initial conditions** — adiabatic only.
- **Halofit / nonlinear P(k)** — absent; the matter power spectrum is linear.
- The dominant late-time (reionization + structure formation) contribution to the y-distortion.
- **Neutrino flavour oscillations in decoupling** — the derived N_eff = 3.0428 is the
  no-oscillation value; the full quantum kinetic treatment that moves it to 3.044 is not done.

## Contributing

Open an issue or pull request on [GitHub](https://github.com/aburousan/Cosmic.jl).
New features should come with tests and a bumped minor version.
