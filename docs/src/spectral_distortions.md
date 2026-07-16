# CMB spectral distortions

The CMB is a blackbody to about one part in 10⁴. It is not a *perfect* one, and the
deviations carry information no anisotropy can: they record how much energy was
injected into the photon bath, and *when*.

The story is set by one competition -- how fast photons can be created or shuffled
in energy, versus how fast the universe is expanding.

## The three eras

**Thermalization (z ≳ 2×10⁶).** Double Compton scattering (γ + e → 2γ + e) and
bremsstrahlung create and destroy photons fast enough to restore a perfect
blackbody after any energy release. Distortions here are erased; the universe has
no memory of them. This is why spectral distortions probe no earlier than
z ~ 2×10⁶ -- a hard wall.

**The μ era (5×10⁴ ≲ z ≲ 2×10⁶).** Compton scattering (which conserves photon
number) still redistributes photons in energy efficiently, but photon-*number*-
changing processes have frozen out. Energy injected here brings the photons to a
Bose-Einstein spectrum with a non-zero chemical potential μ -- a *kinetic*
equilibrium that is not a blackbody. The signature is a μ-distortion.

**The y era (z ≲ 5×10⁴).** Compton scattering can no longer fully redistribute
energy. Injected energy leaves the low-frequency photons under-heated and the
high-frequency ones over-heated relative to a blackbody: a y-distortion, the same
spectral shape as the thermal Sunyaev-Zel'dovich effect from galaxy clusters.

Between the two eras is a continuous transition (the "i-type" or residual
distortion), which is what actually distinguishes energy-release *histories* rather
than just totals.

## Green's function method (Chluba 2013)

Because the Boltzmann collision term is linear in small distortions, the distortion
today from an energy-release history dQ/dz is a convolution:

    ΔI(ν) = ∫ dz  G_th(ν, z)  (dQ/dz)/(ργ)

G_th is the *thermalization Green's function*: the spectrum today produced by a
delta-function energy release at redshift z. Chluba computed it once, numerically,
from the full two-photon/Compton/bremsstrahlung problem, and it is what
`Greens.v1.0.3` tabulates. Given G_th, any heating history is a single integral.

A useful and accurate shortcut collapses G_th onto just two numbers via visibility
weights (Chluba 2016):

    μ ≈ 1.401 ∫ dz  𝒥_μ(z) (dQ/dz)/ργ,     𝒥_μ(z) = exp[-(z/z_μ)^{5/2}]
    y ≈ (1/4) ∫ dz  𝒥_y(z) (dQ/dz)/ργ,     𝒥_y(z) = 1 - 𝒥_μ(z)

with z_μ ≈ 1.98×10⁶ (the thermalization redshift). The 1.401 is the standard
μ-to-energy conversion (μ = 1.401 ΔQ/Q_γ in the μ era), and the 1/4 likewise for y.

## Sources of energy release

- **Silk damping of small-scale acoustic modes.** As perturbations on scales that
  enter the horizon before recombination oscillate, photon diffusion erases them,
  and the dissipated energy heats the photon bath. This is a *guaranteed* ΛCDM
  signal, fixed by the primordial power spectrum:

      (dQ/dz)/ργ ≈ ∫ (dk/k) Δ²_ℛ(k) · [d/dz of the damping envelope]

  It probes P(k) at k ~ 10²–10⁴ /Mpc, far smaller scales than the CMB or galaxy
  surveys reach, so a μ measurement is a constraint on inflation.

- **Adiabatic cooling of baryons.** The baryons cool as (1+z)² and draw heat *out*
  of the photons via Compton coupling -- a small negative distortion, the one
  guaranteed signal of the opposite sign.

- **Particle decay / annihilation, primordial black holes, dark matter.** Any
  exotic energy injection with the right timing shows up here; this is why
  distortions are a clean probe of the dark sector.

## Current limits and forecasts

FIRAS (COBE, 1996): |μ| < 9×10⁻⁵, |y| < 1.5×10⁻⁵ (95%). The ΛCDM prediction is
μ ~ 2×10⁻⁸ and y ~ 4×10⁻⁹ -- three to four orders of magnitude below FIRAS, and the
target of proposed spectrometers (PIXIE, and the data files `PIXIE_*.dat` here).

## Deriving μ and y from scratch (the thermalization solver)

Reading Chluba's precomputed `G_th` table is exact but not *derived*. `thermalization.jl`
instead solves the photon Boltzmann equation directly, so μ and y come out of the
physics rather than a file:

    dn(ω₀)/dτ = ∫dω [(ω²/ω₀²)P(ω→ω₀)n(ω)(1+n₀)
                          - P(ω₀→ω)n₀(1+n(ω))]              exact Compton
                + (1/x³) [ 1 - n(e^{x_e} - 1) ] · K(x)          double Compton + bremsstrahlung

with τ the Compton optical depth (dτ = σ_T nₑ c dt) and x = hν/kT_ref, T_ref = T₀(1+z).
Here P is the exact Klein–Nishina kernel folded with the relativistic
Maxwell–Jüttner electron distribution. Its tree trace is generated with FORM and
its cancellation-prone thermal integral with Double64; the installed binary is a
cache of those verified calculations, not external physics data. Pure expansion
leaves x fixed, so only the collision terms move photons. The equation is
linearised around Planck (the distortion is ~1e-8 of the blackbody, so evolving
the full n(x) and subtracting loses every digit) and evolved with a stiff solver.

### Emission coefficients (Chluba & Sunyaev 2012, §2.2.2)

    K_DC = (4α/3π) θ² · g_dc(x,θ)
    I₄^pl   = ∫ x⁴ n_Pl(n_Pl+1) dx = 4π⁴/15 ≈ 25.976

    K_BR = [ α λ_e³ θ^{-7/2} e^{-x} / (2π√(6π)) ] · Σᵢ Zᵢ² Nᵢ g_ff,   λ_e = h/(m_e c)

I₄^pl is the **seed-photon phase space**: double Compton emission is *stimulated* by
the ambient blackbody, so the rate carries a factor of ~26 that the headline formula
K_DC = (4α/3π)θ² does not advertise. `g_dc(x,θ)` is evaluated from the
FORM-derived 2→3 matrix element and cached after positivity, soft-factorization,
and Double64 convergence gates; the common Chluba–Sunyaev polynomial is retained
only as a comparison. Likewise `g_ff` is the thermally averaged exact
Karzas–Latter Coulomb-wave result, not a Born logarithm or unity.

### Four things that each, on their own, destroy the answer

None of them announces itself; each returns a plausible-looking number (or exactly
zero) rather than an error.

1. **The discrete kernel must obey the H theorem.** A pointwise interpolation of
   the redistribution integral can conserve neither photon number nor entropy;
   at low temperature its interpolation error can even reverse the thermal
   energy flow or create growing modes. The production operator uses the
   detailed-balance finite-element weak form. Every exact transition contributes
   a negative-semidefinite rank-one pair exchange at its interpolated sub-cell
   destination. Photon number, the Bose–Einstein null mode, common-temperature
   equilibrium, and non-positive entropy production then hold algebraically;
   grid doubling checks the remaining frequency discretization error.

2. **Compton conserves photon *number*, not *energy*.** It trades energy with the
   electrons. Holding T_e = T_ref dumps the distortion's energy into an infinite
   electron bath -- 100% of it was destroyed at z = 1e6. Electrons have negligible
   heat capacity, so T_e is solved from the exact discrete energy ledger. Compton
   redistribution and DC/BR emission share this same ledger; assigning the
   injected energy independently to both would double count it.

   **This term is what lets the injected y-distortion Comptonize into a Bose-Einstein
   μ. Without it, no μ forms at all.**

3. **Detailed balance must survive tabulation.** Forward and reverse cache entries
   obey the exact thermal relation, Planck is stationary, and a constant
   chemical-potential perturbation is a null mode. These are tested after table
   interpolation, not only on the analytic kernel before caching.

4. **x_c, and where the grid may end.** The linearised emission is
   -K(x)(e^x-1)/x³ · Δn, which goes as **1/x², not 1/x³**, so the critical frequency is
   x_c = √(Λ/κ) ≈ 4e-3 - 9e-3 across the μ era. The grid must end below x_c but not far
   below: the 1/x² stiffness breaks the solver by xmin = 1e-4, while the answer is
   already converged at xmin = 1e-3.

### Validation, and which reference is the right one

Green's function (delta-function energy injection at z, read off μ and y). The
coefficient of M(ν) that the projection returns should equal 1.401·J_μ(z)·𝒥(z)
(Chluba 2013 eq. 6):

| z | derived | **Chluba 2013** | CLASS's FIRAS file |
|---|---|---|---|
| 10⁶     | 1.173 | **1.168** | 1.171 |
| 3×10⁵   | 1.414 | **1.388** | 1.433 |
| 10⁵     | 1.354 | **1.312** | *1.896* |
| 5×10⁴   | 0.713 | **0.744** | *1.502* |

Totals for ΛCDM, convolving each branching with the same heating history:

|            | μ         | y         |
|------------|-----------|-----------|
| **derived (this solver, exact kernel, nx = 512)** | **1.630e-8** | **3.731e-9** |
| **Chluba 2013 (published)** | 1.647e-8 | 3.657e-9 |
| CLASS      | 1.941e-8  | 3.467e-9  |

**The solver reproduces Chluba's own published Green's function to 1.0% (μ) and
2.0% (y).** (The earlier Kompaneets-operator version of the solver gave 1.606e-8;
replacing it with the exact Klein–Nishina ⊗ Maxwell–Jüttner kernel moved the
converged answer by well under a percent — the approximation it removed was
already good in this regime, which is exactly what the swap was supposed to show.)

#### The CLASS discrepancy, and a wrong argument I made about it

The derived μ = 1.630e-8 reproduces Chluba's own published Green's function (1.647e-8) to
1.0%, but CLASS gives 1.941e-8 — a ~17% gap. What follows is first the wrong argument,
then the actual resolution.

An earlier version of this document claimed the gap was an artefact: that CLASS's
`FIRAS_branching_ratios.dat` was a PCA decomposition in a detector-specific basis rather
than physical branching ratios, and that chasing the difference was chasing a phantom. The
evidence offered was that its J_μ reaches 1.93 at z ≈ 10⁵, above the photon-number ceiling

    ΔN/N = 3.0 Δ_T - 1.3686 μ = 0   ⟹   Δ_T = 0.4562 μ
    ΔQ   = 4 Δ_T - 1.1106 μ = 0.7142 μ   ⟹   μ = 1.401 ΔQ

**That argument is wrong, and the claim is withdrawn.** The ceiling μ ≤ 1.401 ΔQ only binds
when the distortion has *no temperature-shift component*. CLASS's decomposition has one, and
it is negative (J_T = −0.13 at the peak). With it included the file conserves energy exactly:

    4·J_T + 4·J_y + J_μ/1.401 ≈ 1.00   at every z

so the branching ratios are physical and self-consistent, and a J_μ above 1.401 is not
evidence of anything.

The remaining ~17% is now **resolved, and it is a convention, not an error** — established
by running CLASS 3.3.4 itself (`sd_branching_approx = exact`, FIRAS detector) on this exact
cosmology and comparing term by term:

|  | μ | y |
|---|---|---|
| `distortion_amplitudes` (this code, FIRAS PCA file) | 2.0916e-8 | 3.3753e-9 |
| CLASS 3.3.4 live run | 2.0921e-8 | 3.3760e-9 |
| agreement | **0.02%** | **0.02%** |

All six PCA residual amplitudes μ_k agree to ~0.3% as well. (A bookkeeping note that cost
a detour: CLASS divides the file's J_μ by 1.401 when reading and multiplies by 1.401 after
convolving — likewise a factor 4 for J_T and J_y. The factors are baked into the file
itself, so convolving it as-is, as this code does, lands on the same answer by a shorter
path. The file here is column-for-column identical to CLASS's.)

So the two decompositions coexist, both implemented and both validated against their own
reference: `distortion_amplitudes` uses the Chluba & Jeong (2014) PCA split in the FIRAS
band and reproduces CLASS to 0.02%; `thermalization_distortions` solves the photon
Boltzmann equation and reproduces Chluba (2013)'s visibility split to 2.5%.

That the 17% between them is bookkeeping and not physics is not asserted — it is
**closed by experiment** (and pinned in the test suite): take the thermalization
solver's evolved spectrum, one physical Δn(x), and read μ off it with both estimators.

| estimator, same spectrum (exact kernel, nx = 512) | μ | y |
|---|---|---|
| grid projection onto {G, Y, M} (visibility style) | 1.630e-8 | 3.731e-9 |
| FIRAS-band weighted fit with {G, Y, M, S₁..₆} (CLASS's) | 2.082e-8 | 3.332e-9 |
| CLASS's own pipeline, for reference | 2.092e-8 | 3.376e-9 |

The CLASS-style estimator applied to our own spectrum lands within 0.5% of CLASS
(converging as the frequency grid refines: μ/CLASS = 0.907, 0.974, 0.991, 0.995
at nx = 64, 128, 256, 512) — and its fitted temperature shift comes out *negative*,
which is precisely the negative J_T in CLASS's branching file. The mechanism: much of
ΛCDM's dissipation happens at z ≈ 5×10⁴–3×10⁵, where the spectrum is genuinely between
the μ and y shapes; the FIRAS-band fit books more of that intermediate signal as μ
(compensated by ΔT < 0), the visibility split books it as y and r-type residuals. The
observable ΔI(ν) is identical either way. μ alone is not an observable — quote it only
together with its convention.

## Design note: modified gravity

Every layer reaches the expansion history only through `H_Mpc`, `E`, `ℋ`, and the
species densities, and the thermal/distortion layers reach it only through those
plus `n_e` and `R_baryon`. Nothing hard-codes the Friedmann equation or GR
Poisson law inline. A modified-gravity model therefore enters at two well-defined
seams: (i) a species (or a modified `E`) for the background, and (ii) the
metric-potential closure in `perturbations.jl` (the μ_MG, Σ_MG functions relating
Φ, Ψ to the matter). The thermalization and background-derived quantities recover
their ΛCDM values automatically once those seams are set to GR -- a property to be
tested when the MG layer is added.
