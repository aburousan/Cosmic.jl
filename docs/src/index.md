```@meta
CurrentModule = Cosmic
```

# Cosmic.jl

[Cosmic.jl](https://github.com/aburousan/Cosmic.jl) is a cosmology code in pure
Julia: a Boltzmann solver in the tradition of CLASS and CAMB, together with the
early-universe machinery those codes outsource — Big Bang nucleosynthesis,
neutrino decoupling, CMB spectral distortions — solved within the same package
and the same conventions.

The guiding principle is exactness. Wherever a standard code would switch to a
fitting formula or a truncated expansion, Cosmic.jl either integrates the exact
expression or documents precisely what is cached and how it was generated.
Every physics module is cross-checked against an independent reference
(CLASS, CAMB, FAST-PT, or published closed forms) before it is trusted; the
numbers quoted below come from those comparisons.

## What it computes

The homogeneous universe: any mixture of photons, baryons, cold dark matter,
massless and massive neutrinos (with chemical potentials or arbitrary
distribution functions), warm dark matter, decaying dark matter, spatial
curvature, and dark energy as a cosmological constant, a ``w_0``–``w_a`` fluid,
an arbitrary ``w(a)``, or a quintessence field rolling in a user-supplied
potential. Distances, ages, horizons, and the BBN abundances follow from the
one `cosmology(...)` object.

Recombination through a port of HyRec (effective multi-level atom), with
reionization, exotic energy injection from annihilating or decaying dark
matter, and the derived quantities of CMB analysis: optical depths,
visibility, the sound horizon, ``z_*``, ``z_{\rm drag}``, ``\theta_*``.

The linear perturbations: a stiff-integrated Boltzmann hierarchy in the
conformal Newtonian gauge (no tight-coupling approximation) with
momentum-resolved massive neutrinos, adiabatic and all four regular
isocurvature initial conditions, decaying-dark-matter and dark-energy
perturbation sectors, and tensor modes — in flat, open, and closed geometries.
CMB temperature and polarization agree with CLASS at the 0.3% level; the
lensed spectra (computed by the exact correlation-function method) at 0.3–0.4%
with ``C_\ell^{\phi\phi}`` to better than 1%.

The matter power spectrum, linear and nonlinear: halofit (0.25% against
CLASS), HMcode-2020 (sub-percent against CAMB to ``k = 1\,h``/Mpc), and the
exact 1-loop EFT of large-scale structure with IR resummation, whose loop
integrals match the independent FFTLog code FAST-PT to better than 0.5%.

The primordial universe: power-law, arbitrary, or multi-field-correlated
primordial spectra; an inflaton module that integrates the exact
Mukhanov–Sasaki equation for any potential ``V(\phi)`` (with an optional
effective sound speed ``c_s(N)``); and the scalar-induced gravitational-wave
background from the Kohri–Terada kernel, validated against its published
closed forms.

CMB spectral distortions: the thermalization Green's function with exact
double-Compton and Karzas–Latter free–free coefficients, and ``\mu``, ``y``
amplitudes from any injection history.

## Installation

Cosmic requires Julia 1.9+:

```
pkg> add https://github.com/aburousan/Cosmic.jl
```

## Quick start

```julia
using Cosmic

c   = cosmology()                          # Planck-2018-like defaults, Yp from BBN
rec = recombination(c)                     # HyRec thermal history
cl  = cmb_spectra(c, rec; lmax = 2500)     # unlensed TT/TE/EE
len = lensed_cmb_spectra(c, rec)           # lensed spectra + Cℓ^φφ
P   = matter_power_spectrum(c, rec; z = 0) # linear P(k)
hm  = hmcode_power(P)                      # nonlinear P(k)
```

A feature spectrum from an inflaton potential, and its induced gravitational
waves:

```julia
V(φ)  = 0.5 * 6e-6^2 * φ^2
dV(φ) = 6e-6^2 * φ
s  = inflaton_spectrum(V, dV; φ0 = 16.0, N_pivot = 55)
c  = cosmology(primordial = s)             # the computed spectrum feeds everything
gw = sigw_spectrum(s, c; kmin = 1e4, kmax = 1e9)
sigw_ΔNeff(gw, c)                          # the integrated N_eff of the background
```

Examples live in the
[Examples](https://github.com/aburousan/Cosmic.jl/tree/main/Examples) folder.

## Documentation map

- [Background cosmology and species](background.md)
- [Neutrinos and relativistic relics](neutrinos.md)
- [Thermal history: recombination and reionization](thermal_history.md)
- [Big Bang nucleosynthesis](bbn.md)
- [Linear perturbations](perturbations.md) and the
  [derivations of the perturbation equations](perturbation_equations.md)
- [CMB spectra, lensing, and tensors](cmb.md)
- [The matter power spectrum, linear and nonlinear](matter_power.md)
- [Primordial spectra and inflation](primordial.md)
- [Scalar-induced gravitational waves](sigw.md)
- [Spectral distortions](spectral_distortions.md)
- [API index](api.md)

## Contributing

Open an issue or pull request on the
[GitHub page](https://github.com/aburousan/Cosmic.jl). New features should come
with unit tests, and — in the spirit of the package — with a cross-check
against an independent reference wherever one exists.
