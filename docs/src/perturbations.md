```@meta
CurrentModule = Cosmic
```

# Linear perturbations

`solve_perturbations` integrates the full scalar Boltzmann hierarchy for one
Fourier mode: photons (temperature and polarization), massless neutrinos,
baryons, cold dark matter, momentum-resolved massive neutrinos, and — when the
cosmology contains them — decaying dark matter with its dark-radiation
product, and dark-energy fluid perturbations. The integration is in the
conformal Newtonian gauge with the metric potential ``\phi`` evolved from the
energy constraint (numerically stable outside the horizon, where the
algebraic alternative cancels catastrophically) and ``\psi`` obtained from the
anisotropic stress.

Two design choices differ from the standard codes. There is no tight-coupling
approximation: the pre-recombination era, where the Thomson rate exceeds every
other scale in the problem, is handled by a stiff implicit integrator instead
of an analytic expansion, so there is no switching error to tune. And the
hierarchy truncations were set by convergence against CAMB rather than
tradition — the defaults ``\ell_{\max,\gamma} = 25``, ``\ell_{\max,\nu} = 32``
bring ``P(k)`` to 0.04% between ``k = 0.1`` and ``2\,\mathrm{Mpc}^{-1}``,
where the textbook ``\ell_{\max} = 12`` leaves a 2% error from power reflected
off the top of the tower.

The derivations of the evolution equations, including the curved-space
generalizations, are collected in
[Perturbation equations](perturbation_equations.md).

## Initial conditions

Each mode starts when it is both far outside the horizon and deep in radiation
domination — a per-mode start time, since a fixed early start would be wasted
effort at large scales and wrong at small ones. The adiabatic growing mode is
normalized to unit primordial curvature ``\mathcal{R} = 1``; the system being
linear, the primordial spectrum multiplies afterwards.

The four regular isocurvature modes are implemented alongside: cold-dark-matter
(`:cdi`), baryon (`:bi`), neutrino-density (`:nid`) and neutrino-velocity
(`:niv`) isocurvature, each normalized to unit initial entropy. The
super-horizon solutions are the Bucher–Moodley–Turok expansions written in the
synchronous gauge and rotated to Newtonian through the Einstein constraints.
Pure-mode CMB spectra reproduce the characteristic isocurvature acoustic-peak
positions of CLASS (first CDI peak at ``\ell \approx 300`` where the adiabatic
peak sits at 220), and correlated mixtures are handled by `PrimordialMatrix`
(see [Primordial spectra](primordial.md)).

## Dark-energy perturbations

Any non-``\Lambda`` dark energy is perturbed automatically. The fluid
equations with a rest-frame sound speed ``\hat c_s^2`` (default 1) are evolved
in the variables ``(\delta, Q \equiv (1+w)\theta)``, which are algebraically
identical to the textbook pair away from ``w = -1`` but remain regular there —
necessary for quintessence, which starts frozen at exactly ``w = -1``. For a
canonical scalar field these fluid equations with ``\hat c_s^2 = 1`` are the
exact field perturbations, no Klein–Gordon block required.

Parameterized fluids that genuinely cross ``w = -1`` (the region current
baryon-acoustic-oscillation data mildly favour) leave the single-fluid
description ill-defined at the crossing; there the PPF construction of Fang,
Hu and Lewis takes over automatically (`de_scheme = :auto`), evolving one
regular variable ``\Gamma`` matched to scalar-field behaviour on both sides.
The crossing response of the CMB agrees with CLASS's PPF implementation.

## Gauges

The solution is stored in the Newtonian gauge, and `gauge_shift`, `δ_species`,
`θ_species` translate to synchronous, comoving, uniform-density and
spatially-flat gauges; the gauge-invariant curvature perturbations
`curvature_ℛ`, `curvature_ζ` and the Bardeen potentials are provided directly.
``\mathcal{R}`` conservation outside the horizon is part of the test suite.

## Reference

```@autodocs
Modules = [Cosmic]
Pages   = ["perturbations.jl", "gauges.jl"]
Private = false
```
