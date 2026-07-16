```@meta
CurrentModule = Cosmic
```

# Primordial spectra and inflation

## The primordial hook

The Boltzmann solver runs at unit primordial curvature, so the primordial
spectrum multiplies transfer functions linearly and can be swapped without
re-solving anything. Three levels of generality:

The default power law ``\Delta^2_\mathcal{R}(k) = A_s (k/k_\star)^{n_s-1}``
from the `A_s`, `n_s`, `k_pivot` parameters.

An arbitrary spectrum: `cosmology(primordial = f)` accepts any callable
``k \mapsto \Delta^2_\mathcal{R}(k)`` — features, steps, running, broken scale
invariance, a table — and every observable (CMB, lensing, ``P(k)``, induced
gravitational waves) picks it up through the same hook.

Correlated multi-mode spectra: `PrimordialMatrix` specifies auto- and
cross-power for any set of adiabatic and isocurvature modes,
``\mathcal{P}_{ij}(k) = A_{ij} (k/k_\star)^{n_{ij}-1}``, and the corresponding
`cmb_spectra` method evolves each mode's transfer functions and sums every
ordered pair — fully correlated isocurvature phenomenology, not just pure
modes.

## Inflation

`inflaton_spectrum` computes the spectrum instead of positing it. The
background trajectory of a canonical field in a user potential ``V(\phi)``
(reduced Planck units) is integrated exactly to the end of inflation, and each
comoving mode's Mukhanov–Sasaki equation is solved from Bunch–Davies initial
conditions deep inside the horizon through freeze-out — no slow-roll
truncation anywhere, so any feature the potential imprints appears in
``\Delta^2_\mathcal{R}(k)`` with its true shape. Tensors are evolved the same
way for ``\Delta^2_t(k)`` and ``r``. For the workhorse ``m^2\phi^2`` test the
exact evolution lands on the slow-roll attractor predictions to 0.03% in
``n_s``.

An effective sound speed ``c_s(N)`` (the EFT-of-inflation quadratic action,
``z = a\sqrt{2\varepsilon}/c_s``) can be supplied for models with non-trivial
kinetic structure. The solver refuses ``c_s > 1`` — superluminal propagation
is a causality violation, whatever the literature it appears in — and warns
when ``c_s`` drops low enough to raise strong-coupling concerns. Sharp
transitions in ``c_s`` or the potential are handled by the exact integration,
which makes the code suited to studying how transition *width* affects the
resulting features: observables from "instantaneous" transitions can depend
strongly on the smoothing assumed, and here that sensitivity can be quantified
rather than assumed away.

The returned `InflatonSpectrum` is callable and plugs directly into
`cosmology(primordial = s)` and into the induced-gravitational-wave pipeline
(see [Scalar-induced gravitational waves](sigw.md)).

## Reference

```@autodocs
Modules = [Cosmic]
Pages   = ["inflation.jl"]
Private = false
```
