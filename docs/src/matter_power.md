```@meta
CurrentModule = Cosmic
```

# The matter power spectrum

## Linear

`matter_power_spectrum` runs the Boltzmann hierarchy over a ``k``-grid (with
sampling boosted through the BAO range) and stores the matter transfer
function; `power` evaluates ``P(k)``, `transfer` the transfer function itself,
and `σ_R`/`σ8` the variance in spheres. The gauge subtlety matters at large
scales: `δ_matter` is the Newtonian-gauge density contrast while
`δ_matter_comoving` is the comoving (CAMB-convention) one that behaves like
the textbook growth function outside the horizon.

## Nonlinear

Three routes, in increasing rigour about their own regime of validity:

Halofit (Takahashi calibration with the Bird neutrino correction) as the
fast baseline — the coefficients are a published fit to N-body simulations,
and are labelled as such. Against CLASS's implementation: 0.25%.

HMcode-2020, the current standard for weak-lensing accuracy: spherical-collapse
``\delta_c`` and ``\Delta_v``, a Sheth–Tormen mass function, NFW profiles with
the bloating and two-halo damping parameters of Mead et al., and the BAO
de-wiggling treatment. Sub-percent against CAMB's `mead2020` to
``k = 1\,h/\mathrm{Mpc}``; the remaining few-percent residual at
``k \gtrsim 3\,h/\mathrm{Mpc}`` traces to the halo-concentration formation
history and is being tightened.

The effective field theory of large-scale structure for the mildly nonlinear
regime (``k \lesssim 0.3\,h/\mathrm{Mpc}``) where it is exact perturbation
theory rather than a fit: the standard 1-loop integrals ``P_{22}`` and
``P_{13}`` evaluated by direct quadrature over the linear spectrum (in
variables chosen so the mode-coupling ridge is resolved), the effective
sound-speed counterterm ``c_s^2`` exposed as the free parameter it is, and IR
resummation of the BAO wiggle. The loop integrals match FAST-PT — an entirely
independent FFTLog implementation — to better than 0.5% across
``0.01 < k < 0.5\,h/\mathrm{Mpc}``, and each kernel is verified against its
known analytic limits.

## Reference

```@autodocs
Modules = [Cosmic]
Pages   = ["observables.jl", "nonlinear.jl", "hmcode.jl", "eft.jl"]
Private = false
```
