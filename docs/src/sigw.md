```@meta
CurrentModule = Cosmic
```

# Scalar-induced gravitational waves

Scalar and tensor perturbations decouple at first order, but not at second:
curvature perturbations source gravitational waves through the quadratic terms
of the Einstein equations. For a nearly scale-invariant spectrum the effect is
negligible, but any model that enhances small-scale power — the inflationary
feature scenarios that produce primordial black holes above all — radiates a
stochastic background that pulsar-timing arrays and space interferometers can
reach. The induced background is the unavoidable counterpart of PBH
formation, which makes it both a discovery channel and a constraint.

## Implementation

The calculation is the semianalytic Kohri–Terada one: gravitational waves
generated during radiation domination, the source integrated analytically
against the exact Green's function, leaving a two-dimensional integral of a
universal kernel against the primordial spectrum,

```math
\Omega_{\rm GW,RD}(k) \;=\; \frac{1}{12}\int_0^\infty \!dt \int_{-1}^{1}\!ds\;
\mathcal{K}(t,s)\; \mathcal{P}_\zeta(k u)\,\mathcal{P}_\zeta(k v),
```

with the kernel carrying the resonance structure (a logarithmic amplification
where the source oscillates at twice the potential's frequency) evaluated in
closed form. `Ω_gw_rd` performs this integral for any spectrum callable —
including an `InflatonSpectrum`, so the chain from a potential ``V(\phi)`` to
an observable ``\Omega_{\rm GW}(f)`` has no hand-offs. The implementation
reproduces the published anchors of the method exactly: the scale-invariant
result ``\Omega_{\rm GW} = 0.8222\,A_\zeta^2``, the power-law coefficients
``Q(n_s)``, and the monochromatic closed form (`sigw_monochromatic`), against
which the numerical kernel converges in the delta-function limit.

`sigw_spectrum` maps the radiation-era result to today,
``\Omega_{\rm GW,0} = \Omega_{\rm GW,RD}\,\Omega_{r,0}`` with an optional
``g_*`` correction hook, and returns frequencies in Hz
(``f = 1.55\times10^{-15}\,k\,\mathrm{Mpc}``).

## The integrated bound

A gravitational-wave background is radiation, and its integral over frequency
counts toward ``N_{\rm eff}``:

```math
\Delta N_{\rm eff} = \frac{8}{7}\left(\frac{11}{4}\right)^{4/3}
\frac{1}{\Omega_{\gamma,0}} \int \Omega_{\rm GW,0}\, d\ln k .
```

`sigw_ΔNeff` evaluates this, and because the BBN network lives in the same
package, the bound can be applied self-consistently: feed the result back
through `cosmology(N_eff = 3.044 + ΔN)` and check the helium and deuterium
abundances directly, rather than quoting a detached limit. Spectra with large
peaks that clear a detector's sensitivity curve but fail this integral
constraint are not viable, and checking is one function call.

## Reference

```@autodocs
Modules = [Cosmic]
Pages   = ["sigw.jl"]
Private = false
```
