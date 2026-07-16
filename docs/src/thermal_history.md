```@meta
CurrentModule = Cosmic
```

# Thermal history: recombination and reionization

## Recombination

The ionization history is solved with a port of HyRec's effective multi-level
atom: hydrogen recombination through the exact effective rates (with the
two-photon 2s→1s channel and Lyman-α escape), helium recombination with the
Sobolev treatment, and the matter temperature evolved through Compton coupling
to the radiation. There is no RECFAST-style fudge-function calibration; where
effective rates are tabulated they are the published HyRec tables, shipped
with provenance in `data/hyrec/`.

`recombination(c)` returns a `RecombinationSolution` from which the standard
derived quantities follow: `x_e`, `T_matter`, the Thomson `optical_depth` and
`visibility`, the baryon-drag depth, `z_star`, `z_drag`, the sound horizon
`r_star` and `r_drag`, and the acoustic scale `θ_star`. The Saha functions are
available separately for cross-checks and initial conditions.

## Reionization

The standard tanh-in-``(1+z)^{3/2}`` parameterization with helium double
reionization, matching the CLASS convention. The default `z_reio = 7.67`
reproduces the Planck 2018 ``\tau = 0.0544``; `z_reio_from_τ` inverts the map
when the optical depth is the preferred parameter. Pass `no_reionization` to
switch it off.

## Exotic energy injection

`Injection` adds dark-matter annihilation (`p_ann`, the CLASS
parameterization ``\langle\sigma v\rangle/m_{\rm DM}`` with the injected power
``\rho_{\rm cdm}^2 c^2 p_{\rm ann}``, optionally halo-boosted at late times)
and decay (`Γ_dec`, `f_dec`) into the pre-reionization gas. The extra
ionization and heating widen the visibility function and raise the residual
``x_e`` at ``z \sim 100``–600, which is precisely how the CMB constrains these
processes. The trace nuclear energy release computed by the BBN network can
also be propagated (`trace_nuclear_history`, `nuclear_deposition_rate`),
closing the loop between nucleosynthesis and recombination physics.

## Reference

```@autodocs
Modules = [Cosmic]
Pages   = ["thermodynamics.jl", "hyrec.jl", "nuclear_deposition.jl"]
Private = false
```
