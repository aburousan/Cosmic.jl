```@meta
CurrentModule = Cosmic
```

# Neutrinos and relativistic relics

## The effective number of neutrino species

``N_{\rm eff} = 3.044`` is usually quoted; here it can also be computed.
`neff()` solves the neutrino decoupling era directly — the coupled evolution of
the neutrino and electromagnetic temperatures through ``e^+e^-`` annihilation,
including the finite-temperature QED corrections to the photon–electron
equation of state (`qed_pressure`, `qed_energy_density`) and the
energy-transfer integrals `F_ann` and `F_scat`. The no-oscillation value this
yields, ``N_{\rm eff} = 3.0428``, is the honest output of the physics
implemented; the remaining ``+0.001`` from flavour oscillations requires the
full quantum-kinetic treatment. Passing `N_eff = :derived` to `cosmology()`
uses the computed value instead of the standard 3.044.

## Distribution functions

The momentum distribution ``f_0(q)`` a neutrino is born with enters the physics
in exactly one place, the `NeutrinoDistribution` abstraction, and everything
downstream — background density and pressure, and the momentum-resolved
free-streaming hierarchy — is built from it.

`FermiDirac(ξ)` carries a degeneracy parameter ``\xi = \mu/T`` (a lepton
asymmetry). A nonzero ``\xi`` enhances the relativistic energy density by the
exact factor

```math
\frac{\rho_\nu(\xi)}{\rho_\nu(0)} = 1 + \frac{30}{7}\left(\frac{\xi}{\pi}\right)^2
                                     + \frac{15}{7}\left(\frac{\xi}{\pi}\right)^4,
```

which is what `degenerate_Neff_factor` returns and what a degenerate species
contributes as extra ``\Delta N_{\rm eff}``. Set per-species asymmetries with
`cosmology(ξ_ν = [ξ1, ξ2, ξ3])` or one common value with `ξ_ν = 0.1`.

`TabulatedNu(f0, dlnf0)` accepts an arbitrary distribution — non-thermal
relics, sterile admixtures — as callables for ``f_0(q)`` and
``d\ln f_0/d\ln q``.

## Massive neutrinos in the perturbations

A massive free-streaming species cannot be described by a fluid: particles of
different momenta move at different speeds through the same point, and that
velocity dispersion is what erases small-scale structure. The Boltzmann
hierarchy is therefore solved per momentum bin (`MassiveNuGrid`), with the
background and perturbation sectors built from the same quadrature so the two
are consistent by construction. The warm-dark-matter species of
[the background page](background.md) reuses this machinery unchanged.

## Reference

```@autodocs
Modules = [Cosmic]
Pages   = ["neff.jl"]
Private = false
```
