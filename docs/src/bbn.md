```@meta
CurrentModule = Cosmic
```

# Big Bang nucleosynthesis

Most Boltzmann codes take the helium fraction from an external table. Here the
network is part of the package: the weak ``n \leftrightarrow p`` rates with
their finite-temperature QED and thermal corrections, a reaction network over
the light nuclides with the PRIMAT rate compilation (each of the 63 reaction
rates shipped as a data file with provenance, plus detailed-balance closures),
and the nuclear evolution integrated through the same background the rest of
the code uses — including any exotic components, since the Friedmann equation
is assembled from the species list.

This buys three things. First, consistency: `cosmology()` with the default
`Yp = :bbn` gives every downstream module — recombination, CMB damping tail,
diffusion distortions — a helium abundance that actually corresponds to its
baryon density and expansion history. Second, reach: deuterium, helium-3,
lithium — the observables of BBN — come out of the same call (`bbn`,
`primordial_bbn`), with Monte-Carlo rate-uncertainty propagation in `bbn_mc`.
Third, a gate for new physics: anything that changes the expansion rate or
injects energy (extra radiation from decaying dark matter, a gravitational-wave
background contributing to ``N_{\rm eff}``, lepton asymmetries) feeds through
to abundances that can be checked against observation within the same package.

Two caveats worth knowing. Rate compilations differ: PRIMAT and PArthENoPE
disagree at the 2.6% level on D/H and more on lithium-7, so comparisons should
name the compilation. And the quoted abundances inherit the network's
truncation — the full CNO extension relevant to the lithium problem at the
``10^{-4}`` level is not yet included.

## Reference

```@autodocs
Modules = [Cosmic]
Pages   = ["bbn.jl", "bbn_thermal.jl"]
Private = false
```
