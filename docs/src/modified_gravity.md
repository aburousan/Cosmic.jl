```@meta
CurrentModule = Cosmic
```

# Modified gravity (Horndeski / EFT of dark energy)

Every scalar–tensor theory with second-order equations of motion — Horndeski
gravity — reduces at linear order to a background expansion plus four functions
of time: the kineticity ``\alpha_K``, the braiding ``\alpha_B``, the
Planck-mass run rate ``\alpha_M = d\ln M_*^2/d\ln a``, and the tensor-speed
excess ``\alpha_T`` (Bellini & Sawicki 2014). The GW170817 bound on the speed
of gravitational waves fixes ``\alpha_T = 0``, which is the sub-class Cosmic
implements. Fixing the background to the model already in the `Cosmology` (the
"designer" approach) and modifying only the perturbations isolates the genuine
signatures of the extra scalar: a scale-dependent growth of structure, a
modified late integrated Sachs–Wolfe effect, and a slip between the two metric
potentials.

## The stable basis

Sampling the ``\alpha`` functions directly is treacherous: most choices drift
into ghost or gradient instabilities, and the sound speed reconstructed from
them suffers catastrophic cancellation. Cosmic follows `mochi_class`
(Cataneo & Bellini 2024) and uses the **stable basis** of Kennedy, Lombriser &
Taylor instead, in which the user supplies

```math
D_{\rm kin}(x) > 0, \qquad M_*^2(x) > 0, \qquad c_s^2(x) > 0, \qquad
\alpha_{B0} \equiv \alpha_B(z=0),
```

as callables of ``x = \ln a``. Positivity of the three inputs *is* the
no-ghost/no-gradient stability requirement, satisfied by construction —
``D_{\rm kin} = \alpha_K + \tfrac32\alpha_B^2`` is the de-mixed kinetic term
and ``c_s^2`` the scalar sound speed. The braiding then follows by integrating
the Riccati equation that the exact ``c_s^2`` expression (Bellini & Sawicki
eq. 3.13) defines for ``\alpha_B`` once ``c_s^2`` is prescribed, backward from
``\alpha_B(0) = \alpha_{B0}``; the kineticity is
``\alpha_K = D_{\rm kin} - \tfrac32\alpha_B^2``.

[`stable_basis_solve`](@ref) does this and returns a [`HorndeskiFunctions`](@ref)
carrying the splined ``\alpha_B, \alpha_K, \alpha_M, M_*^2`` and their exact
derivatives. Pass it as `solve_perturbations(c, bg, k; mg = hf)` (or
`cmb_spectra(...; mg = hf)`, `matter_power_spectrum(...; mg = hf)`). With no
`mg` the general-relativistic code path is untouched, bit for bit.

```julia
ΩΛ(x) = # fractional dark-energy density as a function of x = ln a
spec = HorndeskiEFT(D_kin = x -> ΩΛ(x) + 1.5*(0.2ΩΛ(x))^2,
                    M_star2 = x -> exp(0.1 * ∫ΩΛ),   # α_M = 0.1 ΩΛ
                    c_s2 = x -> 0.8, α_B0 = 0.2ΩΛ(0.0))
hf = stable_basis_solve(c, spec)
S  = cmb_spectra(c, rec; mg = hf)
```

## How the scalar is evolved

The Horndeski scalar carries a velocity potential ``v_X`` with its own
second-order equation of motion, reduced from the Bellini & Sawicki
Newtonian-gauge system by algebraically eliminating the second time derivative
of the metric between the pressure equation and the scalar equation — the
kinetic coefficient collects to ``H^2 D_{\rm kin}``, non-zero exactly because
the no-ghost condition holds. The metric potentials ``\phi, \psi`` then follow
algebraically from the momentum and anisotropy constraints (the modified
Poisson sector), and propagate through the existing source functions with no
further change. The scalar stays frozen on the GR equations until its kinetic
term is a meaningful fraction of its peak, mirroring `mochi_class`'s activation
threshold. In the tensor sector the running Planck mass adds Hubble friction,
``\ddot h + (2+\alpha_M)\mathcal{H}\dot h + k^2 h = \text{source}/M_*^2``, while
the wave speed is untouched by the ``\alpha_T = 0`` assumption.

The implementation is validated against `hi_class` v3.2.3: for the
`propto_omega` model the sub-horizon ``P(k)`` enhancement and every ``C_\ell``
agree to ``3\times10^{-5}``, and the gauge-invariant Bardeen potentials agree
to ``\le 2\times10^{-3}`` at all scales, including well outside the horizon.

## Reference

```@autodocs
Modules = [Cosmic]
Pages   = ["modified_gravity.jl"]
Private = false
```
