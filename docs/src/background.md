```@meta
CurrentModule = Cosmic
```

# Background cosmology and species

The homogeneous universe in Cosmic.jl is a `Cosmology`: a container of
*species*, each knowing its own energy density ``\rho(a)/\rho_{c,0}`` and
equation of state ``w(a)``. The Friedmann equation is then just the sum

```math
\frac{H^2(a)}{H_0^2} = \sum_i \frac{\rho_i(a)}{\rho_{c,0}},
```

and every background quantity — distances, ages, horizons — is an integral
over it. Nothing assumes a particular species content, so exotic components
(warm or decaying dark matter, quintessence) enter on exactly the same footing
as photons and baryons.

## Building a cosmology

`cosmology()` assembles the species from familiar parameters. Two choices are
worth spelling out. First, the helium fraction defaults to `Yp = :bbn`: given
``\omega_b`` and ``N_{\rm eff}``, nucleosynthesis fixes ``Y_p``, so treating it
as a free parameter would let a model quietly carry a helium abundance
inconsistent with its own baryon density; the network is solved (and cached)
instead. Second, dark energy closes the budget: ``\Omega_{\rm DE}`` is never an
input but whatever flatness (or the requested ``\Omega_k``) demands.

Massive neutrinos are never approximated by a power law in ``a``. Their
density comes from the exact phase-space integral
``\rho \propto \int dq\, q^2 \sqrt{q^2 + (ma/T_{\nu,0})^2}\, f_0(q)``,
which crosses from radiation-like to matter-like behaviour with no
interpolation.

## Dark energy

Four options, in increasing generality:

- a cosmological constant (the default, `w0 = -1, wa = 0`);
- the CPL fluid `w(a) = w0 + wa (1 - a)`;
- an arbitrary `w_func(a)`, whose density integrates the continuity equation
  numerically;
- a quintessence field: pass a potential `V_scf` (in units of the critical
  density, with the field in reduced Planck masses) and its derivative
  `dV_scf`. The Klein–Gordon equation is solved through the full expansion
  history, with the potential amplitude shot by bisection so the field's
  density today closes the budget exactly. The field starts frozen — the
  early-time attractor under Hubble friction — and its ``w(a)`` is whatever
  the potential produces, not a parameterization.

For anything other than the cosmological constant, the corresponding
perturbations are evolved automatically (see
[Linear perturbations](perturbations.md)); a true ``\Lambda`` carries none,
exactly.

## Dark matter beyond CDM

Warm dark matter (`cosmology(m_wdm = ..., Ω_wdm = ...)`) is a thermal fermion
relic riding the same momentum-resolved machinery as massive neutrinos. Its
temperature ratio to the neutrinos follows from the abundance through the
exact phase-space integral rather than the usual ``93.14\,\mathrm{eV}``
shortcut (which assumes instantaneous-decoupling neutrino temperature and is
about 1% off). The free-streaming suppression of structure it produces is
computed, not fitted; the resulting cutoff reproduces the Viel transfer
function at the few-percent level as an output.

Decaying dark matter (`Γ_dcdm` in km/s/Mpc, `Ω_dcdm_ini` the abundance it
would have today without decay) converts to a dark radiation component through
``\dot\rho_{\rm dcdm} = -3\mathcal{H}\rho_{\rm dcdm} - a\Gamma\rho_{\rm dcdm}``.
No closed form exists — the decay is in proper time, so the history depends on
the expansion it modifies — and the coupled system is solved at construction,
self-consistently with the dark-energy closure. The realised
``\Omega_{\rm dcdm}`` today matches the analytic ``\Omega_{\rm ini}
e^{-\Gamma t_0}`` to five digits, and the full perturbation sector (including
the sourced dark-radiation hierarchy) reproduces the decay response of CLASS
to better than half a percent in the CMB.

Annihilating dark matter affects the thermal history rather than the
background; see [Thermal history](thermal_history.md).

## Distances and times

The usual kinematic toolbox — comoving, transverse comoving, angular-diameter
and luminosity distances, lookback time, age, conformal time — all as direct
integrals over ``E(a)`` with curvature handled exactly in open, flat, and
closed geometries.

## Reference

```@autodocs
Modules = [Cosmic]
Pages   = ["background.jl", "curvature.jl"]
Private = false
```

```@autodocs
Modules = [Cosmic]
Pages   = ["species.jl"]
Private = false
```
