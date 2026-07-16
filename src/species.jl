"""
Energy components of the universe.

Everything the background knows about a component is captured by one function,
`ρ_over_ρc0(species, a)`: the energy density at scale factor `a` in units of
today's critical density. The expansion rate is then just a sum over species,

    E(a)² = Σ_s ρ_over_ρc0(s, a)

Adding a new component -- early dark energy, a sterile neutrino, a scalar
field, an unusual equation of state -- means adding one method, not touching
the rest of the package. This is why the old FlatLCDM/OpenWCDM/... struct
zoo, which baked Ω_m/Ω_r/Ω_Λ in as fields, is gone.

Pressure `w(species, a)` is supplied where it is well defined; the
perturbation equations need it, and the massive-neutrino case is the one that
actually earns its keep.
"""

using QuadGK: quadgk

abstract type Species end

"""
    ρ_over_ρc0(s::Species, a)

Energy density of `s` at scale factor `a`, in units of the critical density
today. Every species must implement this.
"""
function ρ_over_ρc0 end

"""
    w(s::Species, a)

Equation-of-state parameter P/ρ.
"""
function w end

"Density parameter today, Ω_s = ρ_s(a=1)/ρ_c0."
Ω0(s::Species) = ρ_over_ρc0(s, one(Float64))

# --- Pressureless matter ----------------------------------------------------

"""
    ColdDarkMatter(Ω_c)

Collisionless, pressureless matter. ρ ∝ a⁻³.
"""
struct ColdDarkMatter{T<:Real} <: Species
    Ω::T
end
ρ_over_ρc0(s::ColdDarkMatter, a) = s.Ω / a^3
w(::ColdDarkMatter, a) = 0

"""
    Baryons(Ω_b)

Baryonic matter. Dynamically identical to CDM in the background (ρ ∝ a⁻³) but
tracked separately: it is tightly coupled to the photons before recombination,
and the thermodynamics and perturbation layers need Ω_b on its own.
"""
struct Baryons{T<:Real} <: Species
    Ω::T
end
ρ_over_ρc0(s::Baryons, a) = s.Ω / a^3
w(::Baryons, a) = 0

# --- Radiation --------------------------------------------------------------

"""
    Photons(Ω_γ)

The CMB. ρ ∝ a⁻⁴, w = 1/3. Construct from the CMB temperature with
[`photons_from_Tcmb`](@ref) rather than guessing Ω_γ.
"""
struct Photons{T<:Real} <: Species
    Ω::T
end
ρ_over_ρc0(s::Photons, a) = s.Ω / a^4
w(::Photons, a) = 1 // 3

photons_from_Tcmb(Tcmb, h) = Photons(Constants.Ω_γ(Tcmb, h))

"""
    MasslessNeutrinos(Ω_ν)

A relativistic neutrino background, or any other relativistic relic. ρ ∝ a⁻⁴.

The standard normalisation for `N` effectively-massless species is
Ω_ν = N · (7/8) · (4/11)^(4/3) · Ω_γ, which [`massless_neutrinos`](@ref)
applies for you.
"""
struct MasslessNeutrinos{T<:Real} <: Species
    Ω::T
end
ρ_over_ρc0(s::MasslessNeutrinos, a) = s.Ω / a^4
w(::MasslessNeutrinos, a) = 1 // 3

"Ω_ν for `N_eff` massless species given the photon density."
function massless_neutrinos(N_eff, Ω_γ)
    MasslessNeutrinos(N_eff * (7 / 8) * Constants.T_ν_over_T_γ^4 * Ω_γ)
end

# --- Massive neutrinos ------------------------------------------------------

# --- neutrino phase-space distribution --------------------------------------

"""
    NeutrinoDistribution

The momentum distribution f₀(q) (q = p/T_ν in units of the present neutrino
temperature) that a neutrino species is born with. The default is a Fermi–Dirac
gas; a nonzero degeneracy parameter ξ = μ/T encodes a lepton asymmetry, and
[`TabulatedNu`](@ref) carries an arbitrary numerically-specified f₀. Everything
downstream — the background ρ, p and the free-streaming perturbation hierarchy —
is built from `nu_f0` and `nu_dlnf0`, so the distribution enters the physics in
exactly one place.
"""
abstract type NeutrinoDistribution end

"""
    FermiDirac(ξ = 0)

Fermi–Dirac distribution f₀(q) = 1/(e^{q−ξ} + 1) with degeneracy parameter
ξ = μ/T (a chemical potential / lepton asymmetry). ξ = 0 is the standard case.
"""
struct FermiDirac{T<:Real} <: NeutrinoDistribution
    ξ::T
end
FermiDirac() = FermiDirac(0.0)
@inline nu_f0(d::FermiDirac, q) = 1 / (exp(q - d.ξ) + 1)
@inline nu_dlnf0(d::FermiDirac, q) = -q / (1 + exp(-(q - d.ξ)))

"""
    TabulatedNu(f0, dlnf0)

Arbitrary neutrino distribution given by callables `f0(q)` and its logarithmic
derivative `dlnf0(q) = dln f₀/dln q` (non-thermal relics, sterile mixtures, …).
"""
struct TabulatedNu{F,G} <: NeutrinoDistribution
    f0::F
    dlnf0::G
end
@inline nu_f0(d::TabulatedNu, q) = d.f0(q)
@inline nu_dlnf0(d::TabulatedNu, q) = d.dlnf0(q)

"""
    nu_ρ_integral(y, dist = FermiDirac())

The dimensionless energy integral

    F(y) = ∫₀^∞ dx  x² √(x² + y²) f₀(x),     y = m c² / kT_ν(a),

which controls how a massive neutrino's density evolves. For Fermi–Dirac with
ξ = 0, F(0) = 7π⁴/120 ≈ 5.6822 (relativistic, ρ ∝ a⁻⁴) and F(y→∞) → y·(3/2)ζ(3)·2
(ρ → m·n, ρ ∝ a⁻³). Neutrinos are the one component that crosses between the two,
so no ρ ∝ a⁻ⁿ power law describes them. A generic `dist` (degenerate or tabulated)
just replaces f₀.
"""
function nu_ρ_integral(y, dist::NeutrinoDistribution=FermiDirac())
    # x = u/(1-u) maps [0,∞) onto [0,1); quadgk then hits machine precision fast.
    f(u) = begin
        x = u / (1 - u)
        jac = 1 / (1 - u)^2
        x^2 * sqrt(x^2 + y^2) * nu_f0(dist, x) * jac
    end
    quadgk(f, 0.0, 1.0; rtol=1e-12)[1]
end

"Same integral with √(x²+y²) → x²/√(x²+y²): gives the pressure."
function nu_P_integral(y, dist::NeutrinoDistribution=FermiDirac())
    f(u) = begin
        x = u / (1 - u)
        jac = 1 / (1 - u)^2
        x^4 / sqrt(x^2 + y^2) * nu_f0(dist, x) * jac
    end
    quadgk(f, 0.0, 1.0; rtol=1e-12)[1] / 3
end

# thin backward-compatible wrappers (standard Fermi–Dirac)
fermi_dirac_ρ(y) = nu_ρ_integral(y, FermiDirac())
fermi_dirac_P(y) = nu_P_integral(y, FermiDirac())

const F_massless = 7 * π^4 / 120   # nu_ρ_integral(0, FermiDirac(0)), exactly

"""
    degenerate_Neff_factor(ξ)

Radiation-density enhancement of a degenerate neutrino relative to ξ = 0,
ρ_ν(ξ)/ρ_ν(0) = 1 + (30/7)(ξ/π)² + (15/7)(ξ/π)⁴ (exact for Fermi–Dirac). A
single non-zero ξ therefore adds ΔN_eff = (this − 1) per affected species.
"""
degenerate_Neff_factor(ξ) = 1 + (30 / 7) * (ξ / π)^2 + (15 / 7) * (ξ / π)^4

"""
    MassiveNeutrinos(Ω_rel, y1)

One massive neutrino species. `Ω_rel` is the density parameter it *would* have
if it stayed massless, and `y1 = m c²/(k T_ν,0)` is its mass in units of the
neutrino temperature today.

Density and pressure follow from the phase-space integral:

    ρ(a)/ρ_c0 = Ω_rel · F(y1·a) / F(0) / a⁴

The species behaves as radiation while y1·a ≪ 1 and as matter once y1·a ≫ 1,
with the transition handled exactly rather than by interpolating between two
power laws.
"""
struct MassiveNeutrinos{T<:Real,D<:NeutrinoDistribution} <: Species
    Ω_rel::T   # density it would have if massless (with the ξ enhancement folded in)
    y1::T      # m c² / (k T_ν,0)
    dist::D    # momentum distribution f₀(q)
    F_norm::T  # nu_ρ_integral(0, dist) -- the massless normalisation for this f₀
end
function MassiveNeutrinos(Ω_rel::Real, y1::Real, dist::NeutrinoDistribution=FermiDirac())
    Ωr, y = promote(float(Ω_rel), float(y1))
    Fn = (dist isa FermiDirac && iszero(dist.ξ)) ? F_massless : nu_ρ_integral(0.0, dist)
    MassiveNeutrinos(Ωr, y, dist, oftype(Ωr, Fn))
end

ρ_over_ρc0(s::MassiveNeutrinos, a) = s.Ω_rel * nu_ρ_integral(s.y1 * a, s.dist) / s.F_norm / a^4
P_over_ρc0(s::MassiveNeutrinos, a) = s.Ω_rel * nu_P_integral(s.y1 * a, s.dist) / s.F_norm / a^4
w(s::MassiveNeutrinos, a) = P_over_ρc0(s, a) / ρ_over_ρc0(s, a)

"""
    massive_neutrino(m_eV, Tcmb, h; N_eff_share = 1.0, ξ = 0, dist = nothing)

Build a massive neutrino species of mass `m_eV`.

`N_eff_share` is how much of N_eff this species carries (1.0 for one of three
standard neutrinos). It rescales the effective temperature by
`N_eff_share^(1/4)`, folding the non-instantaneous decoupling correction
(N_eff = 3.044, not 3) into a per-species temperature.

`ξ` is a degeneracy parameter (μ/T, a lepton asymmetry): it uses a
[`FermiDirac`](@ref)`(ξ)` distribution and enhances the relativistic density by
[`degenerate_Neff_factor`](@ref)`(ξ)`. Pass a custom [`TabulatedNu`](@ref) as
`dist` for a fully arbitrary phase space (then `ξ` is ignored).
"""
function massive_neutrino(m_eV, Tcmb, h; N_eff_share=1.0, ξ=0.0, dist=nothing)
    Ωγ = Constants.Ω_γ(Tcmb, h)
    d = dist === nothing ? FermiDirac(ξ) : dist
    # ρ_rel enhanced by the distribution relative to the standard Fermi–Dirac gas.
    enh = (d isa FermiDirac) ? degenerate_Neff_factor(d.ξ) :
          nu_ρ_integral(0.0, d) / F_massless
    Ω_rel = N_eff_share * (7 / 8) * Constants.T_ν_over_T_γ^4 * Ωγ * enh
    # Neutrino temperature today in eV, including the N_eff rescaling.
    T_ν0_eV = Constants.T_ν_over_T_γ * Tcmb * Constants.K_to_eV * N_eff_share^(1 / 4)
    MassiveNeutrinos(Ω_rel, m_eV / T_ν0_eV, d)
end

"""
    warm_dark_matter(m_keV, Ω_wdm, Tcmb, h)

A thermal warm-dark-matter relic: a fermion that froze out relativistic with two
degrees of freedom at its own temperature T_w < T_ν, cool enough to be all of
(or part of) the dark matter today, warm enough that its residual free
streaming erases small-scale structure. Exactly the same phase-space physics as
a massive neutrino — the full momentum-resolved Boltzmann hierarchy — with the
temperature ratio β = T_w/T_ν fixed by the abundance:

    Ω_wdm h² = β³ · m / (93.14 eV)    ⇒    β = (Ω_wdm h² · 93.14 eV / m)^{1/3},

the standard thermal-relic normalisation (the same one CLASS uses for its
`T_ncdm`). Momenta are measured in units of T_w,0, so the distribution is a
clean Fermi–Dirac and `y1 = m/T_w,0`; the relativistic density scales as β⁴.

The free-streaming cutoff this produces can be cross-checked against the
Bode/Viel transfer-function fits — but here it is *computed*, not fitted.
"""
function warm_dark_matter(m_keV, Ω_wdm, Tcmb, h)
    m_eV = 1e3 * float(m_keV)
    Ωγ = Constants.Ω_γ(Tcmb, h)
    T_ν0_eV = Constants.T_ν_over_T_γ * Tcmb * Constants.K_to_eV
    # Today the relic is deeply non-relativistic: ρ = m·n with
    # n ∝ (βT_ν)³·(3ζ(3)/2)/π² for 2 fermion dof, so
    #   Ω_wdm = β³ · m/T_ν0 · (7/8)(T_ν/T_γ)⁴ Ω_γ · [F(∞)/F(0)] with
    #   F(y)/F(0) → y · (3ζ(3)/2)/(7π⁴/120).
    # Solve for β from the code's own constants (not the 93.14 eV literature
    # shortcut, which assumes the instantaneous-decoupling T_ν and is ~1% off
    # against the exact phase-space integral used here).
    coef = (7 / 8) * Constants.T_ν_over_T_γ^4 * Ωγ *
           (1.5 * Constants.ζ3 * 120 / (7 * π^4)) * m_eV / T_ν0_eV
    β = cbrt(Ω_wdm / coef)
    β < 1 || error("warm_dark_matter: T_w/T_ν = $(round(β,digits=3)) ≥ 1 — " *
                   "this mass/abundance needs a hotter-than-neutrino relic; " *
                   "lower Ω_wdm or raise m.")
    Ω_rel = β^4 * (7 / 8) * Constants.T_ν_over_T_γ^4 * Ωγ
    MassiveNeutrinos(Ω_rel, m_eV / (β * T_ν0_eV), FermiDirac(0.0))
end

# --- Decaying dark matter ----------------------------------------------------

"""
    DecayingCDM(Γ, lnρ) / DecayRadiation(lnρ)

Cold dark matter decaying into a relativistic dark radiation with rate Γ (in
1/Mpc, per proper time). Unlike every other species these have no closed-form
density: ρ_dcdm ∝ a⁻³ e^{−Γt} needs t(a), which depends on the full expansion
history — including the decay products themselves. `cosmology()` therefore
solves the coupled system

    d ln ρ_dcdm/dx = −3 − aΓ/ℋ
    d ln ρ_dr/dx   = −4 + (aΓ/ℋ) ρ_dcdm/ρ_dr,      x = ln a, ℋ = aH,

self-consistently (fixed-point in the dark-energy closure) and stores the two
histories as interpolants carried by this species pair.
"""
struct DecayingCDM{I} <: Species
    Γ::Float64            # decay rate, 1/Mpc
    lnρ::I                # x = ln a  ↦  ln(ρ/ρ_c0)
end
struct DecayRadiation{I} <: Species
    lnρ::I
end
ρ_over_ρc0(s::DecayingCDM, a) = exp(s.lnρ(log(a)))
ρ_over_ρc0(s::DecayRadiation, a) = exp(s.lnρ(log(a)))
w(::DecayingCDM, a) = 0
w(::DecayRadiation, a) = 1 // 3

# --- Curvature --------------------------------------------------------------

"""
    Curvature(Ω_k)

Not a fluid, but it enters E(a)² the same way, with ρ_eff ∝ a⁻². Sign
convention: Ω_k > 0 is an open (negatively curved) universe.
"""
struct Curvature{T<:Real} <: Species
    Ω::T
end
ρ_over_ρc0(s::Curvature, a) = s.Ω / a^2
w(::Curvature, a) = -1 // 3

# --- Dark energy ------------------------------------------------------------

abstract type AbstractDarkEnergy <: Species end

"""
    CosmologicalConstant(Ω_Λ)

w = -1 exactly, ρ constant.
"""
struct CosmologicalConstant{T<:Real} <: AbstractDarkEnergy
    Ω::T
end
ρ_over_ρc0(s::CosmologicalConstant, a) = s.Ω
w(::CosmologicalConstant, a) = -1

"""
    W0WaDarkEnergy(Ω_de, w0, wa)

CPL / Chevallier-Polarski-Linder parameterisation, w(a) = w0 + wa(1 - a).
Integrating the continuity equation gives

    ρ(a)/ρ₀ = a^(-3(1 + w0 + wa)) · exp(-3 wa (1 - a))
"""
struct W0WaDarkEnergy{T<:Real} <: AbstractDarkEnergy
    Ω::T
    w0::T
    wa::T
end
function W0WaDarkEnergy(Ω::Real, w0::Real, wa::Real)
    W0WaDarkEnergy(promote(float(Ω), float(w0), float(wa))...)
end

function ρ_over_ρc0(s::W0WaDarkEnergy, a)
    s.Ω * a^(-3 * (1 + s.w0 + s.wa)) * exp(-3 * s.wa * (1 - a))
end
w(s::W0WaDarkEnergy, a) = s.w0 + s.wa * (1 - a)

"""
    GeneralDarkEnergy(Ω_de, w_func)

Dark energy with an arbitrary equation of state `w_func(a)`. The density comes
from integrating the continuity equation numerically,

    ρ(a)/ρ₀ = exp( 3 ∫_a^1 [1 + w(a')] da'/a' )

which is slower than the closed forms above but imposes no restriction on the
shape of w(a) -- early dark energy, quintessence tracker fields, oscillating
models all just work.
"""
struct GeneralDarkEnergy{T<:Real,F} <: AbstractDarkEnergy
    Ω::T
    w_func::F
end

function ρ_over_ρc0(s::GeneralDarkEnergy, a)
    integrand(x) = (1 + s.w_func(x)) / x
    s.Ω * exp(3 * quadgk(integrand, a, one(a); rtol=1e-10)[1])
end
w(s::GeneralDarkEnergy, a) = s.w_func(a)

"""
    QuintessenceDE(Ω, lnρ, w_itp)

Dark energy from a canonical scalar field rolling in a potential V(φ). The
background Klein–Gordon equation is solved in `cosmology()` (with the potential
amplitude shot to close the budget) and the resulting ρ(a), w(a) histories are
carried here as interpolants in x = ln a.

Perturbations come for free and **exactly**: a canonical scalar field is
precisely a fluid with rest-frame sound speed ĉ_s² = 1 (e.g. arXiv:1004.5509
§2), so the validated `fld` machinery in `perturbations.jl` — which picks up
any non-Λ `AbstractDarkEnergy` with cs2_de = 1 by default — evolves the exact
quintessence perturbations with no Klein–Gordon perturbation block needed.
"""
struct QuintessenceDE{I1,I2} <: AbstractDarkEnergy
    Ω::Float64            # realised Ω today
    lnρ::I1               # x = ln a ↦ ln(ρ/ρ_c0)
    w_itp::I2             # x = ln a ↦ w(a)
end
ρ_over_ρc0(s::QuintessenceDE, a) = exp(s.lnρ(log(a)))
w(s::QuintessenceDE, a) = s.w_itp(log(a))
