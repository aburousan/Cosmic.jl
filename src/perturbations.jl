"""
Linear perturbations: the Einstein-Boltzmann hierarchy.

Conformal Newtonian gauge, following Ma & Bertschinger (1995). The equations as
implemented -- including the derivation of the adiabatic initial conditions --
are written out in `docs/src/perturbation_equations.md`. Read that alongside this
file: a sign error in a Boltzmann hierarchy is invisible until it surfaces as a
wrong power spectrum, so the derivation is kept where it can be checked.

Each Fourier mode evolves independently (the equations are linear over a
homogeneous background), so "solving" means picking a k and integrating one ODE
system from deep in radiation domination to today.

Time variable is x = ln a: the interesting physics -- horizon crossing, equality,
recombination -- is spread over eight decades in a, and a linear grid would spend
all its effort where nothing happens.
"""

using OrdinaryDiffEq: ODEProblem, solve, KenCarp4
using Interpolations: linear_interpolation, Line
import QuadGK

# --- Background cache -------------------------------------------------------
#
# The right-hand side is evaluated tens of thousands of times per mode, and it
# needs ℋ(a), η(a), the Thomson rate, and the baryon sound speed each time.
# Computing η(a) from its integral inside the RHS -- which is the obvious thing
# to write -- runs a fresh adaptive quadrature per evaluation and costs gigabytes.
# Everything here depends only on the background, not on k, so it is tabulated
# once and shared across every mode.

struct BackgroundCache{C,R,I1,I2,I3,I4}
    cosmo::C
    rec::R
    ℋ::I1        # x = ln a  ->  aH        [1/Mpc]
    η::I2        # x         ->  conformal time [Mpc]
    κ̇::I3        # x         ->  a nₑ σ_T  [1/Mpc]
    c_s2::I4      # x         ->  baryon sound speed squared
end

function BackgroundCache(c::Cosmology, rec::RecombinationSolution;
    x_min=log(1e-10), x_max=0.0, n=6000)
    xs = range(x_min, x_max; length=n)
    as = exp.(xs)

    ℋs = [ℋ(c, a) for a in as]
    ηs = [conformal_time(c, a) for a in as]
    κ̇s = [-τ_dot(rec, a) for a in as]

    # Baryon sound speed (CLASS `thermodynamics.c`, index_th_cb2):
    #
    #     c_s² = k_B T_b / (μ m_H c²) · [ 1 - (1/3) dlnT_b/dlna ]
    #
    # Two things here are easy to get wrong, and both matter at high k because
    # c_s² sets the baryon Jeans scale:
    #
    #  * **μ, the mean molecular weight.** ρ_b/n_total = μ m_H with
    #    μ = (1 + 3.9715 f_He)/(1 + f_He + x_e) -- helium weighs its actual mass in
    #    units of m_H, not its nucleon count. The plasma is fully ionized before
    #    recombination, so every electron counts as a particle: μ ≈ 0.59. After
    #    recombination the electrons are gone and μ ≈ 1.23. Using μ = 1 makes c_s²
    #    ~1.7x too small in the era that actually matters.
    #
    #  * **The bracket is not a constant.** It equals 4/3 only while T_b ∝ 1/a
    #    (tightly Compton-coupled). Once the gas decouples thermally it cools as
    #    1/a² and the bracket becomes 5/3. Taking it from the actual T_b(a) rather
    #    than hard-coding 4/3 keeps both regimes right.
    fHe = f_He(c)
    c_s2s = similar(as)
    for (i, a) in enumerate(as)
        z = redshift(a)
        T_b = T_matter(rec, z)
        xe = x_e(rec, z)
        μ = (1 + 3.9715 * fHe) / (1 + fHe + xe)
        # dlnT_b/dlna by central difference in ln a
        δ = 1e-3
        Tp = T_matter(rec, redshift(a * exp(δ)))
        Tm = T_matter(rec, redshift(a * exp(-δ)))
        dlnT_dlna = (log(Tp) - log(Tm)) / (2δ)
        c_s2s[i] = Constants.k_B_SI * T_b /
                   (μ * Constants.m_H_SI * Constants.c_SI^2) * (1 - dlnT_dlna / 3)
    end

    BackgroundCache(c, rec,
        linear_interpolation(xs, ℋs; extrapolation_bc=Line()),
        linear_interpolation(xs, ηs; extrapolation_bc=Line()),
        linear_interpolation(xs, κ̇s; extrapolation_bc=Line()),
        linear_interpolation(xs, c_s2s; extrapolation_bc=Line()),
    )
end

# --- Massive-neutrino momentum grid -----------------------------------------

"""
    MassiveNuGrid

Gauss-Legendre nodes in comoving momentum q, with the Fermi-Dirac weight folded in.

Massive neutrinos are the one species that cannot be reduced to a fluid. A fluid
has a single velocity at each point. A free-streaming gas of *massive* particles
has particles of different momenta moving at different speeds through the same
point, and that spread in velocity is precisely what erases small-scale structure.
So the perturbation has to be tracked as a function of momentum -- a whole
hierarchy per momentum bin -- which is why this is expensive and why it is the last
thing to get implemented in any Boltzmann code.
"""
struct MassiveNuGrid
    q::Vector{Float64}
    w::Vector{Float64}       # Gauss-Legendre weights
    f0::Vector{Float64}      # 1/(e^q + 1)
    dlnf0::Vector{Float64}   # dln f₀/dln q = -q/(1 + e^{-q})
end

function MassiveNuGrid(nq::Int; qmax=25.0, dist::NeutrinoDistribution=FermiDirac())
    q, w = QuadGK.gauss(nq, 0.0, qmax)
    f0 = [nu_f0(dist, qi) for qi in q]
    dlnf0 = [nu_dlnf0(dist, qi) for qi in q]
    MassiveNuGrid(collect(q), collect(w), f0, dlnf0)
end

"Background density integral ∫dq q²ε f₀ on the grid, in units where ρ ∝ this/a⁴."
function _Iρ(G::MassiveNuGrid, am)
    s = 0.0
    @inbounds for i in eachindex(G.q)
        ε = sqrt(G.q[i]^2 + am^2)
        s += G.w[i] * G.q[i]^2 * ε * G.f0[i]
    end
    s
end

# --- State vector layout ----------------------------------------------------
#
# 1              δ_c
# 2              θ_c
# 3              δ_b
# 4              θ_b
# 5              φ                       (evolved, not reconstructed -- see below)
# iγ .. iγ+lγ    F_γ0 .. F_γ,lmax_γ      photon temperature
# iG .. iG+lγ    G_γ0 .. G_γ,lmax_γ      photon polarization
# iν .. iν+lν    F_ν0 .. F_ν,lmax_ν      massless neutrinos
# iM ..          Ψ_ℓ(q) for each massive species, each q node: the momentum-
#                resolved hierarchy. Indexed by _imν below.

struct Layout
    lmax_γ::Int
    lmax_ν::Int
    lmax_m::Int
    nq::Int
    nmν::Int          # number of massive neutrino species
    iγ::Int
    iG::Int
    iν::Int
    iM::Int
    # decaying DM + its dark radiation. The block sits at the very end of the
    # state so that a cosmology without dcdm has exactly the layout it always
    # had. idc = δ_dcdm, idc+1 = θ_dcdm; idr..idr+lmax_dr = density-weighted
    # dark-radiation moments F_ℓ (CLASS convention: F carries ρ_dr a⁴/ρ_c0).
    has_dcdm::Bool
    lmax_dr::Int
    idc::Int
    idr::Int
    # dark-energy fluid perturbations (w ≠ −1): δ_fld, θ_fld at the very end
    has_fld::Bool
    ifld::Int
    n::Int
end

function Layout(lmax_γ, lmax_ν; lmax_m=8, nq=0, nmν=0, has_dcdm=false, lmax_dr=17,
    has_fld=false)
    iγ = 6
    iG = iγ + (lmax_γ + 1)
    iν = iG + (lmax_γ + 1)
    iM = iν + (lmax_ν + 1)
    n_mν = iM - 1 + nmν * nq * (lmax_m + 1)
    idc = n_mν + 1
    idr = idc + 2
    n_dc = has_dcdm ? (idr + lmax_dr) : n_mν
    ifld = n_dc + 1
    n = has_fld ? (ifld + 1) : n_dc
    Layout(lmax_γ, lmax_ν, lmax_m, nq, nmν, iγ, iG, iν, iM,
        has_dcdm, lmax_dr, idc, idr, has_fld, ifld, n)
end

"Index of Ψ_ℓ(q_iq) for massive species `s`."
@inline _imν(L::Layout, s, iq, ℓ) =
    L.iM + ((s - 1) * L.nq + (iq - 1)) * (L.lmax_m + 1) + ℓ

"""
    PerturbationSolution

The solved hierarchy for one wavenumber `k` (1/Mpc).
"""
struct PerturbationSolution{C,R,S}
    cosmo::C
    rec::R
    k::Float64
    sol::S
    lmax_γ::Int
    lmax_ν::Int
    lmax_m::Int
    nq::Int
    nmν::Int
end

_layout(p::PerturbationSolution) =
    Layout(p.lmax_γ, p.lmax_ν; lmax_m=p.lmax_m, nq=p.nq, nmν=p.nmν,
        has_dcdm=get_species(p.cosmo, DecayingCDM) !== nothing,
        has_fld=_fld_species(p.cosmo) !== nothing)

# Decaying-DM parameters threaded through the RHS: the decay rate and the two
# background ln-density histories (needed at every step for the source terms).
_dcdm_params(c::Cosmology) = begin
    dc = get_species(c, DecayingCDM)
    dc === nothing ? nothing :
    (Γ=dc.Γ, lnρd=dc.lnρ, lnρr=get_species(c, DecayRadiation).lnρ)
end

# The dark-energy species whose perturbations must be evolved: any non-Λ fluid.
# A cosmological constant has no perturbations (exactly — (1+w) = 0).
_fld_species(c::Cosmology) = begin
    de = get_species(c, AbstractDarkEnergy)
    de isa CosmologicalConstant ? nothing : de
end

_dwda(s::W0WaDarkEnergy, a) = -s.wa
# relative step: an absolute one underflows a < h to negative scale factors
# (the solver starts at a ~ 1e-9)
_dwda(s::GeneralDarkEnergy, a; ε=1e-4) =
    (w(s, a * (1 + ε)) - w(s, a * (1 - ε))) / (2ε * a)
_dwda(s::QuintessenceDE, a; ε=1e-4) =
    (w(s, a * (1 + ε)) - w(s, a * (1 - ε))) / (2ε * a)

# Does w(a) cross −1 anywhere the perturbations are evolved? (A canonical
# scalar field cannot — 1+w = φ̇²/ρ ≥ 0 — so quintessence always says no.)
_crosses_w_minus1(s::QuintessenceDE) = false
function _crosses_w_minus1(s::AbstractDarkEnergy)
    lo, hi = false, false
    for x in range(log(1e-9), 0.0; length=200)
        opw = 1 + w(s, exp(x))
        opw > 1e-8 && (hi = true)
        opw < -1e-8 && (lo = true)
        lo && hi && return true
    end
    false
end

# --- Metric ------------------------------------------------------------------
#
# φ is *evolved*, not reconstructed. The tempting alternative -- eliminate φ̇
# between the energy and momentum constraints and get φ algebraically -- is
# numerically fatal outside the horizon: there Δ and 3ℋΘ/k² agree to several
# significant figures and cancel, so the answer is the difference of two large
# nearly-equal numbers. Evolving the energy constraint has no such cancellation.
#
#   φ̇ = [ -(3/2) H₀²a² Δ  -  k²φ  -  3ℋ²ψ ] / (3ℋ)
#   ψ = φ - (9/2)(H₀²a²/k²) Σ_σ
#
# ψ is safe to get algebraically: the anisotropic stress Σ_σ is a small
# correction, not a cancelling one.
#
# Sanity check of the sign: superhorizon in radiation domination, δ_tot → -2ψ and
# ψ → φ, so (3/2)H₀²a²Δ = (3/2)ℋ²δ_tot = -3ℋ²ψ and φ̇ → 0. φ is constant, as it
# must be.

"""
    _massive_ν_moments(u, L, G, k, a, mν)

Integrate the momentum-resolved hierarchy back down to the three numbers the
Einstein equations actually want: δρ, (ρ̄+P̄)θ, (ρ̄+P̄)σ, all in units of ρ_c0.

    δρ      ∝ (1/a⁴) ∫dq q² ε    f₀ Ψ₀
    (ρ̄+P̄)θ ∝ (k/a⁴) ∫dq q³      f₀ Ψ₁
    (ρ̄+P̄)σ ∝ (2/3a⁴) ∫dq q⁴/ε   f₀ Ψ₂

normalised by the same ∫dq q²ε f₀ the *background* uses, so the two are consistent
by construction rather than by coincidence.
"""
@inline function _massive_ν_moments(u, L, G::MassiveNuGrid, k, a, mν)
    δρ = 0.0
    ρPθ = 0.0
    ρPσ = 0.0
    for (s, m) in enumerate(mν)
        am = m.y1 * a
        norm = m.Ω_rel / F_massless / a^4       # matches species.jl exactly
        Iρ = 0.0
        Iδ = 0.0
        Iθ = 0.0
        Iσ = 0.0
        @inbounds for iq in 1:L.nq
            q = G.q[iq]
            ε = sqrt(q^2 + am^2)
            wf = G.w[iq] * G.f0[iq]
            Ψ0 = u[_imν(L, s, iq, 0)]
            Ψ1 = u[_imν(L, s, iq, 1)]
            Ψ2 = u[_imν(L, s, iq, 2)]
            Iρ += wf * q^2 * ε
            Iδ += wf * q^2 * ε * Ψ0
            Iθ += wf * q^3 * Ψ1
            Iσ += wf * q^4 / ε * Ψ2
        end
        δρ += norm * Iδ
        ρPθ += norm * k * Iθ
        ρPσ += norm * (2 / 3) * Iσ
    end
    (δρ, ρPθ, ρPσ)
end

@inline function _metric(u, L, k, K, a, ℋ_, H0, ργ, ρν, ρb, ρc, G, mν,
    ρdc=0.0, ρdr=0.0, ρfld=0.0, wfld=-1.0, ppf=false, cs2f=1.0, Pmν=0.0)
    δ_γ = u[L.iγ]
    θ_γ = 3k / 4 * u[L.iγ+1]
    s2 = _curvature_streaming(K, k, 2)
    σ_γ = u[L.iγ+2] / (2s2)
    δ_ν = u[L.iν]
    θ_ν = 3k / 4 * u[L.iν+1]
    σ_ν = u[L.iν+2] / (2s2)
    φ = u[5]

    pre = 1.5 * H0^2 * a^2 / k^2
    Δ = ρc * u[1] + ρb * u[3] + ργ * δ_γ + ρν * δ_ν
    Θtot = ρc * u[2] + ρb * u[4] + (4 / 3) * ργ * θ_γ + (4 / 3) * ρν * θ_ν
    Σσ = (4 / 3) * ργ * σ_γ + (4 / 3) * ρν * σ_ν

    if L.nmν > 0
        dρ_m, ρPθ_m, ρPσ_m = _massive_ν_moments(u, L, G, k, a, mν)
        Δ += dρ_m
        Θtot += ρPθ_m
        Σσ += ρPσ_m
    end

    if L.has_dcdm
        # dcdm is a plain matter fluid; the dark radiation is carried as CLASS's
        # density-weighted moments F_ℓ = (physical moment)·ρ_dr a⁴/ρ_c0, so its
        # Einstein-sum entries are δρ = F0/a⁴, (ρ+p)θ = kF1·(3/4)(4/3)/a⁴ = ...
        # δρ_dr = F0/a⁴, (ρ+p)θ_dr = k F1/a⁴, (ρ+p)σ_dr = (2/3) F2/a⁴.
        Δ += ρdc * u[L.idc] + u[L.idr] / a^4
        Θtot += ρdc * u[L.idc+1] + k * u[L.idr+1] / a^4
        Σσ += (2 / 3) * u[L.idr+2] / a^4
    end

    γp = 0.0
    if L.has_fld
        if ppf
            # PPF (Fang–Hu–Lewis 0808.3125; CLASS perturbations.c 7180-7265,
            # Newtonian branch): one evolved variable Γ = u[ifld]; δρ_fld and
            # (ρ+p)θ_fld are algebraic in Γ and the everything-else sums (which
            # is why this block sits after all other species). Units: CLASS
            # densities are H0²·(ρ/ρ_c0), hence the explicit H0² factors.
            # (ρ̃+p̃)_tot without the fld; Pmν = Σ_mν(ρ+P) passed from the caller
            ρpt = (4 / 3) * (ργ + ρν + ρdr) + ρb + ρc + ρdc + Pmν
            θt = Θtot / ρpt
            x = (0.4 * k / ℋ_)^2 * cs2f
            Γf = x > 1e3 ? 0.0 : u[L.ifld]
            S = pre * ρfld * (1 + wfld) * θt / ℋ_
            γp = x > 1e3 ? 0.0 : ℋ_ * (S / (1 + x) - (1 + x) * Γf)
            Y = 4.5 * a^2 * H0^2 * ρpt / k^2
            ρpθf = ρfld * (1 + wfld) * θt -
                   (2 / 3) * k^2 * ℋ_ / (a^2 * H0^2) / (1 + Y) *
                   (S - γp / ℋ_ - Γf)
            s2 = _curvature_streaming(K, k, 2)
            δρf = -(2 / 3) * k^2 * s2^2 * Γf / (a^2 * H0^2) - 3 * ℋ_ * ρpθf / k^2
            Δ += δρf
            Θtot += ρpθf
        else
            # exact fluid in (δ, Q = (1+w)θ): no anisotropic stress,
            # (ρ+p)θ = ρ·Q directly.
            Δ += ρfld * u[L.ifld]
            Θtot += ρfld * u[L.ifld+1]
        end
    end

    ψ = φ - 3 * pre * Σσ
    φ̇ = (-1.5 * H0^2 * a^2 * Δ - (k^2 - 3K) * φ - 3 * ℋ_^2 * ψ) / (3 * ℋ_)
    (φ, ψ, φ̇, θ_γ, σ_γ, δ_γ, θ_ν, σ_ν, δ_ν, γp)
end

function _rhs!(du, u, p, x)
    bg, k, K, L, G, mν, dcdm, fld = p
    c = bg.cosmo
    a = exp(x)

    ℋ_ = bg.ℋ(x)
    η = bg.η(x)
    κ̇ = bg.κ̇(x)
    c_s2 = bg.c_s2(x)
    H0 = Constants.H0_in_invMpc(c.h)

    ργ = Ω_γ(c) / a^4
    ρν = _Ω_or_zero(c, MasslessNeutrinos) / a^4
    ρb = Ω_b(c) / a^3
    ρc = Ω_c(c) / a^3
    ρdc = dcdm === nothing ? 0.0 : exp(dcdm.lnρd(x))
    ρdr = dcdm === nothing ? 0.0 : exp(dcdm.lnρr(x))
    ρfld = fld === nothing ? 0.0 : ρ_over_ρc0(fld.de, a)
    wfld = fld === nothing ? -1.0 : w(fld.de, a)
    use_ppf = fld !== nothing && fld.ppf
    cs2f = fld === nothing ? 1.0 : fld.cs2
    # background (ρ+P) of the massive neutrinos, needed only by the PPF S-term
    Pmν = (use_ppf && !isempty(mν)) ?
          sum(ρ_over_ρc0(m, a) + P_over_ρc0(m, a) for m in mν) : 0.0

    φ, ψ, φ̇, θ_γ, σ_γ, δ_γ, θ_ν, σ_ν, δ_ν, γp =
        _metric(u, L, k, K, a, ℋ_, H0, ργ, ρν, ρb, ρc, G, mν, ρdc, ρdr,
            ρfld, wfld, use_ppf, cs2f, Pmν)

    δ_c, θ_c, δ_b, θ_b = u[1], u[2], u[3], u[4]
    iγ, iG, iν = L.iγ, L.iG, L.iν
    F_γ2 = u[iγ+2]
    Π = F_γ2 + u[iG] + u[iG+2]     # couples temperature to polarization

    R = 4 * ργ / (3 * ρb)          # photon-to-baryon momentum ratio
    inv = 1 / ℋ_                   # d/dx = (1/ℋ) d/dη
    s(ell) = _curvature_streaming(K, k, ell)
    cotend = _curvature_cot_closure(K, k, η)

    du[1] = inv * (-θ_c + 3 * φ̇)
    du[2] = inv * (-ℋ_ * θ_c + k^2 * ψ)
    du[3] = inv * (-θ_b + 3 * φ̇)
    du[4] = inv * (-ℋ_ * θ_b + c_s2 * k^2 * δ_b + k^2 * ψ + κ̇ * R * (θ_γ - θ_b))
    du[5] = inv * φ̇

    # ---- Photons: temperature ----------------------------------------------
    du[iγ] = inv * (-(4 / 3) * θ_γ + 4 * φ̇)
    dθ_γ = k^2 * (δ_γ / 4 - s(2)^2 * σ_γ) + k^2 * ψ + κ̇ * (θ_b - θ_γ)
    du[iγ+1] = inv * (4 / (3k)) * dθ_γ                       # F_γ1 = (4/3k) θ_γ
    du[iγ+2] = inv * ((8 / 15) * s(2) * θ_γ - (3 / 5) * k * s(3) * u[iγ+3] +
                       κ̇ * (-F_γ2 + Π / 10))
    for ℓ in 3:(L.lmax_γ-1)
        du[iγ+ℓ] = inv * (k / (2ℓ + 1) *
            (ℓ * s(ℓ) * u[iγ+ℓ-1] - (ℓ + 1) * s(ℓ + 1) * u[iγ+ℓ+1]) -
            κ̇ * u[iγ+ℓ])
    end
    # Truncation (M&B eq. 65): let power leave the top of the hierarchy rather
    # than reflect back down it.
    let ℓ = L.lmax_γ
        du[iγ+ℓ] = inv * (k * (s(ℓ) * u[iγ+ℓ-1] - (ℓ + 1) * cotend * u[iγ+ℓ]) -
                           κ̇ * u[iγ+ℓ])
    end

    # ---- Photons: polarization ---------------------------------------------
    du[iG] = inv * (-k * u[iG+1] + κ̇ * (-u[iG] + Π / 2))
    for ℓ in 1:(L.lmax_γ-1)
        src = (ℓ == 2) ? Π / 10 : 0.0
        du[iG+ℓ] = inv * (k / (2ℓ + 1) *
                          (ℓ * s(ℓ) * u[iG+ℓ-1] - (ℓ + 1) * s(ℓ + 1) * u[iG+ℓ+1]) +
                          κ̇ * (-u[iG+ℓ] + src))
    end
    let ℓ = L.lmax_γ
        du[iG+ℓ] = inv * (k * (s(ℓ) * u[iG+ℓ-1] - (ℓ + 1) * cotend * u[iG+ℓ]) -
                           κ̇ * u[iG+ℓ])
    end

    # ---- Massless neutrinos: the same hierarchy, but free-streaming (κ̇ = 0) --
    du[iν] = inv * (-(4 / 3) * θ_ν + 4 * φ̇)
    du[iν+1] = inv * (4 / (3k)) * (k^2 * (δ_ν / 4 - s(2)^2 * σ_ν) + k^2 * ψ)
    du[iν+2] = inv * ((8 / 15) * s(2) * θ_ν - (3 / 5) * k * s(3) * u[iν+3])
    for ℓ in 3:(L.lmax_ν-1)
        du[iν+ℓ] = inv * (k / (2ℓ + 1) *
            (ℓ * s(ℓ) * u[iν+ℓ-1] - (ℓ + 1) * s(ℓ + 1) * u[iν+ℓ+1]))
    end
    let ℓ = L.lmax_ν
        du[iν+ℓ] = inv * k * (s(ℓ) * u[iν+ℓ-1] - (ℓ + 1) * cotend * u[iν+ℓ])
    end

    # ---- Massive neutrinos: one hierarchy per momentum bin -------------------
    # The factor q/ε is the particle's *velocity*. When q ≫ am̃ it is 1 and these
    # reduce exactly to the massless equations above; when q ≪ am̃ it tends to 0,
    # the free-streaming term shuts off, and the neutrinos start to cluster like
    # matter. Nothing switches over by hand -- the transition is in q/ε, and it
    # happens at a different time for each momentum bin. That spread is the whole
    # reason a fluid approximation fails here.
    if L.nmν > 0
        for (jν, m) in enumerate(mν)
            am = m.y1 * a
            @inbounds for iq in 1:L.nq
                q = G.q[iq]
                ε = sqrt(q^2 + am^2)
                dlnf0 = G.dlnf0[iq]
                qkε = q * k / ε

                i0 = _imν(L, jν, iq, 0)
                du[i0] = inv * (-qkε * u[i0+1] - φ̇ * dlnf0)
                du[i0+1] = inv * (qkε / 3 * (u[i0] - 2 * s(2) * u[i0+2]) -
                                  (ε * k / (3q)) * ψ * dlnf0)
                for ℓ in 2:(L.lmax_m-1)
                    du[i0+ℓ] = inv * (qkε / (2ℓ + 1) *
                                      (ℓ * s(ℓ) * u[i0+ℓ-1] -
                                       (ℓ + 1) * s(ℓ + 1) * u[i0+ℓ+1]))
                end
                let ℓ = L.lmax_m
                    du[i0+ℓ] = inv * (qkε * s(ℓ) * u[i0+ℓ-1] -
                                      (ℓ + 1) * k * cotend * u[i0+ℓ])
                end
            end
        end
    end

    # ---- Decaying DM and its dark radiation ---------------------------------
    # Transcribed from CLASS 3.3.4 perturbations.c (~9305-9355), Newtonian gauge
    # (metric_continuity = −3φ', metric_euler = k²ψ, metric_shear = 0). The dr
    # moments are density-weighted, F_ℓ = (physical)·ρ_dr a⁴/ρ_c0, so the decay
    # source enters linearly through f' = aΓρ_dcdm a⁴ (per conformal time).
    if L.has_dcdm && dcdm !== nothing
        idc, idr = L.idc, L.idr
        aΓ = a * dcdm.Γ
        f = ρdr * a^4
        fp = aΓ * ρdc * a^4                       # df/dη = aΓ ρ_dcdm a⁴

        du[idc] = inv * (-u[idc+1] + 3 * φ̇ - aΓ * ψ)
        du[idc+1] = inv * (-ℋ_ * u[idc+1] + k^2 * ψ)

        du[idr] = inv * (-k * u[idr+1] + 4 * φ̇ * f + fp * (u[idc] + ψ))
        du[idr+1] = inv * (k / 3 * u[idr] - (2 / 3) * k * s(2)^2 * u[idr+2] +
                           (4 / (3k)) * k^2 * ψ * f + (fp / k) * u[idc+1])
        du[idr+2] = inv * ((8 / 15) * (3 / 4) * k * u[idr+1] -
                           (3 / 5) * k * (s(3) / s(2)) * u[idr+3])
        du[idr+3] = inv * (k / 7 * (3 * s(3) * s(2) * u[idr+2] -
                                    4 * s(4) * u[idr+4]))
        for ℓ in 4:(L.lmax_dr-1)
            du[idr+ℓ] = inv * (k / (2ℓ + 1) *
                               (ℓ * s(ℓ) * u[idr+ℓ-1] - (ℓ + 1) * s(ℓ + 1) * u[idr+ℓ+1]))
        end
        let ℓ = L.lmax_dr
            du[idr+ℓ] = inv * (k * (s(ℓ) * u[idr+ℓ-1] -
                                    (1 + ℓ) * cotend * u[idr+ℓ]))
        end
    end

    # ---- Dark-energy fluid ---------------------------------------------------
    # The CLASS fld equations (perturbations.c 9358-9390; verified against
    # Ma-Bertschinger + the rest-frame cs² relation, and against 1004.5509 eqs
    # 3.5-3.6) rewritten in the regular variables (δ, Q ≡ (1+w)θ). The (δ,θ)
    # form is singular at w = −1 through cs²k²δ/(1+w) — fatal for quintessence,
    # which *starts* frozen at exactly w = −1. In (δ,Q) the algebra is identical
    # for w ≠ −1 (a linear change of variables) and every term is regular at
    # w = −1: the lone w′θ = w′Q/(1+w) piece vanishes identically for a frozen
    # field (w pinned at −1 ⇒ w′ = 0), so it is set to zero there.
    if L.has_fld && fld !== nothing
        ifld = L.ifld
        if use_ppf
            # PPF: only Γ is evolved (γp = Γ′ computed in _metric); the second
            # slot is unused in this scheme.
            du[ifld] = inv * γp
            du[ifld+1] = 0.0
        else
            wp = _dwda(fld.de, a) * a * ℋ_            # dw/dη
            δf = u[ifld]
            Qf = u[ifld+1]                             # Q = (1+w)θ
            opw = 1 + wfld
            W = abs(opw) > 1e-10 ? wp * Qf / opw : 0.0 # w′θ, frozen ⇒ 0
            du[ifld] = inv * (-Qf + 3 * opw * φ̇ -
                              3 * ℋ_ * (cs2f - wfld) * δf -
                              9 * ℋ_^2 * (cs2f - wfld) * Qf / k^2 - 3 * ℋ_ * W / k^2)
            du[ifld+1] = inv * (-(1 - 3 * cs2f) * ℋ_ * Qf +
                                k^2 * cs2f * δf + opw * k^2 * ψ + W)
        end
    end

    nothing
end

# --- Initial conditions -----------------------------------------------------

"""
    _isocurvature_ic!(u, c, bg, k, L, x0, mν, mode)

Fill the state vector `u` with one of the four regular isocurvature modes, each
normalised to a unit initial entropy perturbation:

  - `:cdi` — CDM density isocurvature (S_c = δ_c − ¾δ_γ = 1),
  - `:bi`  — baryon density isocurvature,
  - `:nid` — neutrino density isocurvature,
  - `:niv` — neutrino velocity isocurvature.

The amplitudes are the Bucher–Moodley–Turok (astro-ph/9904231) leading-order
super-horizon solutions, valid at leading order in kη and ωη with ω ≡ a ρ_m/√ρ_r
(so ω = H₀ Ω_m/√Ω_r, the parameter that governs a(η) through the
radiation-to-matter transition). They are written first in the synchronous gauge
— exactly as in CLASS `perturbations.c` — and then rotated to the conformal
Newtonian gauge that Cosmic integrates in, via the shift α = (h'+6η')/2k², whose
value follows from the two Einstein constraints:

    α = [η_syn + (3/2)(ℋ²/k²s₂²)(δ_tot + 3ℋ θ_tot/k²)] / ℋ

with δ_tot, θ_tot the total (radiation+matter) density and velocity. Each species
then transforms as δ → δ − (1+w)·3ℋα and θ → θ + k²α; the shear and the ℓ=3
neutrino moment are gauge invariant. The massive-neutrino tower is reshaped from
the massless answer by the caller, identically to the adiabatic case.
"""
function _isocurvature_ic!(u, c::Cosmology, bg::BackgroundCache, k, L::Layout, x0,
    mν, mode::Symbol)
    η = bg.η(x0)
    a = exp(x0)
    ℋ = bg.ℋ(x0)
    K = spatial_curvature_K(c)
    s2sq = 1 - 3K / k^2
    s2 = _curvature_streaming(K, k, 2)

    _ρ(S) = (sp = get_species(c, S); sp === nothing ? 0.0 : ρ_over_ρc0(sp, a))
    ργ = _ρ(Photons)
    ρν = _ρ(MasslessNeutrinos) + sum(ρ_over_ρc0(m, a) for m in mν; init=0.0)
    ρr = ργ + ρν
    ρb = _ρ(Baryons)
    ρcdm = _ρ(ColdDarkMatter)
    ρm = ρb + ρcdm
    H0 = Constants.H0_in_invMpc(c.h)

    fnu = ρν / ρr
    fg = ργ / ρr
    fb = ρm > 0 ? ρb / ρm : 0.0
    fcdm = ρm > 0 ? ρcdm / ρm : 0.0
    rmr = ρm / ρr                       # ρ_m/ρ_r at the start time
    om = a * H0 * ρm / sqrt(ρr)         # ω parameter (1/Mpc)

    kt2 = (k * η)^2
    kt3 = (k * η)^3
    ot = om * η

    δcdm = δb = δg = δur = θg = θb = θur = shear_ur = eta = 0.0
    l3_ur = 0.0
    if mode === :cdi || mode === :bi
        f = mode === :cdi ? fcdm : fb
        δg = f * om * η * (-2 / 3 + ot / 4)
        θg = -f * om * kt2 / 12
        δur = δg
        θur = θg
        θb = θg
        shear_ur = -f * kt2 * η * om / 6 / (2fnu + 15)
        eta = -f * om * η * (1 / 6 - ot / 16)
        if mode === :cdi
            δcdm = 1.0 + 0.75 * δg
            δb = 0.75 * δg
        else
            δb = 1.0 + 0.75 * δg
            δcdm = 0.75 * δg
        end
    elseif mode === :nid
        δg = (fnu / fg) * (-1 + kt2 / 6)
        θg = -(fnu / fg) * k^2 * η * (0.25 - (fb / fg) * (3 / 16) * ot)
        δb = (fnu / fg) / 8 * kt2
        θb = θg
        δcdm = -(fnu * fb / fg / 80) * kt2 * ot
        δur = 1 - kt2 / 6
        θur = k^2 * η / 4
        shear_ur = kt2 / (4fnu + 15) / 2
        eta = -(fnu / (4fnu + 15) / 6) * kt2
    elseif mode === :niv
        δg = k * η * (fnu / fg) * (1 - (3 / 16) * fb * (2 + fg) / fg * ot)
        θg = (fnu / fg) * 0.75 * k *
             (-1 + 0.75 * (fb / fg) * ot + (3 / 16) * ot^2 * fb / fg^2 * (fg - 3fb) + kt2 / 6)
        δb = 0.75 * δg
        θb = θg
        δcdm = -(9 / 64) * fnu * fb / fg * k * η * ot
        δur = -k * η * (1 + (3 / 16) * fb * fnu / fg * ot)
        θur = 0.75 * k * (1 - (1 / 6) * kt2 * (4fnu + 9) / (4fnu + 5))
        shear_ur = k * η / (4fnu + 15) * (1 + 3 * ot * fnu / (4fnu + 15))
        eta = fnu * k * η *
              (-1 / (4fnu + 5) + (-(3 / 64) * fb / fg + (15 / 4) / (4fnu + 15) / (4fnu + 5) * ot))
    else
        error("unknown isocurvature mode: $mode (expected :cdi, :bi, :nid, or :niv)")
    end

    # synchronous -> conformal Newtonian gauge
    delta_tot = (fg * δg + fnu * δur + rmr * (fb * δb + fcdm * δcdm)) / (1 + rmr)
    velocity_tot = ((4 / 3) * (fg * θg + fnu * θur) + rmr * fb * θb) / (1 + rmr)
    α = (eta + 1.5 * ℋ^2 / k^2 / s2sq * (delta_tot + 3 * ℋ / k^2 * velocity_tot)) / ℋ
    φ = eta - ℋ * α

    δg -= 4ℋ * α
    θg += k^2 * α
    δb -= 3ℋ * α
    θb += k^2 * α
    δcdm -= 3ℋ * α
    θcdm = k^2 * α
    δur -= 4ℋ * α
    θur += k^2 * α

    u[1] = δcdm
    u[2] = θcdm
    u[3] = δb
    u[4] = θb
    u[5] = φ
    u[L.iγ] = δg
    u[L.iγ+1] = (4 / (3k)) * θg
    u[L.iν] = δur
    u[L.iν+1] = (4 / (3k)) * θur
    u[L.iν+2] = 2 * s2 * shear_ur
    L.lmax_ν >= 3 && (u[L.iν+3] = l3_ur)
    return u
end

"""
    initial_conditions(c, bg, k, L, x0; ic = :adiabatic)

Adiabatic growing mode, deep in radiation domination and far outside the horizon,
normalised to a primordial curvature perturbation ℛ = 1. Since the system is
linear, the primordial spectrum is applied afterwards rather than here.

Derived in `docs/src/perturbation_equations.md`. With f_ν = ρ_ν/(ρ_γ+ρ_ν):

    ψ = 1/(3/2 + 2f_ν/5),  φ = (1 + 2f_ν/5) ψ
    δ_γ = δ_ν = -2ψ,       δ_c = δ_b = -(3/2)ψ
    θ (all species) = (1/2) k² η ψ
    σ_ν = (1/15) k² η² ψ
"""
function initial_conditions(c::Cosmology, bg::BackgroundCache, k, L::Layout, x0,
    G::MassiveNuGrid, mν; ic::Symbol=:adiabatic, fld_ppf::Bool=false)
    η = bg.η(x0)
    a = exp(x0)

    Ωγ = Ω_γ(c)
    # At the start every neutrino is ultra-relativistic, massive ones included, so
    # they all count toward the radiation fraction that sets f_ν.
    Ων = _Ω_or_zero(c, MasslessNeutrinos) +
         sum(m.Ω_rel for m in mν; init=0.0)
    f_ν = Ων / (Ωγ + Ων)

    u = zeros(L.n)
    K = spatial_curvature_K(c)
    if ic !== :adiabatic
        _isocurvature_ic!(u, c, bg, k, L, x0, mν, ic)
    elseif iszero(K)
        ψ = 1 / (1.5 + 0.4 * f_ν)
        φ = (1 + 0.4 * f_ν) * ψ
        θ = 0.5 * k^2 * η * ψ
        δ_ν = -2ψ
        σ_ν = (1 / 15) * k^2 * η^2 * ψ

        u[1] = -1.5ψ; u[2] = θ
        u[3] = -1.5ψ; u[4] = θ
        u[5] = φ
        u[L.iγ] = -2ψ
        u[L.iγ+1] = (4 / (3k)) * θ
        u[L.iν] = δ_ν
        u[L.iν+1] = (4 / (3k)) * θ
        u[L.iν+2] = 2 * σ_ν
    else
        # Regular adiabatic solution in curved FRW, derived in synchronous
        # gauge and transformed to conformal Newtonian gauge.  These are the
        # leading non-vanishing terms in k*eta; moving x0 earlier is the
        # controlled convergence parameter, just as in CLASS.
        s2 = _curvature_streaming(K, k, 2)
        s2sq = s2^2
        kt2 = (k * η)^2
        kt3 = (k * η)^3
        fg = 1 - f_ν
        δγs = -kt2 * s2sq / 3
        θγs = -k * kt3 * s2sq / 36
        θνs = -k * kt3 * (4f_ν + 11 + 12s2sq) * s2sq /
            (36 * (4f_ν + 15))
        σν = kt2 * (3s2sq - 1) / (45 + 12f_ν)
        Fν3 = 2kt3 / (7 * (12f_ν + 45))
        ηsyn = 1 - kt2 * (5 + 4s2sq * f_ν) / (12 * (15 + 4f_ν))
        Hconf = bg.ℋ(x0)
        δtot = δγs
        θtot = (4 / 3) * (fg * θγs + f_ν * θνs)
        α = (ηsyn + 1.5 * Hconf^2 / (k^2 * s2sq) *
             (δtot + 3Hconf * θtot / k^2)) / Hconf

        δγ = δγs - 4Hconf * α
        θγ = θγs + k^2 * α
        δν = δγs - 4Hconf * α
        θν = θνs + k^2 * α
        δm = 0.75δγs - 3Hconf * α
        θm = k^2 * α
        u[1] = δm; u[2] = θm
        u[3] = δm; u[4] = θγ
        u[5] = ηsyn - Hconf * α
        u[L.iγ] = δγ
        u[L.iγ+1] = 4θγ / (3k)
        u[L.iν] = δν
        u[L.iν+1] = 4θν / (3k)
        u[L.iν+2] = 2s2 * σν
        L.lmax_ν >= 3 && (u[L.iν+3] = Fν3)
    end

    δ_ν = u[L.iν]
    θ_ν = 3k * u[L.iν+1] / 4
    σ_ν = u[L.iν+2] / (2 * _curvature_streaming(K, k, 2))

    # Massive neutrinos start relativistic and adiabatic, so their distribution
    # perturbation is the massless answer reshaped by the Fermi-Dirac weight:
    #   Ψ_ℓ = -(coefficient) · dln f₀/dln q
    for (s, m) in enumerate(mν)
        am = m.y1 * a
        for iq in 1:L.nq
            q = G.q[iq]
            ε = sqrt(q^2 + am^2)
            dlnf0 = G.dlnf0[iq]
            u[_imν(L, s, iq, 0)] = -0.25 * δ_ν * dlnf0
            u[_imν(L, s, iq, 1)] = -(ε / (3 * q * k)) * θ_ν * dlnf0
            u[_imν(L, s, iq, 2)] = -0.5 *
                _curvature_streaming(K, k, 2) * σ_ν * dlnf0
            if L.lmax_m >= 3
                u[_imν(L, s, iq, 3)] = -0.25 * u[L.iν+3] * dlnf0
            end
        end
    end

    # Decaying DM starts as ordinary cold matter (Γ/ℋ ≪ 1 this early, for any Γ
    # expressed in km/s/Mpc), so it tracks the cdm slots; its dark radiation
    # tracks the massless neutrinos, rescaled to the density-weighted moments by
    # f_dr = ρ_dr a⁴ — a minuscule number this early, and the sourced attractor
    # takes over immediately regardless.
    if L.has_dcdm
        dc = _dcdm_params(c)
        f = exp(dc.lnρr(x0)) * a^4
        u[L.idc] = u[1]
        u[L.idc+1] = u[2]
        u[L.idr] = u[L.iν] * f
        u[L.idr+1] = u[L.iν+1] * f
        u[L.idr+2] = u[L.iν+2] * f
    end

    if L.has_fld && !fld_ppf
        # adiabatic relation δ_i/(1+w_i) common to all species; DE is utterly
        # subdominant this early, so the attractor erases any O((kη)²) detail.
        # Velocity slot carries Q = (1+w)θ. (PPF instead starts Γ = 0, i.e. the
        # slots stay at their zero default — same as CLASS.)
        de = _fld_species(c)
        wi = w(de, a)
        u[L.ifld] = (1 + wi) * u[1]
        u[L.ifld+1] = (1 + wi) * u[2]
    end

    u
end

# --- Driver -----------------------------------------------------------------

"""
    solve_perturbations(c, bg, k; lmax_γ = 12, lmax_ν = 12)

Integrate the Boltzmann hierarchy for one mode `k` (1/Mpc).

Defaults for the hierarchy truncations are set by convergence against CAMB, not by
tradition: ℓ_max_γ = 12 (a common textbook choice) leaves a 2% error in P(k) at
k = 1/Mpc, because power reflected off the top of the truncated hierarchy leaks back
into the low multipoles that source the metric. ℓ_max_γ = 25, ℓ_max_ν = 32 brings
P(k) to 0.04% of CAMB across 0.1 < k < 2.

The start time is chosen per mode so that the mode is comfortably outside the
horizon (kη ≈ 1e-2) and deep in radiation domination -- both are what make the
adiabatic initial conditions valid. A fixed early start for every k would waste
effort on large scales and be outright wrong for small ones.

Before recombination the system is stiff: κ̇ exceeds k and ℋ by orders of
magnitude, which is exactly what "tight coupling" means. Boltzmann codes usually
swap in an analytic tight-coupling expansion at that point for speed; we instead
pay for an implicit solver, which is slower but has no extra approximation to
get wrong.
"""
function solve_perturbations(c::Cosmology, bg::BackgroundCache, k::Real;
    lmax_γ=25, lmax_ν=32, lmax_m=12, nq=20,
    a_end=1.0, reltol=1e-8, abstol=1e-10, ic::Symbol=:adiabatic, cs2_de=1.0,
    de_scheme::Symbol=:auto)

    K = spatial_curvature_K(c)
    mν = Tuple(get_all_species(c, MassiveNeutrinos))
    nmν = length(mν)
    # The free-streaming grid samples f₀; use the neutrinos' own distribution
    # (degenerate or tabulated), which is exact when they share one f₀.
    ν_dist = isempty(mν) ? FermiDirac() : first(mν).dist
    G = MassiveNuGrid(nq; dist=ν_dist)
    dcdm = _dcdm_params(c)
    de = _fld_species(c)
    fld = if de === nothing
        nothing
    else
        ppf = if de_scheme === :ppf
            true
        elseif de_scheme === :fluid
            false
        else
            # :auto — the exact (δ,Q) fluid unless w(a) crosses −1 somewhere,
            # where the single-fluid description itself breaks down and the
            # PPF construction takes over (CLASS behaviour).
            _crosses_w_minus1(de)
        end
        (; de, cs2=float(cs2_de), ppf)
    end
    L = Layout(lmax_γ, lmax_ν; lmax_m, nq, nmν, has_dcdm=dcdm !== nothing,
        has_fld=fld !== nothing)
    k = float(k)

    # The initial conditions demand *both* that the mode is outside the horizon
    # (kη ≪ 1) and that we are deep in radiation domination -- the adiabatic
    # growing mode is derived under both assumptions. So take whichever condition
    # bites first.
    #
    # Keying on the horizon alone is wrong for small k: a mode with k = 1e-4 does
    # not reach kη = 1e-2 until a ≈ 2e-4, which is already past equality, and RD
    # initial conditions imposed there are simply false.
    x_horizon = log(1e-9)
    for x in range(log(1e-9), log(1e-2); length=600)
        k * bg.η(x) > 1e-2 && break
        x_horizon = x
    end
    x0 = min(x_horizon, log(1e-5))     # 1e-5 ≪ a_eq ≈ 3e-4
    x0 = max(x0, log(1e-10))

    u0 = initial_conditions(c, bg, k, L, x0, G, mν; ic,
        fld_ppf=fld !== nothing && fld.ppf)
    prob = ODEProblem(_rhs!, u0, (x0, log(a_end)), (bg, k, K, L, G, mν, dcdm, fld))
    sol = solve(prob, KenCarp4(); reltol, abstol)

    PerturbationSolution(c, bg.rec, k, sol, lmax_γ, lmax_ν, lmax_m, nq, nmν)
end

"Convenience: build the background cache and solve in one call."
function solve_perturbations(c::Cosmology, rec::RecombinationSolution, k::Real; kws...)
    solve_perturbations(c, BackgroundCache(c, rec), k; kws...)
end

# --- Accessors --------------------------------------------------------------

δ_cdm(p::PerturbationSolution, a) = p.sol(log(a))[1]
θ_cdm(p::PerturbationSolution, a) = p.sol(log(a))[2]
δ_baryon(p::PerturbationSolution, a) = p.sol(log(a))[3]
θ_baryon(p::PerturbationSolution, a) = p.sol(log(a))[4]
δ_photon(p::PerturbationSolution, a) = p.sol(log(a))[_layout(p).iγ]
δ_neutrino(p::PerturbationSolution, a) = p.sol(log(a))[_layout(p).iν]

"""
    Φ(p, a)

Spatial curvature perturbation. Constant outside the horizon in both radiation
and matter domination, and famously drops by a factor 9/10 across the transition
between them.
"""
Φ(p::PerturbationSolution, a) = p.sol(log(a))[5]

"""
    Ψ(p, a)

Newtonian potential. Equal to Φ up to the anisotropic stress, which is carried by
free-streaming neutrinos and, after decoupling, photons.
"""
function Ψ(p::PerturbationSolution, a)
    c = p.cosmo
    L = _layout(p)
    u = p.sol(log(a))
    H0 = Constants.H0_in_invMpc(c.h)
    ℋ_ = ℋ(c, a)
    ργ = Ω_γ(c) / a^4
    ρν = _Ω_or_zero(c, MasslessNeutrinos) / a^4
    ρb = Ω_b(c) / a^3
    ρc = Ω_c(c) / a^3
    G = MassiveNuGrid(L.nq == 0 ? 1 : L.nq)
    mν = Tuple(get_all_species(c, MassiveNeutrinos))
    # Reuse the metric routine so the massive-neutrino anisotropic stress is
    # included consistently rather than duplicated here.
    dc = _dcdm_params(c)
    ρdc = dc === nothing ? 0.0 : exp(dc.lnρd(log(a)))
    ρdr = dc === nothing ? 0.0 : exp(dc.lnρr(log(a)))
    defl = _fld_species(c)
    ρfld = defl === nothing ? 0.0 : ρ_over_ρc0(defl, a)
    wfld = defl === nothing ? -1.0 : w(defl, a)
    _, ψ, _ = _metric(u, L, p.k, spatial_curvature_K(c), a, ℋ_, H0,
        ργ, ρν, ρb, ρc, G, mν, ρdc, ρdr, ρfld, wfld)
    ψ
end

"Photon temperature multipole Θ_ℓ = F_γℓ/4."
Θ(p::PerturbationSolution, ℓ::Int, a) = p.sol(log(a))[_layout(p).iγ+ℓ] / 4

"""
    δ_matter(p, a)

Total matter density contrast in the **conformal Newtonian gauge**, ρ-weighted over
CDM and baryons. Baryons are ~16% of the matter and do not trace CDM inside the
sound horizon, so this is not δ_cdm.

For the matter power spectrum you almost certainly want
[`δ_matter_comoving`](@ref) instead -- see the warning there.
"""
function δ_matter(p::PerturbationSolution, a)
    c = p.cosmo
    u = p.sol(log(a))
    ρb, ρc = Ω_b(c), Ω_c(c)
    (ρc * u[1] + ρb * u[3]) / (ρc + ρb)
end

"""
    δ_matter_comoving(p, a)

Total matter density contrast in the **comoving (synchronous) gauge**,

    δ_sync = δ_Newtonian + 3ℋ θ_c / k²

This -- not the Newtonian-gauge δ -- is what "the matter power spectrum" means by
convention, and what CAMB and CLASS report.

The distinction is not academic. A density contrast is gauge dependent, and the
two definitions agree only well inside the horizon, where the correction term
3ℋθ/k² is negligible. As k approaches ℋ they diverge badly: at k = 5e-4 /Mpc the
Newtonian-gauge P(k) is **70% larger** than the comoving one. Using the wrong one
looks fine at k ≳ 0.01 and is quietly, badly wrong on the largest scales -- the
ones surveys are increasingly trying to measure.

The synchronous gauge is defined by CDM being at rest, which is why the correction
is built from θ_c specifically.
"""
function δ_matter_comoving(p::PerturbationSolution, a)
    c = p.cosmo
    u = p.sol(log(a))
    ρb, ρc = Ω_b(c), Ω_c(c)
    δ_N = (ρc * u[1] + ρb * u[3]) / (ρc + ρb)
    θ_c = u[2]
    δ_N + 3 * ℋ(c, a) * θ_c / p.k^2
end
