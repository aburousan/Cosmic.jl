"""
Primordial black holes: the mass fraction of the universe collapsing into black
holes when enhanced curvature perturbations re-enter the horizon during
radiation domination, and the resulting present-day abundance and mass function.

The calculation follows Young, Musco & Byrnes (arXiv:1904.00984) — the
treatment that keeps the pieces the abundance is exponentially sensitive to:

  * the **non-linear relation** between the (Gaussian) curvature perturbation
    and the top-hat-smoothed density contrast, δ_m = δ_l − (3/8)δ_l² — using
    the linear relation overstates the abundance by orders of magnitude;
  * **peaks theory** (rather than Press–Schechter) for the number density of
    rare peaks;
  * **critical collapse**, M = 𝒦 M_H (δ_m − δ_c)^γ with γ ≃ 0.36, so each
    horizon mass produces a spectrum of black-hole masses;
  * the simulation-calibrated **threshold** δ_c (0.41–2/3 depending on the
    density profile; 0.55 is the representative value for smooth peaks) with
    the type-I branch δ_l < 4/3 enforced exactly.

The linear component of the smoothed density contrast has variance

    σ²(r_m) = (16/81) ∫ dk/k (k r_m)⁴ W̃²(k r_m) T²(k r_m) 𝒫_ζ(k),

with W̃ the top-hat window and T the linear radiation-era transfer function,
and the peak-scale moment μ² carries an extra k². The collapsing fraction at
one horizon scale is the ν-integral of eq. (30) of the paper, and the total
present-day abundance integrates over horizon masses with the (M_eq/M_H)^{1/2}
growth weight. Everything takes an arbitrary 𝒫_ζ callable — including an
`InflatonSpectrum` — so the same primordial object feeds the CMB, the induced
gravitational waves, and the black holes.

Black holes lighter than ~5×10¹⁴ g have evaporated by today; they are excluded
from `f_pbh` (never report a present-day dark-matter fraction for them).
"""

using QuadGK: quadgk

# top-hat window and linear RD transfer function (paper eqs. 21-22)
@inline _pbh_W(x) = 3 * (sin(x) - x * cos(x)) / x^3
@inline _pbh_T(x) = (y = x / sqrt(3); 3 * (sin(y) - y * cos(y)) / y^3)
# the sinc window that shows up when the curvature field itself (not its
# gradient) is volume-averaged — the ζ-average in the compaction variances
@inline _pbh_Ws(x) = x == 0 ? 1.0 : sin(x) / x

"""
    pbh_moments(Pζ, r_m; rtol = 1e-6) -> (σ², μ²)

The variance σ² and peak-scale moment μ² of the linear top-hat-smoothed density
contrast at horizon scale `r_m` (Mpc), for dimensionless primordial spectrum
`Pζ(k)` (paper eqs. 20 and 28). The integrand support is k r_m ~ O(1); the
window+transfer product cuts it off well before k r_m ~ 40.
"""
function pbh_moments(Pζ, r_m; rtol=1e-6)
    # Integrate in ln k over unit-width segments: primordial spectra in the PBH
    # context are often sharply peaked (USR features), and a wide adaptive
    # interval can straddle such a peak with its first quadrature nodes and
    # silently miss it. Unit segments in ln k bound the structure any segment
    # can hide.
    lo, hi = log(1e-4 / r_m), log(40.0 / r_m)
    segs = collect(range(lo, hi; length=ceil(Int, hi - lo) + 1))
    fσ(u) = (k = exp(u); x = k * r_m; x^4 * _pbh_W(x)^2 * _pbh_T(x)^2 * Pζ(k))
    fμ(u) = (k = exp(u); fσ(u) * k^2)
    σ2 = (16 / 81) * quadgk(fσ, segs...; rtol)[1]
    μ2 = (16 / 81) * quadgk(fμ, segs...; rtol)[1]
    (σ2, μ2)
end

"""
    pbh_β(Pζ, r_m; δc = 0.55, K = 4.0, γ = 0.36, rtol = 1e-6)

The mass fraction of the universe collapsing to black holes at horizon scale
`r_m`, weighted by the critical-collapse mass ratio (paper eq. 30):

    β = ∫ dν (𝒦/3π) (νσ − (3/8)(νσ)² − δ_c)^γ (μ r_m/σ)³ ν³ e^{−ν²/2}

between ν_c = δ_{c,l−}/σ (the critical linear amplitude, eq. 23) and 4/(3σ)
(the type-I upper bound). Returns 0 when the threshold is unreachable.
"""
function pbh_β(Pζ, r_m; δc=0.55, K=4.0, γ=0.36, rtol=1e-6)
    σ2, μ2 = pbh_moments(Pζ, r_m; rtol)
    σ = sqrt(σ2)
    # critical linear amplitude, type-I branch (eq. 23)
    δc < 2 / 3 || throw(ArgumentError("δc must be < 2/3 (the type-I maximum)"))
    δcl = (4 / 3) * (1 - sqrt((2 - 3 * δc) / 2))
    νc = δcl / σ
    νmax = 4 / (3σ)
    νc ≥ νmax && return 0.0
    # (μ/(aH σ))³ with aH = 1/r_m
    peakcube = (sqrt(μ2) * r_m / σ)^3
    integrand(ν) = begin
        δm = ν * σ - (3 / 8) * (ν * σ)^2
        δm ≤ δc && return 0.0
        (K / (3π)) * (δm - δc)^γ * peakcube * ν^3 * exp(-ν^2 / 2)
    end
    quadgk(integrand, νc, νmax; rtol)[1]
end

# --- local non-Gaussianity: threshold statistics on the compaction function --
#
# Everything above assumes a Gaussian curvature field. The abundance is
# exponentially sensitive to the tail of that field, so any primordial
# non-Gaussianity — a curvaton, a USR feature — moves β by orders of magnitude
# and cannot be neglected for a serious PBH prediction. We follow Ferrante,
# Franciolini, Iovino & Urbano (arXiv:2211.01728), who compute β exactly for a
# generic local relation ζ = F(ζ_G) by working with the compaction function
# rather than the density contrast directly.
#
# The linear compaction is 𝒞₁ = 𝒞_G · F'(ζ_G) and the full (non-linear)
# compaction is 𝒞 = 𝒞₁ − (1/4Φ)𝒞₁² (their eqs. 44, 47, 48), with Φ = 2/3 in
# radiation so that 1/4Φ = 3/8 — the same non-linear map already used above,
# now carrying the intrinsic NG through F'. Both 𝒞_G and ζ_G are Gaussian and
# correlated; the abundance is the two-dimensional tail integral (eq. 56)
#
#     β = ∫_𝒟 𝒦(𝒞 − 𝒞_th)^γ P_G(𝒞_G, ζ_G) d𝒞_G dζ_G,
#
# over the type-I domain 𝒟 = {𝒞 > 𝒞_th ∧ 𝒞₁ < 2Φ}. When F(ζ_G) = ζ_G the map
# F' ≡ 1 and this collapses to the Gaussian compaction-threshold result.

# bivariate PDF of the linear compaction 𝒞_G and the curvature ζ_G, eq. (57):
# a correlated Gaussian with variances σ_c², σ_r² and correlation γ_cr.
@inline function _pbh_pg(CG, ζG, σc, σr, γcr)
    z = CG / σc - γcr * ζG / σr
    inv(2π * σc * σr * sqrt(1 - γcr^2)) *
    exp(-ζG^2 / (2σr^2)) * exp(-z^2 / (2 * (1 - γcr^2)))
end

"""
    pbh_compaction_variances(Pζ, r_m; Φ = 2/3, rtol = 1e-6) -> (σ_c², σ_cr², σ_r²)

The three second moments of the compaction/curvature pair at horizon scale
`r_m` (Mpc), for dimensionless primordial spectrum `Pζ(k)` — Ferrante et al.
eqs. (50)–(52):

    σ_c²  = (4Φ²/9) ∫ dk/k (k r_m)⁴ W² T² 𝒫_ζ,
    σ_cr² = (2Φ/3)  ∫ dk/k (k r_m)² W Wₛ T² 𝒫_ζ,
    σ_r²  =         ∫ dk/k          Wₛ² T² 𝒫_ζ,

with `W` the top-hat window, `Wₛ` the sinc (volume-average) window, `T` the RD
transfer function, and `Φ = 2/3` in radiation. `σ_c²` is by construction the
same integral as the σ² returned by [`pbh_moments`](@ref).
"""
function pbh_compaction_variances(Pζ, r_m; Φ=2 / 3, rtol=1e-6)
    lo, hi = log(1e-4 / r_m), log(40.0 / r_m)
    segs = collect(range(lo, hi; length=ceil(Int, hi - lo) + 1))
    fc(u) = (k = exp(u); x = k * r_m; x^4 * _pbh_W(x)^2 * _pbh_T(x)^2 * Pζ(k))
    fcr(u) = (k = exp(u); x = k * r_m; x^2 * _pbh_W(x) * _pbh_Ws(x) * _pbh_T(x)^2 * Pζ(k))
    fr(u) = (k = exp(u); x = k * r_m; _pbh_Ws(x)^2 * _pbh_T(x)^2 * Pζ(k))
    σc2 = (4Φ^2 / 9) * quadgk(fc, segs...; rtol)[1]
    σcr2 = (2Φ / 3) * quadgk(fcr, segs...; rtol)[1]
    σr2 = quadgk(fr, segs...; rtol)[1]
    (σc2, σcr2, σr2)
end

# derivative F'(ζ_G) of the local map for the power-series parametrisation
# ζ = ζ_G + (3/5)f_NL ζ_G² + (9/25)g_NL ζ_G³ (eq. 3). The constant −σ² sometimes
# written into the f_NL term drops out of the derivative, so it never enters here.
@inline _pbh_dF(ζG, f_NL, g_NL) = 1 + (6 / 5) * f_NL * ζG + (27 / 25) * g_NL * ζG^2

"""
    pbh_β_ng(Pζ, r_m; f_NL = 0, g_NL = 0, dFdζ = nothing, δc = 0.55, K = 4.0,
             γ = 0.36, Φ = 2/3, ζmax = 12.0, rtol = 1e-6)

The PBH mass fraction at horizon scale `r_m` (Mpc) for a spectrum `Pζ(k)` with
local non-Gaussianity, via threshold statistics on the compaction function
(Ferrante et al. eq. 56). The curvature map ζ = F(ζ_G) enters only through its
derivative: pass a callable `dFdζ(ζ_G)` for the exact form (e.g. the curvaton
`F` of their eq. 4), or leave it `nothing` to use the `f_NL`/`g_NL` power
series. `δc` is the compaction threshold 𝒞_th (< Φ = 2/3), `ζmax` the ζ_G
integration half-width in units of σ_r.

Reduces to the Gaussian compaction-threshold abundance when `f_NL = g_NL = 0`
(or `dFdζ ≡ 1`). Returns 0 when the threshold cannot be reached.

The threshold `δc` (= 𝒞_th) is held at its Gaussian, simulation-calibrated value,
following Ferrante et al., who show the *statistics* dominate the NG response.
The threshold itself also drifts under NG — a few percent for |f_NL| ≲ 𝒪(5)
(Kehagias, Musco & Riotto, arXiv:1906.07135), set by the sim-calibrated
shape→threshold relation. That drift is not modelled here; pass an NG-appropriate
`δc` if you want to fold it in.
"""
function pbh_β_ng(Pζ, r_m; f_NL=0.0, g_NL=0.0, dFdζ=nothing, δc=0.55, K=4.0,
    γ=0.36, Φ=2 / 3, ζmax=12.0, rtol=1e-6)
    δc < Φ || throw(ArgumentError("δc (= 𝒞_th) must be < Φ = $(Φ), the maximum compaction"))
    σc2, σcr2, σr2 = pbh_compaction_variances(Pζ, r_m; Φ, rtol)
    (σc2 ≤ 0 || σr2 ≤ 0) && return 0.0
    σc, σr = sqrt(σc2), sqrt(σr2)
    γcr = σcr2 / (σc * σr)
    dF(ζG) = dFdζ === nothing ? _pbh_dF(ζG, f_NL, g_NL) : dFdζ(ζG)

    # the type-I 𝒞₁ window: from the lower threshold root up to the parabola's
    # apex at 2Φ (eqs. 44, 58). Integrating in 𝒞₁ keeps the limits fixed for
    # every ζ_G — the 1/|F'| Jacobian and P_G(𝒞₁/F', ζ_G) carry the NG.
    C1lo = 2Φ * (1 - sqrt(1 - δc / Φ))
    C1hi = 2Φ
    Φ4 = 1 / (4Φ)
    inner(ζG) = begin
        d = dF(ζG)
        d == 0 && return 0.0
        quadgk(C1 -> begin
                C = C1 - Φ4 * C1^2
                C ≤ δc && return 0.0
                K * (C - δc)^γ / abs(d) * _pbh_pg(C1 / d, ζG, σc, σr, γcr)
            end, C1lo, C1hi; rtol)[1]
    end
    quadgk(inner, -ζmax * σr, ζmax * σr; rtol)[1]
end

const _M_sun_g = 1.98841e33
const _M_evap_g = 5e14                 # evaporated by today below this

"""
    PBHAbundance

Present-day primordial-black-hole abundance: horizon wavenumbers `k` (1/Mpc),
horizon masses `M_H` and black-hole masses `M` (solar masses), the formation
fraction `β(M_H)`, the mass function `f(M) = dΩ_PBH/dln M / Ω_cdm`, and the
total `f_pbh = Ω_PBH/Ω_cdm` (counting only non-evaporated masses).
"""
struct PBHAbundance
    k::Vector{Float64}
    M_H::Vector{Float64}
    β::Vector{Float64}
    M::Vector{Float64}
    f::Vector{Float64}
    f_pbh::Float64
    f_evaporated::Float64   # formation-weighted fraction below the evaporation floor
end

"""
    pbh_abundance(Pζ, c::Cosmology; kmin, kmax, nk = 80, δc = 0.55, K = 4.0,
                  γ = 0.36)

Present-day PBH abundance for primordial spectrum `Pζ` over horizon scales
`kmin..kmax` (1/Mpc). Horizon masses are anchored to the cosmology's own
matter–radiation equality, M_H(k) = M_eq (k/k_eq)⁻² with
M_eq = 2.8×10¹⁷ M_⊙ (paper eq. 25/31 convention), and the present abundance is

    Ω_PBH = ∫ dln M_H (M_eq/M_H)^{1/2} β(M_H).

The mass function is reported against the *black-hole* mass
M = 𝒦 M_H (δ_peak − δ_c)^γ evaluated at the β-weighted mean overdensity.
Masses below the evaporation floor (5×10¹⁴ g ≈ 2.5×10⁻¹⁹ M_⊙) never enter
`f_pbh`; their formation fraction is reported separately as `f_evaporated`.

Pass `f_NL`, `g_NL`, or a curvature map derivative `dFdζ` to include primordial
local non-Gaussianity; the formation fractions then come from [`pbh_β_ng`](@ref)
(compaction-threshold statistics) instead of the Gaussian peaks-theory `pbh_β`.
"""
function pbh_abundance(Pζ, c::Cosmology; kmin, kmax, nk=80, δc=0.55, K=4.0,
    γ=0.36, f_NL=0.0, g_NL=0.0, dFdζ=nothing, rtol=1e-6)
    k_eq = _k_equality(c)
    M_eq = 2.8e17
    ks = exp.(range(log(kmin), log(kmax); length=nk))
    MH = M_eq .* (ks ./ k_eq) .^ -2
    # any primordial NG routes through the compaction-threshold integral; the
    # Gaussian case keeps the lighter peaks-theory β
    ng = f_NL != 0 || g_NL != 0 || dFdζ !== nothing
    βs = ng ?
         [pbh_β_ng(Pζ, 1 / k; f_NL, g_NL, dFdζ, δc, K, γ, rtol) for k in ks] :
         [pbh_β(Pζ, 1 / k; δc, K, γ, rtol) for k in ks]

    # representative black-hole mass per horizon scale: critical collapse at the
    # typical forming amplitude (δ_m one σ above threshold, capped at type-I max)
    Ms = similar(MH)
    for (i, k) in enumerate(ks)
        σ = sqrt(pbh_moments(Pζ, 1 / k; rtol)[1])
        δtyp = min(δc + σ, 2 / 3)
        Ms[i] = K * MH[i] * (δtyp - δc)^γ
    end

    # Ω_PBH today (eq. 31), splitting evaporated from surviving masses
    M_floor = _M_evap_g / _M_sun_g
    Ω = 0.0
    Ωev = 0.0
    lnMH = log.(MH)
    for i in 1:(nk-1)
        w = 0.5 * (sqrt(M_eq / MH[i]) * βs[i] + sqrt(M_eq / MH[i+1]) * βs[i+1]) *
            abs(lnMH[i+1] - lnMH[i])
        if Ms[i] > M_floor
            Ω += w
        else
            Ωev += w
        end
    end
    Ωcdm = Ω_c(c)
    f = [sqrt(M_eq / MH[i]) * βs[i] / Ωcdm for i in 1:nk]
    PBHAbundance(ks, MH, βs, Ms, f, Ω / Ωcdm, Ωev / Ωcdm)
end

# k at matter-radiation equality from the cosmology's own background
function _k_equality(c::Cosmology)
    zeq = z_equality(c)
    aeq = 1 / (1 + zeq)
    aeq * H_Mpc(c, aeq)
end