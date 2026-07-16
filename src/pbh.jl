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
"""
function pbh_abundance(Pζ, c::Cosmology; kmin, kmax, nk=80, δc=0.55, K=4.0,
    γ=0.36, rtol=1e-6)
    k_eq = _k_equality(c)
    M_eq = 2.8e17
    ks = exp.(range(log(kmin), log(kmax); length=nk))
    MH = M_eq .* (ks ./ k_eq) .^ -2
    βs = [pbh_β(Pζ, 1 / k; δc, K, γ, rtol) for k in ks]

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