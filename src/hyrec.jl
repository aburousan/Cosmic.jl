# Hydrogen recombination at HyRec-2020 accuracy (SWIFT mode).
#
# RECFAST is a three-level atom with a fudge factor: the multilevel cascade, the
# l-resolved high-n states, two-photon decays from n ≥ 3 and Lyman-α feedback are all
# bundled into F = 1.125 and a Gaussian correction to K. HyRec-2020 (Ali-Haïmoud &
# Hirata 2011; SWIFT mode: Lee & Ali-Haïmoud 2020, arXiv:2007.14114) does it properly:
#
#  * An *effective four-level atom* (1s, 2s, 2p, continuum). Everything above n = 2 is
#    folded into effective recombination coefficients α_2s(T_M, T_R), α_2p(T_M, T_R) and
#    an effective 2p ↔ 2s transfer rate, computed once from the full l-resolved
#    multilevel atom (n up to 500) and tabulated. That is `Alpha_inf.dat`/`R_inf.dat`.
#  * The radiative-transfer effects near Lyman-α (two-photon decays from higher n,
#    Raman scattering, diffusion, feedback) enter as a single correction ΔK/K to the
#    Lyman-α escape rate, fitted from the full radiative-transfer solve and expanded to
#    first order around a fiducial cosmology in (ω_cb, ω_b(1−Y_He), N_eff). That is
#    `fit_swift.dat`, valid for 1775 K < T_R < 4415 K -- exactly the window where
#    radiative transfer matters. Outside it ΔK/K = 0 and the model is the effective
#    multilevel atom, which is the correct limit on both sides.
#
# The photoionization rates are not tabulated: they follow from α by detailed balance,
#    β_l = α_l(T_R, T_R) · (2πμ_e T_R/h²)^{3/2} e^{−E_2/T_R} / (2l+1),
# which is the same discipline as the BBN reverse rates -- forward rates are physics,
# reverse rates are thermodynamics.
#
# Everything here works in HyRec's native units (temperatures in eV, densities in
# cm⁻³) and mirrors its 4-point Lagrange interpolation exactly, so the two
# implementations can be compared digit by digit.

const _HYREC_DIR = joinpath(@__DIR__, "..", "data", "hyrec")

# hydrogen.h constants (fsR = meR = 1: no varying fundamental constants)
const _HY_EI = 13.598286071938324      # H ionization energy, eV (reduced mass)
const _HY_L2s1s = 8.2206               # 2s -> 1s two-photon rate, 1/s
const _HY_kB = 8.617343e-5             # eV/K
const _HY_SAHA = 3.016103031869581e21  # (2π μ_e/h²)^{3/2} in eV^{-3/2} cm^{-3}
const _HY_LYA = 4.662899067555897e15   # 8π/(3λ_Lyα³) in cm^{-3}
const _HY_NTR = 100
const _HY_NTM = 40
const _HY_lnTR_MIN, _HY_lnTR_MAX = log(0.004), log(0.4)
const _HY_RATIO_MIN = 0.1

struct HyrecTables
    logα::Array{Float64,3}      # (4, NTM, NTR)
    logR2p2s::Vector{Float64}   # (NTR)
    dlnTR::Float64
    dratio::Float64
    swift_T::Float64            # first T of the swift grid (K)
    swift_dT::Float64           # grid step (K)
    swift::Matrix{Float64}      # (265, 4): fiducial + 3 derivatives
end

const _HYREC = Ref{Union{Nothing,HyrecTables}}(nothing)

function _hyrec_tables()
    if _HYREC[] === nothing
        A = readdlm(joinpath(_HYREC_DIR, "Alpha_inf.dat"))
        logα = Array{Float64}(undef, 4, _HY_NTM, _HY_NTR)
        for i in 1:_HY_NTR, j in 1:_HY_NTM, l in 1:4
            logα[l, j, i] = log(A[(i-1)*_HY_NTM+j, l])
        end
        R = vec(readdlm(joinpath(_HYREC_DIR, "R_inf.dat")))
        S = readdlm(joinpath(_HYREC_DIR, "fit_swift.dat"))
        _HYREC[] = HyrecTables(logα, log.(R),
            (_HY_lnTR_MAX - _HY_lnTR_MIN) / (_HY_NTR - 1),
            (1.0 - _HY_RATIO_MIN) / (_HY_NTM - 1),
            S[1, 1], S[2, 1] - S[1, 1], S[:, 2:5])
    end
    _HYREC[]
end

# HyRec's 4-point Lagrange weights on a uniform grid: i0 is the second of the four
# points, frac ∈ [0,1) the offset past it. Matches rec_interp1d / interpolate_rates.
@inline function _lag4(frac)
    (frac * (frac - 1) * (2 - frac) / 6,
        (1 + frac) * (1 - frac) * (2 - frac) / 2,
        (1 + frac) * frac * (2 - frac) / 2,
        (1 + frac) * frac * (frac - 1) / 6)
end

@inline function _locate(x, x0, dx, n)
    i = floor(Int, (x - x0) / dx)
    i = clamp(i, 1, n - 3)
    (i, (x - x0) / dx - i)
end

"""
Effective rates at radiation temperature `TR` and ratio `TM/TR` (both eV-based):
returns (α_2s, α_2p, Δα_2s, Δα_2p, β_2s, β_2p, R_2p2s) in cm³/s and 1/s.
"""
function _hyrec_rates(TR, TM_over_TR)
    t = _hyrec_tables()
    # T_RATIO is min(TM/TR, TR/TM); the table's l-offset picks the branch.
    ratio = TM_over_TR > 1 ? 1 / TM_over_TR : TM_over_TR
    off = TM_over_TR > 1 ? 2 : 0
    ratio = clamp(ratio, _HY_RATIO_MIN, 1.0)

    iR, fR = _locate(log(TR), _HY_lnTR_MIN, t.dlnTR, _HY_NTR)
    w2 = _lag4(fR)
    iM, fM = _locate(ratio, _HY_RATIO_MIN, t.dratio, _HY_NTM)
    w1 = _lag4(fM)

    # Scalar per level, no Float64-typed temporaries: the stiff solver differentiates
    # this whole function with ForwardDiff, and a zeros(2) here would reject the Duals.
    function level(l)
        αeq = exp(sum(t.logα[l, _HY_NTM, iR-1+k] * w2[k] for k in 1:4))
        β = αeq * _HY_SAHA * TR * sqrt(TR) * exp(-0.25 * _HY_EI / TR) / (2l - 1)
        tmp = ntuple(k -> sum(t.logα[l+off, iM-1+k, iR-1+m] * w2[m] for m in 1:4), 4)
        α = exp(sum(tmp[k] * w1[k] for k in 1:4))
        (α, α - αeq, β)
    end
    α1, Δα1, β1 = level(1)
    α2, Δα2, β2 = level(2)
    R2p2s = exp(sum(t.logR2p2s[iR-1+k] * w2[k] for k in 1:4))
    (α1, α2, Δα1, Δα2, β1, β2, R2p2s)
end

"""
SWIFT radiative-transfer correction ΔK/K at radiation temperature `TR_K` (Kelvin), for
a cosmology offset from the fit's fiducial by `dω_cb`, `dω_bH` = Δ[ω_b(1−Y_He)] and
`dNeff`, all pre-scaled by (T_fid/T₀)³ as in HyRec.
"""
function _swift_DKK(TR_K, dω_cb, dω_bH, dNeff)
    t = _hyrec_tables()
    n = size(t.swift, 1)
    (TR_K <= t.swift_T || TR_K >= t.swift_T + t.swift_dT * (n - 1)) && return 0.0
    i, f = _locate(TR_K, t.swift_T, t.swift_dT, n)
    w = _lag4(f)
    v(col) = sum(t.swift[i-1+k, col] * w[k] for k in 1:4)
    v(1) + dω_cb * v(2) + dω_bH * v(3) + dNeff * v(4)
end

"""
    hyrec_dxHII_dlna(xe, xHII, nH_cm3, H, TM_eV, TR_eV, dω_cb, dω_bH, dNeff)

d x_HII / d ln a from the effective four-level atom with the SWIFT Lyman-α
correction -- HyRec-2020's default hydrogen model. Negative while recombining.

Also returns the Lyman-α non-return fraction ¼(1−C_2s) + ¾(1−C_2p): the share of
2p/2s excitations that end in ionization rather than cascading back to the ground
state, which is what an exotic Lyman-α photon needs to know to ionize anything.
"""
function hyrec_dxHII_dlna(xe, xHII, nH, H, TM, TR, dω_cb, dω_bH, dNeff)
    α2s, α2p, Δα2s, Δα2p, β2s, β2p, R2p2s = _hyrec_rates(TR, TM / TR)

    RLya = _HY_LYA * H / nH / (1 - xHII)
    DKK = _swift_DKK(TR / _HY_kB, dω_cb, dω_bH, dNeff)
    RLya_fit = RLya / (1 + DKK)

    γ2s = β2s + 3 * R2p2s + _HY_L2s1s
    γ2p = β2p + R2p2s + RLya_fit
    C2s = (_HY_L2s1s + 3 * R2p2s * RLya_fit / γ2p) / (γ2s - 3 * R2p2s^2 / γ2p)
    C2p = (RLya_fit + R2p2s * _HY_L2s1s / γ2s) / (γ2p - 3 * R2p2s^2 / γ2s)

    s = _HY_SAHA * TR * sqrt(TR) * exp(-_HY_EI / TR) / nH
    Dxe2 = xe * xHII - s * (1 - xHII)

    dx = -nH / H * ((s * (1 - xHII) * Δα2s + α2s * Dxe2) * C2s +
                    (s * (1 - xHII) * Δα2p + α2p * Dxe2) * C2p)
    (dx, 0.25 * (1 - C2s) + 0.75 * (1 - C2p))
end
