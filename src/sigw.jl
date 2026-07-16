"""
Scalar-induced gravitational waves (SIGW): the stochastic background sourced at
second order by the primordial curvature spectrum, the smoking-gun counterpart
of primordial-black-hole formation. Everything here is the exact semianalytic
Kohri–Terada calculation (arXiv:1804.08577) — their eqs. (6), (17), (27) — with
no fits: the two-dimensional (t,s) integral over the primordial spectrum is done
by adaptive quadrature against any P_ζ(k), including a computed inflaton
spectrum (`InflatonSpectrum`) or the `Cosmology.primordial` callable.

During radiation domination the induced Ω_GW(k) freezes (GWs redshift like the
radiation background), so the observable today is

    Ω_GW,0(k) h² = Ω_GW,RD(k) · Ω_r,0 h² · (g_*(T_k)/g_{*,0})(g_{*s,0}/g_{*s}(T_k))^{4/3}

with the g_* factor an optional user hook (defaults to 1; the standard constant
approximation for T ≫ 100 GeV multiplies by ≈ 0.39 — pass `g_factor` for
precision work). Frequency map: f[Hz] = k[Mpc⁻¹] · c/(2π Mpc) = 1.5457e-15 · k.

The integrated background is radiation and counts toward N_eff:

    ΔN_eff = (8/7)(11/4)^{4/3} · (1/Ω_γ,0) ∫ Ω_GW,0(k) dln k,

which `sigw_ΔNeff` evaluates — wire it against the BBN/CMB bound before
believing any large peak — an integrated background this large is excluded regardless of any detector curve.
"""

using QuadGK: quadgk

# --- the Kohri–Terada RD kernel ----------------------------------------------

# x²·Ī²_RD(t,s,x→∞): eq. (27) multiplied by the x² of Ω_GW = (1/24)(k/ℋ)²P̄_h
# (ℋ = 1/η in RD, so (k/ℋ)² = x² exactly cancels the kernel's 1/x²).
@inline function _kt_x2I2(t, s)
    P = -5 + s^2 + t * (2 + t)                 # = 2(u²+v²−3) in the (t,s) variables
    pref = 288 * P^2 / ((1 - s + t)^6 * (1 + s + t)^6)
    osc = -(t - s + 1) * (t + s + 1) + 0.5 * P * log(abs((-2 + t * (2 + t)) / (3 - s^2)))
    res = t > sqrt(3) - 1 ? (π^2 / 4) * P^2 : 0.0
    pref * (res + osc^2)
end

# [t(2+t)(s²−1)/((1−s+t)(1+s+t))]²  — the projection factor of eq. (17)
@inline function _kt_proj(t, s)
    (t * (2 + t) * (s^2 - 1) / ((1 - s + t) * (1 + s + t)))^2
end

"""
    Ω_gw_rd(Pζ, k; tmax = 500, rtol = 1e-5)

The induced-GW density parameter **during radiation domination** at wavenumber
`k` (1/Mpc), for the dimensionless primordial spectrum `Pζ(k)` (any callable).
Kohri–Terada eqs. (6)+(17)+(27):

    Ω_GW,RD(k) = (1/12) ∫₀^∞ dt ∫₋₁¹ ds  proj(t,s) · x²Ī²(t,s) · P_ζ(ku)P_ζ(kv)

with u = (t+s+1)/2, v = (t−s+1)/2 (the 1/12 = (1/24)·2, the 2 from the (t,s)
Jacobian form). The kernel's log singularity and the resonance step both sit at
t = √3−1, which is passed to the quadrature as a breakpoint.
"""
function Ω_gw_rd(Pζ, k; tmax=500.0, rtol=1e-5)
    inner(t) = quadgk(s -> begin
            u = (t + s + 1) / 2
            v = (t - s + 1) / 2
            _kt_proj(t, s) * _kt_x2I2(t, s) * Pζ(k * u) * Pζ(k * v)
        end, -1.0, 1.0; rtol)[1]
    I = quadgk(inner, 0.0, sqrt(3) - 1, 2.0, 20.0, tmax; rtol)[1]
    I / 12
end

"""
    sigw_monochromatic(k̃)

Kohri–Terada eq. (29): the exact closed-form Ω_GW,RD(k)/A_ζ² for a monochromatic
spectrum P_ζ = A_ζ δ(ln k/k*), as a function of k̃ = k/k*. The log-resonance at
k̃ = 2/√3 and the cutoff at k̃ = 2 are physical. Used as the analytic anchor for
the numerical kernel.
"""
function sigw_monochromatic(k̃)
    (k̃ ≤ 0 || k̃ ≥ 2) && return 0.0
    y = 3 * k̃^2 - 2
    L = 4 + y * log(abs(1 - 4 / (3 * k̃^2)))
    res = 2 * sqrt(3) - 3 * k̃ > 0 ? π^2 * y^2 : 0.0
    (3 / 64) * ((4 - k̃^2) / 4)^2 * k̃^2 * y^2 * (res + L^2)
end

# --- today's spectrum + N_eff gate -------------------------------------------

"""
    SIGWSpectrum

Induced-GW background evaluated today: wavenumbers `k` (1/Mpc), frequencies
`f` (Hz), and `Ωh²` = Ω_GW,0(k)h². Built by [`sigw_spectrum`](@ref).
"""
struct SIGWSpectrum
    k::Vector{Float64}
    f::Vector{Float64}
    Ωh²::Vector{Float64}
end

const _f_of_k_Hz = Constants.c_SI / (2π * Constants.Mpc_SI)   # 1.546e-15 Hz per 1/Mpc

"""
    sigw_spectrum(Pζ, c::Cosmology; kmin, kmax, nk = 60, g_factor = k -> 1.0)

Today's induced-GW spectrum Ω_GW,0(k)h² for primordial spectrum `Pζ` (callable;
pass `k -> primordial_power(c, k)` to use the cosmology's own, or an
[`InflatonSpectrum`](@ref) directly). `g_factor(k)` is the relative g_*
correction (defaults to 1; ≈ 0.39 for modes entering above the electroweak
scale — supply your own for precision).
"""
function sigw_spectrum(Pζ, c::Cosmology; kmin, kmax, nk=60, g_factor=k -> 1.0,
    rtol=1e-5)
    Ωr = Ω_r(c) * c.h^2                 # radiation density today ×h²
    ks = exp.(range(log(kmin), log(kmax); length=nk))
    Ω = [Ω_gw_rd(Pζ, k; rtol) * Ωr * g_factor(k) for k in ks]
    SIGWSpectrum(ks, _f_of_k_Hz .* ks, Ω)
end

"""
    sigw_ΔNeff(s::SIGWSpectrum, c::Cosmology)

The integrated background as extra relativistic species,
ΔN_eff = (8/7)(11/4)^{4/3} ∫ Ω_GW,0 dln k / Ω_γ,0. Check it against the
BBN/CMB bound (≲ 0.3) before trusting large peaks — and for a real gate, feed
it back into `cosmology(N_eff = 3.044 + ΔN_eff)` and rerun BBN.
"""
function sigw_ΔNeff(s::SIGWSpectrum, c::Cosmology)
    lnk = log.(s.k)
    Ω = s.Ωh² ./ c.h^2
    acc = 0.0
    for i in 1:(length(lnk)-1)
        acc += 0.5 * (Ω[i] + Ω[i+1]) * (lnk[i+1] - lnk[i])
    end
    (8 / 7) * (11 / 4)^(4 / 3) * acc / Ω_γ(c)
end