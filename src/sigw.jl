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
using Random: MersenneTwister

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

# --- primordial non-Gaussianity in the induced background --------------------
#
# When the curvature field is local non-Gaussian, ζ = ζ_g + F_NL(ζ_g²−⟨ζ_g²⟩) +
# G_NL ζ_g³ (with F_NL = 3/5 f_NL, G_NL = 9/25 g_NL), the induced tensor spectrum
# picks up a four-point structure beyond the Gaussian two-point convolution. We
# follow the diagrammatic decomposition of arXiv:2308.07155, whose kernel is the
# same RD kernel as the Gaussian piece, now evaluated off-diagonal.
#
# What is implemented here, exactly: the joint kernel I(u,v)I(u',v') (their
# eq. 17), the linear-g_NL correction (eq. 50, a closed rescaling of the Gaussian
# spectrum by the curvature variance), and the *complete* f_NL² correction — the
# disconnected "hybrid" (eq. 28) plus the two connected diagrams "Z" (eq. 30) and
# "C" (eq. 37), the latter two by importance-sampled Monte Carlo. The total f_NL²
# induced background is `sigw_hybrid + sigw_Z_connected + sigw_C_connected`. The
# entire *disconnected* NG sector is also closed: the reducible f_NL⁴ term
# (`sigw_f4_reducible`, eq. 42) and the loop g_NL tower (`sigw_gnl_reducible`,
# eqs. 50/54/66). Still remaining: the *connected* higher-order diagrams —
# f_NL⁴ planar/non-planar (eqs. 44, 48) and g_NL² tri/ring, g_NL³ ring3
# (eqs. 56–68) — each a 6–8D integral; their expressions live in the same
# reference and reuse `_sigw_II` below.

# B-factor of the RD kernel: the non-oscillatory log part (2308.07155 eq. 17).
@inline _sigw_Bfac(u, v) =
    (u^2 + v^2 - 3) * log(abs(((u - v)^2 - 3) / ((u + v)^2 - 3))) + 4u * v

"""
    _sigw_II(u, v, u′, v′)

The time-averaged product of RD kernels I(u,v)·I(u′,v′) (arXiv:2308.07155
eq. 17). The π² term is the resonant contribution switched on above u+v = √3.
Its diagonal `_sigw_II(u,v,u,v)` is the I²(u,v) of the Gaussian spectrum.
"""
@inline function _sigw_II(u, v, u′, v′)
    s = u^2 + v^2 - 3
    s′ = u′^2 + v′^2 - 3
    pref = 9 * s * s′ / (1024 * u^3 * v^3 * u′^3 * v′^3)
    proj = (4u^2 - (u^2 - v^2 + 1)^2) * (4u′^2 - (u′^2 - v′^2 + 1)^2)
    r3 = sqrt(3.0)
    res = (u + v > r3 && u′ + v′ > r3) ? π^2 * s * s′ : 0.0
    pref * proj * (_sigw_Bfac(u, v) * _sigw_Bfac(u′, v′) + res)
end
@inline _sigw_I2(u, v) = _sigw_II(u, v, u, v)

# v-limits inside [|1−u|, 1+u] intersected with the support window [klo, khi]
# (the arguments uk, vk must land where 𝒫_ζ ≠ 0), with the u+v=√3 resonance
# handed to the quadrature as a breakpoint. Returns nothing if the window is empty.
function _sigw_vsegs_supp(u, klo, khi)
    lo = max(abs(1 - u), klo)
    hi = min(1 + u, khi)
    hi <= lo && return nothing
    vr = sqrt(3.0) - u
    lo < vr < hi ? (lo, vr, hi) : (lo, hi)
end

# per-e-fold breakpoints across [a, b] so a sharply peaked 𝒫_ζ is never straddled
# by a single wide quadrature interval.
_sigw_ubks(a, b) = exp.(range(log(a), log(b); length=max(4, ceil(Int, log(b / a)) + 1)))

# effective support: the sub-window of [kmin, kmax] where 𝒫_ζ exceeds `rel`·peak.
# Integrating the kernel over the full nominal window wastes refinement in the
# deep tails, where 𝒫_ζ underflows toward zero while the kernel's 1/v⁶ grows —
# a spurious dynamic range that has no bearing on the integral.
function _sigw_support(Pζ, kmin, kmax; rel=1e-8, n=400)
    ks = exp.(range(log(kmin), log(kmax); length=n))
    ps = Pζ.(ks)
    pmax = maximum(ps)
    pmax <= 0 && return (kmin, kmax)
    idx = findall(≥(rel * pmax), ps)
    isempty(idx) ? (kmin, kmax) : (ks[first(idx)], ks[last(idx)])
end

"""
    sigw_omega_g_kernel(Pζ, k; ksupp, rtol = 1e-4)

The Gaussian induced-GW density during RD written in the (u,v) kernel variables
(arXiv:2308.07155 eq. 24), Ω = (1/3)∫du∫dv I²(u,v)/(u²v²) 𝒫_ζ(uk)𝒫_ζ(vk). The
integration is confined to the curvature-spectrum support `ksupp = (kmin, kmax)`,
where 𝒫_ζ(uk)𝒫_ζ(vk) ≠ 0. Numerically identical to [`Ω_gw_rd`](@ref); kept as the
cross-check that fixes the kernel normalisation the NG terms inherit.
"""
function sigw_omega_g_kernel(Pζ, k; ksupp, rtol=1e-4)
    es = _sigw_support(Pζ, ksupp[1], ksupp[2])
    klo, khi = es[1] / k, es[2] / k
    inner(u) = (segs = _sigw_vsegs_supp(u, klo, khi);
    segs === nothing ? 0.0 :
    quadgk(v -> _sigw_I2(u, v) / (u^2 * v^2) * Pζ(u * k) * Pζ(v * k), segs...; rtol)[1])
    quadgk(inner, _sigw_ubks(klo, khi)...; rtol)[1] / 3
end

# phase-space convolution 𝒞(x) = ∫du₁∫dv₁ 𝒫(u₁x)𝒫(v₁x)/(u₁²v₁²): the kernel-free
# inner loop of the disconnected f_NL² term (2308.07155 eq. 28), confined to the
# support u₁x, v₁x ∈ [kmin, kmax].
function _sigw_conv(Pζ, x; ksupp, rtol=1e-4)
    klo, khi = ksupp[1] / x, ksupp[2] / x
    inner(u1) = begin
        lo = max(abs(1 - u1), klo)
        hi = min(1 + u1, khi)
        hi <= lo ? 0.0 :
        quadgk(v1 -> Pζ(u1 * x) * Pζ(v1 * x) / (u1^2 * v1^2), lo, hi; rtol)[1]
    end
    quadgk(inner, _sigw_ubks(klo, khi)...; rtol)[1]
end

# 𝒞(x) is the same across every outer node of the hybrid integral, so tabulate it
# once on a log-x grid and interpolate (log-x, linear value). It is peaked near
# the spectrum support and falls to zero on either side, so values outside
# [xlo, xhi] are taken as 0.
function _sigw_conv_table(Pζ, xlo, xhi; n=48, ksupp, rtol=1e-3)
    lx = collect(range(log(xlo), log(xhi); length=n))
    vals = [_sigw_conv(Pζ, exp(l); ksupp, rtol) for l in lx]
    function conv(x)
        (x ≤ xlo || x ≥ xhi) && return 0.0
        l = log(x)
        j = clamp(searchsortedlast(lx, l), 1, n - 1)
        w = (l - lx[j]) / (lx[j+1] - lx[j])
        (1 - w) * vals[j] + w * vals[j+1]
    end
    conv
end

# the RD kernel integral in the fast (t,s) variables of `Ω_gw_rd`, but for an
# arbitrary product g(ku, kv) in place of 𝒫_ζ(ku)𝒫_ζ(kv). The kernel is even in s,
# so any asymmetry between the two arguments averages out over the domain. With
# g = 𝒫_ζ·𝒫_ζ this reproduces `Ω_gw_rd` exactly; the hybrid NG term uses g = 𝒞·𝒫_ζ.
function _sigw_omega_ts(g, k; tmax=500.0, rtol=1e-4)
    inner(t) = quadgk(s -> begin
            u = (t + s + 1) / 2
            v = (t - s + 1) / 2
            _kt_proj(t, s) * _kt_x2I2(t, s) * g(k * u, k * v)
        end, -1.0, 1.0; rtol)[1]
    quadgk(inner, 0.0, sqrt(3) - 1, 2.0, 20.0, tmax; rtol)[1] / 12
end

"""
    sigw_hybrid(Pζ, k; f_NL, ksupp, conv = nothing, rtol = 1e-4)

The disconnected ("hybrid") f_NL² contribution to Ω_GW during RD
(arXiv:2308.07155 eq. 28):

    Ω_hybrid = (2/3)F_NL² ∫du∫dv I²(u,v)/(u²v²) 𝒫_ζ(vk) 𝒞(uk) = 2F_NL² · Ω_ts[𝒞·𝒫_ζ],

with F_NL = 3/5 f_NL and 𝒞 the phase-space convolution above. Because the term
shares the Gaussian kernel, it is evaluated through the same fast (t,s)
integrator (`_sigw_omega_ts`) rather than the stiff (u,v) form. `ksupp =
(kmin, kmax)` is the curvature-spectrum support; pass a prebuilt `conv` closure
(from [`_sigw_conv_table`](@ref)) to reuse the tabulation across wavenumbers.
This is the disconnected part of the f_NL² correction; the connected pieces are
[`sigw_Z_connected`](@ref) and [`sigw_C_connected`](@ref) — the three sum to the
full f_NL² term.
"""
function sigw_hybrid(Pζ, k; f_NL, ksupp, conv=nothing, rtol=1e-4)
    F = 3 / 5 * f_NL
    es = _sigw_support(Pζ, ksupp[1], ksupp[2])
    cf = conv === nothing ? _sigw_conv_table(Pζ, es[1], es[2]; ksupp=es, rtol) : conv
    2 * F^2 * _sigw_omega_ts((a, b) -> cf(a) * Pζ(b), k; rtol)
end

# inverse-CDF sampler of ln k weighted by 𝒫_ζ over [klo, khi]: returns the log-k
# grid, its normalised cumulative ∫𝒫 dln k, and the total A = ∫𝒫 dln k. Used to
# importance-sample the spectrum legs of the connected term so the sharply peaked
# 𝒫_ζ factors cancel against the sampling density instead of wrecking the variance.
function _sigw_lnk_cdf(Pζ, klo, khi; n=800)
    lnk = collect(range(log(klo), log(khi); length=n))
    p = [Pζ(exp(l)) for l in lnk]
    cdf = zeros(n)
    for i in 2:n
        cdf[i] = cdf[i-1] + 0.5 * (p[i] + p[i-1]) * (lnk[i] - lnk[i-1])
    end
    A = cdf[n]
    A > 0 && (cdf ./= A)
    (lnk, cdf, A)
end
@inline function _sigw_draw(lnk, cdf, r)
    j = clamp(searchsortedfirst(cdf, r), 2, length(cdf))
    c0, c1 = cdf[j-1], cdf[j]
    w = c1 > c0 ? (r - c0) / (c1 - c0) : 0.0
    exp(lnk[j-1] + w * (lnk[j] - lnk[j-1]))
end
@inline _sigw_cth(u, v) = clamp((1 + u^2 - v^2) / (2u), -1.0, 1.0)

"""
    sigw_Z_connected(Pζ, k; f_NL, ksupp, N = 8_000_000, seed = 12345) -> (Ω_Z, err)

The connected "Z" contribution to the f_NL² induced-GW density during RD
(arXiv:2308.07155 eqs. 29–30):

    Ω_Z = (2/3π)F_NL² ∫du dv du' dv' dφ₁ cos2φ₁ · I(u,v)I(u',v') · uu'/(v²v'²w³) · 𝒫_ζ(vk)𝒫_ζ(v'k)𝒫_ζ(wk),

with w = w₀₁₂ = |k−q−q'|/k (their eq. 32) and F_NL = 3/5 f_NL. This is a genuine
five-dimensional integral with an oscillatory cos2φ₁ and no factorisation, so it
is evaluated by importance-sampled Monte Carlo: the two spectrum legs v, v' are
drawn ∝ 𝒫_ζ (cancelling the peaked factors), giving

    Ω_Z = (4/3)F_NL² A² · E[cos2φ₁ · I(u,v)I(u',v') · uu'/(vv'w³) · 𝒫_ζ(wk) · W_u W_u'],

A = ∫𝒫_ζ dln k. Returns the estimate and its 1σ Monte-Carlo error; use
`Threads.nthreads() > 1` for speed. This is *one* of the two connected f_NL²
diagrams — the companion "C" term (their eq. 37) is not yet included; its free
momentum has no direct spectrum leg and needs an angle-solving sampler.
"""
function sigw_Z_connected(Pζ, k; f_NL, ksupp, N=8_000_000, seed=12345)
    F = 3 / 5 * f_NL
    es = _sigw_support(Pζ, ksupp[1], ksupp[2])
    lnk, cdf, A = _sigw_lnk_cdf(Pζ, es[1], es[2])
    A > 0 || return (0.0, 0.0)
    nt = max(Threads.nthreads(), 1)
    S = zeros(nt)
    S2 = zeros(nt)
    cn = zeros(Int, nt)
    Threads.@threads for t in 1:nt
        rng = MersenneTwister(seed + t)
        s = 0.0
        s2 = 0.0
        c = 0
        for _ in 1:(N ÷ nt)
            v = _sigw_draw(lnk, cdf, rand(rng)) / k
            vp = _sigw_draw(lnk, cdf, rand(rng)) / k
            Wu = (1 + v) - abs(1 - v)
            u = abs(1 - v) + Wu * rand(rng)
            Wup = (1 + vp) - abs(1 - vp)
            up = abs(1 - vp) + Wup * rand(rng)
            φ = 2π * rand(rng)
            ct = _sigw_cth(u, v)
            st = sqrt(1 - ct^2)
            ctp = _sigw_cth(up, vp)
            stp = sqrt(1 - ctp^2)
            w2 = 1 + u^2 + up^2 + 2u * up * (st * stp * cos(φ) + ct * ctp) - 2u * ct - 2up * ctp
            w2 <= 0 && continue
            w = sqrt(w2)
            val = cos(2φ) * _sigw_II(u, v, up, vp) * (u * up / (v * vp * w^3)) *
                  Pζ(w * k) * Wu * Wup
            isfinite(val) || continue
            s += val
            s2 += val^2
            c += 1
        end
        S[t] = s
        S2[t] = s2
        cn[t] = c
    end
    Nt = sum(cn)
    Nt == 0 && return (0.0, 0.0)
    m = sum(S) / Nt
    var = max(sum(S2) / Nt - m^2, 0.0)
    pref = (4 / 3) * F^2 * A^2
    (pref * m, pref * sqrt(var / Nt))
end

"""
    sigw_C_connected(Pζ, k; f_NL, ksupp, N = 32_000_000, seed = 4242) -> (Ω_C, err)

The connected "C" contribution to the f_NL² induced-GW density during RD
(arXiv:2308.07155 eqs. 36–37):

    Ω_C = (2/3π)F_NL² ∫du dv du' dv' dφ₁ cos2φ₁ · I(u,v)I(u',v') · uv/(u'²v'²w³) · 𝒫_ζ(u'k)𝒫_ζ(v'k)𝒫_ζ(wk),

with w = w₁₂ = |q−q'|/k (their eq. 38). Here the free momentum q = uk carries no
direct spectrum leg, so the spectrum-leg importance sampling used for
[`sigw_Z_connected`](@ref) is not enough; instead the third leg w is drawn ∝ 𝒫_ζ
and the azimuth φ₁ is *solved* from its definition (|cosφ₁| ≤ 1 gates the sample),
using |dw/dφ₁| = uu'sθsθ'|sinφ₁|/w. Sampling u', v', w ∝ 𝒫_ζ and u, v uniformly,

    Ω_C = (4/3π)F_NL² A³ · E[cos2φ₁ · I(u,v)I(u',v') · v W_u W_v / (u'²v' w sθ sθ' |sinφ₁|)],

the expectation over all trials (rejected samples count as 0). Returns the
estimate and its 1σ Monte-Carlo error. Validated against the monochromatic
collapse (independent 2D quadrature) as the spectrum narrows.
"""
function sigw_C_connected(Pζ, k; f_NL, ksupp, N=32_000_000, seed=4242)
    F = 3 / 5 * f_NL
    es = _sigw_support(Pζ, ksupp[1], ksupp[2])
    lnk, cdf, A = _sigw_lnk_cdf(Pζ, es[1], es[2])
    A > 0 || return (0.0, 0.0)
    ulo = (es[1] / k) * 0.3
    uhi = (es[2] / k) * 3.0
    Wu = uhi - ulo
    nt = max(Threads.nthreads(), 1)
    S = zeros(nt)
    S2 = zeros(nt)
    tot = zeros(Int, nt)
    Threads.@threads for t in 1:nt
        rng = MersenneTwister(seed + t)
        s = 0.0
        s2 = 0.0
        tt = 0
        for _ in 1:(N ÷ nt)
            tt += 1
            up = _sigw_draw(lnk, cdf, rand(rng)) / k
            vp = _sigw_draw(lnk, cdf, rand(rng)) / k
            w = _sigw_draw(lnk, cdf, rand(rng)) / k
            u = ulo + Wu * rand(rng)
            Wv = (1 + u) - abs(1 - u)
            v = abs(1 - u) + Wv * rand(rng)
            ct = _sigw_cth(u, v)
            st = sqrt(1 - ct^2)
            ctp = _sigw_cth(up, vp)
            stp = sqrt(1 - ctp^2)
            (st == 0 || stp == 0) && continue
            cφ = (u^2 + up^2 - w^2 - 2u * up * ct * ctp) / (2u * up * st * stp)
            abs(cφ) > 1 && continue
            sφ = sqrt(1 - cφ^2)
            sφ == 0 && continue
            val = cos(2 * acos(cφ)) * _sigw_II(u, v, up, vp) * v * Wu * Wv /
                  (up^2 * vp * w * st * stp * sφ)
            isfinite(val) || continue
            s += val
            s2 += val^2
        end
        S[t] = s
        S2[t] = s2
        tot[t] = tt
    end
    Nt = sum(tot)
    Nt == 0 && return (0.0, 0.0)
    m = sum(S) / Nt
    var = max(sum(S2) / Nt - m^2, 0.0)
    pref = (4 / (3π)) * F^2 * A^3
    (pref * m, pref * sqrt(var / Nt))
end

"""
    sigw_f4_reducible(Pζ, k; f_NL, ksupp, conv = nothing, rtol = 1e-4)

The disconnected ("reducible") f_NL⁴ contribution to Ω_GW during RD
(arXiv:2308.07155 eq. 42): with *both* external spectrum legs replaced by the
phase-space convolution 𝒞,

    Ω_reducible = (F_NL⁴/3) ∫du dv I²(u,v)/(u²v²) 𝒞(uk) 𝒞(vk) = F_NL⁴ · Ω_ts[𝒞·𝒞],

F_NL = 3/5 f_NL. Like the hybrid it shares the Gaussian kernel and is evaluated
through the fast (t,s) integrator, reusing the convolution table. This is the
disconnected part of the f_NL⁴ correction; the connected "planar" (eq. 44) and
"non-planar" (eq. 48) diagrams are not yet included (see the module note).
"""
function sigw_f4_reducible(Pζ, k; f_NL, ksupp, conv=nothing, rtol=1e-4)
    F = 3 / 5 * f_NL
    es = _sigw_support(Pζ, ksupp[1], ksupp[2])
    cf = conv === nothing ? _sigw_conv_table(Pζ, es[1], es[2]; ksupp=es, rtol) : conv
    F^4 * _sigw_omega_ts((a, b) -> cf(a) * cf(b), k; rtol)
end

"""
    sigw_gnl_reducible(Pζ; kmin, kmax, g_NL, rtol = 1e-6)

The full reducible (loop) g_NL rescaling of the Gaussian spectrum, summing the
closed-form disconnected g_NL diagrams (arXiv:2308.07155 eqs. 50, 54, 66):

    Ω_g_NL,red(k) = (12 G A + 54 G²A² + 108 G³A³) · Ω_g(k),

with G = 9/25 g_NL and A = ∫𝒫_ζ dln k. Returns that multiplicative factor. The
linear term alone is [`sigw_gnl_factor`](@ref) (what `sigw_spectrum_ng` uses); the
connected g_NL² (tri, ring) and g_NL³ (ring3) diagrams are not included here.
"""
function sigw_gnl_reducible(Pζ; kmin, kmax, g_NL, rtol=1e-6)
    A = quadgk(lnp -> Pζ(exp(lnp)), log(kmin), log(kmax); rtol)[1]
    G = (9 / 25) * g_NL
    GA = G * A
    12 * GA + 54 * GA^2 + 108 * GA^3
end

"""
    sigw_gnl_factor(Pζ; kmin, kmax, g_NL, rtol = 1e-6)

The linear-g_NL correction is a pure rescaling of the Gaussian spectrum
(arXiv:2308.07155 eq. 50): Ω_{g_NL}(k) = 12 G_NL A · Ω_g(k), with G_NL = 9/25 g_NL
and A = ∫dln p 𝒫_ζ(p) the variance of the Gaussian curvature. This returns the
multiplicative factor 12 G_NL A to apply to the Gaussian Ω_GW.
"""
function sigw_gnl_factor(Pζ; kmin, kmax, g_NL, rtol=1e-6)
    A = quadgk(lnp -> Pζ(exp(lnp)), log(kmin), log(kmax); rtol)[1]
    12 * (9 / 25) * g_NL * A
end

"""
    sigw_spectrum_ng(Pζ, c::Cosmology; kmin, kmax, nk = 40, f_NL = 0, g_NL = 0,
                     connected = false, f4 = false, gnl_loops = false,
                     N = 16_000_000, rtol = 1e-3, g_factor = k -> 1.0)

Today's induced-GW spectrum Ω_GW,0(k)h² including the primordial-NG corrections.
The default sums the Gaussian piece, the disconnected f_NL² ("hybrid") term, and
the linear-g_NL rescaling. The flags add the further validated contributions:

  * `connected` — the connected f_NL² diagrams [`sigw_Z_connected`](@ref) +
    [`sigw_C_connected`](@ref) (Monte Carlo with `N` samples each per k),
    completing the full f_NL² correction;
  * `f4` — the disconnected f_NL⁴ [`sigw_f4_reducible`](@ref);
  * `gnl_loops` — the full loop g_NL tower [`sigw_gnl_reducible`](@ref)
    (12GA+54G²A²+108G³A³) in place of the linear rescaling.

`kmin`/`kmax` are the curvature-spectrum support. Still not included: the
connected higher-order diagrams (f_NL⁴ planar/non-planar, g_NL² tri/ring,
g_NL³ ring3 — see the module note). With `f_NL = g_NL = 0` this returns exactly
the Gaussian [`sigw_spectrum`](@ref).
"""
function sigw_spectrum_ng(Pζ, c::Cosmology; kmin, kmax, nk=40, f_NL=0.0, g_NL=0.0,
    connected=false, f4=false, gnl_loops=false, N=16_000_000, rtol=1e-3,
    g_factor=k -> 1.0)
    Ωr = Ω_r(c) * c.h^2
    ksupp = _sigw_support(Pζ, kmin, kmax)   # trim to where 𝒫_ζ actually lives
    gfac = g_NL == 0 ? 0.0 :
           gnl_loops ? sigw_gnl_reducible(Pζ; kmin, kmax, g_NL) :
           sigw_gnl_factor(Pζ; kmin, kmax, g_NL)
    conv = (f_NL == 0 && !f4) ? nothing :
           _sigw_conv_table(Pζ, ksupp[1], ksupp[2]; ksupp, rtol)
    ks = exp.(range(log(kmin), log(kmax); length=nk))
    Ω = map(ks) do k
        # Gaussian piece via the fast validated (t,s) integrator; the (u,v) kernel
        # form gives the identical number but is kept only as the cross-check.
        Ωg = Ω_gw_rd(Pζ, k; rtol)
        Ωng = Ωg * (1 + gfac)
        if f_NL != 0
            Ωng += sigw_hybrid(Pζ, k; f_NL, ksupp, conv, rtol)
            if connected
                Ωng += sigw_Z_connected(Pζ, k; f_NL, ksupp, N)[1]
                Ωng += sigw_C_connected(Pζ, k; f_NL, ksupp, N)[1]
            end
            f4 && (Ωng += sigw_f4_reducible(Pζ, k; f_NL, ksupp, conv, rtol))
        end
        Ωng * Ωr * g_factor(k)
    end
    SIGWSpectrum(ks, _f_of_k_Hz .* ks, Ω)
end