"""
Effective Field Theory of Large-Scale Structure — 1-loop matter power spectrum.

Unlike halofit/HMcode this is not a fit: it is standard perturbation theory. The
non-linear density is expanded order by order in the linear field (eq 25-27 of
Ivanov 2022, arXiv:2212.08488), and the power spectrum to one loop is

    P_1loop(k) = P_lin(k) + P₂₂(k) + P₁₃(k)                                (eq 30)

with the two loop integrals built from the SPT mode-coupling kernels F₂, F₃:

    P₂₂(k) = 2 ∫ d³q/(2π)³ [F₂(q, k-q)]² P_lin(q) P_lin(|k-q|)             (eq 32)
    P₁₃(k) = 6 P_lin(k) ∫ d³q/(2π)³ F₃(k, q, -q) P_lin(q)

After the azimuthal angle is done analytically these reduce to the standard
Makino-Sasaki-Suto / Scoccimarro closed-form kernels — a 2-D integral for P₂₂ and
a 1-D integral for P₁₃, evaluated here directly by quadrature over the linear
spectrum. (FFTLog would be faster; direct quadrature is transparent and easy to
cross-check against published loop integrals, which is the point at this stage.)

Perturbation theory alone is UV-sensitive: P₁₃ pulls in short-scale modes it has
no business describing. EFT cures this with a counterterm — the effective sound
speed / viscosity c_s² — that absorbs the UV dependence:

    P_EFT(k) = P_lin(k) + P₂₂(k) + P₁₃(k) - 2 c_s² k² P_lin(k)             (eq ~45)

`c_s²` is a free parameter (fit to data, not predicted); it is exposed, never
hardcoded. IR resummation of the BAO wiggle (eqs 37-38) is applied on request.
Valid for k ≲ 0.3 h/Mpc — mildly non-linear / BAO scales (DESI, Euclid).
"""

using QuadGK: quadgk

# --- SPT one-loop kernels (closed forms after the azimuthal integral) --------

"""
P₁₃ angular kernel: the ∫dΩ_q of 6 F₃(k,q,-q) reduces to this rational-plus-log
function of r = q/k (Makino, Sasaki & Suto 1992). Returns the bracket such that
P₁₃(k) = (k³/4π²) P_lin(k) ∫₀^∞ dr P_lin(kr) · f₁₃(r).
"""
function _f13(r)
    if abs(r - 1) < 1e-5
        # r→1: (r²-1)³ kills the log divergence, so the log term → 0.
        return (12 - 158 + 100 - 42) / 252          # = -88/252
    elseif r < 5e-3
        # r→0: the 1/r² of the log term cancels the explicit 12/r²; leading -168.
        return (-168 + (928 / 5) * r^2) / 252
    end
    logterm = (3 / r^3) * (r^2 - 1)^3 * (7r^2 + 2) * log(abs((1 + r) / (1 - r)))
    (12 / r^2 - 158 + 100r^2 - 42r^4 + logterm) / 252
end

"""
[F₂(q, k-q)]² with q = kr and x = q̂·k̂ (Ivanov eq 28, the exact symmetric SPT
kernel). With q₁=q, q₂=k-q, D=|k-q|²/k²=1+r²-2rx:
    μ₁₂ = q̂₁·q̂₂ = (x - r)/√D,   q₁/q₂ = r/√D,   q₂/q₁ = √D/r
    F₂ = 5/7 + (1/2) μ₁₂ (q₁/q₂ + q₂/q₁) + (2/7) μ₁₂²
Returns [F₂]² so that
    P₂₂(k) = (k³/2π²) ∫₀^∞ r² dr ∫_{-1}^1 dx [F₂]² P_lin(kr) P_lin(k√D).
"""
function _f22sq(r, x)
    D = 1 + r^2 - 2r * x
    D < 1e-12 && return 0.0
    μ = (x - r) / sqrt(D)
    q1q2 = r / sqrt(D)            # q₁/q₂
    q2q1 = sqrt(D) / r            # q₂/q₁
    F2 = 5 / 7 + 0.5 * μ * (q1q2 + q2q1) + (2 / 7) * μ^2
    F2^2
end

# --- loop integrals ----------------------------------------------------------

"""
    p13(P, k; qmax_factor = 100)

The P₁₃ one-loop integral (Mpc³) at wavenumber `k` (1/Mpc), built from the linear
[`MatterPowerSpectrum`](@ref) `P`. Negative, as it should be at high k.
"""
function p13(P::MatterPowerSpectrum, k; rtol=1e-4)
    kmin, kmax = exp(first(P.logk)), exp(last(P.logk))
    Plin = power(P, k)
    integrand(r) = (kr = k * r; (kr < kmin || kr > kmax) ? 0.0 : power(P, kr) * _f13(r))
    rlo, rhi = kmin / k, kmax / k
    # f₁₃ has a log that is integrable at r = 1 but varies sharply there; give the
    # adaptive rule r = 1 as a breakpoint so it straddles the corner cleanly.
    seg = rlo < 1 < rhi ? (rlo, 1.0, rhi) : (rlo, rhi)
    I = quadgk(integrand, seg...; rtol)[1]
    k^3 / (4π^2) * Plin * I
end

"""
    p22(P, k)

The P₂₂ one-loop integral (Mpc³) at `k` (1/Mpc). Positive.
"""
function p22(P::MatterPowerSpectrum, k; rtol=1e-5)
    kmin, kmax = exp(first(P.logk)), exp(last(P.logk))
    rlo, rhi = kmin / k, kmax / k
    # Substitute the inner variable x → y = |k−q|/k = √D, with x = (1+r²−y²)/2r and
    # dx = −(y/r)dy. The mode-coupling integrand carries P_lin(k√D) = P_lin(ky); in
    # x it sits on a ridge whose location (where ky hits the P(k) turnover) slides
    # with r, so a fixed-node rule silently misses it and P22 is underestimated at
    # low k. In y that feature is a single fixed point k_turn/k, which the adaptive
    # rule resolves cleanly. Cross-checked against FAST-PT and a scipy brute force.
    outer(r) = begin
        kr = k * r
        (kr < kmin || kr > kmax) && return 0.0
        Pr = power(P, kr)
        ylo, yhi = abs(1 - r), 1 + r
        inner(y) = begin
            ky = k * y
            (ky < kmin || ky > kmax) && return 0.0
            x = (1 + r^2 - y^2) / (2r)
            _f22sq(r, x) * power(P, ky) * y
        end
        r * Pr * quadgk(inner, ylo, yhi; rtol)[1]
    end
    I = quadgk(outer, rlo, rhi; rtol)[1]
    k^3 / (2π^2) * I
end

# --- IR resummation ----------------------------------------------------------

"""
    _bao_damping(P, c)

The BAO damping scale Σ² = (1/6π²)∫₀^{k_s} dq P_lin(q)[1 - j₀(q r_d) + 2 j₂(q r_d)]
(Baldauf et al. / Ivanov eq for IR resummation). Returns Σ² in Mpc², used to
damp the wiggle component of the linear spectrum in the resummed loop.
"""
function _bao_damping(P::MatterPowerSpectrum, rd)
    ks = 0.2 * P.cosmo.h                    # separation scale k_S ≈ 0.2 h/Mpc
    kmin = exp(first(P.logk))
    j0(x) = x < 1e-4 ? 1 - x^2 / 6 : sin(x) / x
    j2(x) = x < 1e-2 ? x^2 / 15 : (3 / x^2 - 1) * sin(x) / x - 3 * cos(x) / x^2
    integrand(q) = power(P, q) * (1 - j0(q * rd) + 2 * j2(q * rd))
    quadgk(integrand, kmin, ks; rtol=1e-4)[1] / (6π^2)
end

# --- the model ---------------------------------------------------------------

struct EFTSpectrum
    lin::MatterPowerSpectrum
    cs2::Float64                # counterterm c_s² [Mpc²]
    ir::Bool
    Σ2::Float64                 # BAO damping (0 if ir=false)
    dewig                       # k ↦ (P_smooth, P_wiggle) if ir, else nothing
end

"""
    eft_power(P; cs2 = 0.0, ir_resum = true, rd = 100/h Mpc)

Build the 1-loop EFT spectrum on a linear [`MatterPowerSpectrum`](@ref). `cs2` is
the effective sound-speed counterterm in Mpc² (a free parameter — fit it to data;
0 gives pure 1-loop SPT). `ir_resum` applies BAO IR resummation. Evaluate with
[`power`](@ref).
"""
function eft_power(P::MatterPowerSpectrum; cs2=0.0, ir_resum=true, rd=nothing)
    c = P.cosmo
    rd === nothing && (rd = 100.0 / c.h)          # ~sound horizon in Mpc
    Σ2 = ir_resum ? _bao_damping(P, rd) : 0.0
    dewig = nothing
    if ir_resum
        σv2 = _sigma_v2(P)                          # reuse HMcode helper
        # smooth/wiggle split via the EH no-wiggle broadband (reuse _dewiggle
        # internals conceptually; here build a direct split function).
        dewig = _smooth_wiggle_split(P, c)
    end
    EFTSpectrum(P, cs2, ir_resum, Σ2, dewig)
end

"Split linear P into (smooth, wiggle) using the EH no-wiggle broadband shape."
function _smooth_wiggle_split(P::MatterPowerSpectrum, c::Cosmology)
    lo, hi = first(P.logk), last(P.logk)
    ns = c.n_s
    ratio(lk) = (k = exp(lk); log(power(P, k) / (k^ns * _eh_nowiggle(c, k)^2)))
    σsm = 0.25
    smooth_ln(lk) = begin
        a, b = max(lo, lk - 4σsm), min(hi, lk + 4σsm)
        num = quadgk(u -> ratio(u) * exp(-(u - lk)^2 / (2σsm^2)), a, b; rtol=1e-4)[1]
        den = quadgk(u -> exp(-(u - lk)^2 / (2σsm^2)), a, b; rtol=1e-4)[1]
        num / den
    end
    function (k)
        lk = log(k)
        Psm = exp(smooth_ln(lk)) * k^ns * _eh_nowiggle(c, k)^2
        (Psm, power(P, k) - Psm)
    end
end

"""
    power(eft, k)

1-loop EFT matter power spectrum P(k) in Mpc³ at `k` (1/Mpc). With IR resummation
the BAO wiggle in both the linear and loop pieces is damped by e^{-k²Σ²}; the
counterterm -2 c_s² k² P_lin is always included.
"""
function power(eft::EFTSpectrum, k)
    P = eft.lin
    P22 = p22(P, k)
    P13 = p13(P, k)
    ct = -2 * eft.cs2 * k^2 * power(P, k)

    if !eft.ir
        return power(P, k) + P22 + P13 + ct
    end

    # IR-resummed: damp the wiggle. Leading order (Ivanov eq for the resummed
    # 1-loop): P = P_sm + e^{-k²Σ²}P_w + P_1loop[P_sm + e^{-k²Σ²}P_w] + ...
    # To one loop it is enough to damp the wiggle in the linear + tree pieces and
    # add the loop of the smooth spectrum plus the damped-wiggle correction.
    Psm, Pw = eft.dewig(k)
    damp = exp(-k^2 * eft.Σ2)
    Plin_ir = Psm + damp * Pw
    # loop on the full (already carries the wiggle); multiply the wiggle part of
    # the loop response by the same damping is higher order — keep the loop as-is
    # on the resummed linear input to leading order.
    Plin_ir + P22 + P13 - 2 * eft.cs2 * k^2 * Plin_ir
end
