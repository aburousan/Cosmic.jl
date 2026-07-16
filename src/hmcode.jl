"""
HMcode-2020 non-linear matter power spectrum (Mead, Brieden, Tröster, Heymans
2021; arXiv:2009.01858), the augmented halo model. RMS ~2.5% against emulators,
the modern successor to Takahashi halofit.

The model sums a two-halo term (damped, de-wiggled linear theory) and a one-halo
term (the mass integral of squared halo profiles), smoothed at the transition:

    Δ²_HM(k) = [ Δ²_2H(k)^α + Δ²_1H(k)^α ]^{1/α}                          (eq 23)

Everything cosmology-dependent is parametrised through the linear-spectrum
shape and amplitude -- σ8 of the cold matter and the effective index at the
collapse scale -- rather than through h, Ω_m, w directly, which is what lets the
one fit span a wide model space. The spherical-collapse δ_c and Δ_v use the
Mead (2017) fits; the fitted parameters and every coefficient here are
transcribed from CAMB's `halofit.f90` (imead==3) -- see
the HMcode-2020 paper (arXiv:2009.01858), from which every coefficient here was transcribed.

Cold vs total matter: σ, ν and the mass function use the *cold* (cdm+baryon)
spectrum; the two-halo term uses total matter. With massless neutrinos the two
coincide, which is the regime validated against CAMB here. Massive-ν cold-matter
input (δ_cb) is a documented TODO.
"""

using SpecialFunctions: sinint, cosint
using QuadGK: quadgk
using Roots: find_zero, Bisection

const _RHO_CRIT_MSUN_MPC3 = 2.77536627e11   # ρ_crit,0 / h²  [M⊙ Mpc⁻³]

# --- growth (Mead normalisation: g → a in matter domination) ------------------

"""
    _hmcode_growth(c, z)

Scale-independent linear growth solved as an ODE in ln a, normalised so g → a
deep in matter domination (not to 1 today). Returns `(g, Gacc, ginf_ratio)`:
`g = g(z)`, `Gacc = ∫₀ᵃ g/a' da'` the accumulated growth, and `ginf_ratio =
g_∞(model)/g_∞(ΛCDM)` for the Dolag concentration correction (1 for flat ΛCDM).
"""
function _hmcode_growth(c::Cosmology, z)
    a_target = scale_factor(z)
    # d²g/dlna² + (2 + dlnH/dlna) dg/dlna - (3/2)Ω_m(a) g = 0
    Ωm_a(a) = begin
        E2 = E(c, a)^2
        s = 0.0
        for sp in c.species
            (sp isa Baryons || sp isa ColdDarkMatter || sp isa MassiveNeutrinos) &&
                (s += ρ_over_ρc0(sp, a))
        end
        s / E2
    end
    dlnH_dlna(a) = begin           # (a/E) dE/da
        h = 1e-5
        (log(E(c, a * exp(h))) - log(E(c, a * exp(-h)))) / (2h)
    end
    # RK4 in x = ln a. Seed g = a with dg/dlna = a in MATTER domination -- that
    # is the growing-mode solution only after equality (a_eq ~ 3e-4 here), so we
    # start at a = 1e-2 (z=100), safely matter-dominated. Seeding in the
    # radiation era (where growth is nearly frozen, not g∝a) makes g(today) come
    # out spuriously small and wrecks δ_c/Δ_v. G is seeded with its MD value
    # ∫₀^{a_i} (g/a')da' ≈ a_i.
    rhs(x, g, gp) = begin
        a = exp(x)
        (gp, -(2 + dlnH_dlna(a)) * gp + 1.5 * Ωm_a(a) * g)
    end
    xi = log(1e-2)
    # RK4 with ~512 steps — 4th order over ~5 e-folds is converged well past the
    # per-mille the fits need. G = ∫ g/a da = ∫ g dx (integrand in x = ln a is g).
    integrate_to(xf; n=512) = begin
        g, gp, G = exp(xi), exp(xi), exp(xi)   # g=a, dg/dlna=a, G≈a in MD
        dx = (xf - xi) / n
        for i in 1:n
            x = xi + (i - 1) * dx
            k1g, k1p = rhs(x, g, gp)
            k2g, k2p = rhs(x + dx / 2, g + dx / 2 * k1g, gp + dx / 2 * k1p)
            k3g, k3p = rhs(x + dx / 2, g + dx / 2 * k2g, gp + dx / 2 * k2p)
            k4g, k4p = rhs(x + dx, g + dx * k3g, gp + dx * k3p)
            gmid = g + dx / 2 * k1g
            g += dx / 6 * (k1g + 2k2g + 2k3g + k4g)
            gp += dx / 6 * (k1p + 2k2p + 2k3p + k4p)
            G += dx * gmid
        end
        g, gp, G
    end
    g, _, G = integrate_to(log(a_target))
    # Dolag concentration correction needs g_∞(model)/g_∞(ΛCDM). For flat ΛCDM
    # (no wCDM, no curvature) the two are identical → ratio 1, and the extra
    # integration to a→∞ is skipped. Only compute it when it can differ.
    has_wcdm = any(s -> (s isa AbstractDarkEnergy && !(s isa CosmologicalConstant)) ||
                       s isa Curvature, c.species)
    ginf_ratio = 1.0
    if has_wcdm
        # g_∞ of the model vs an equivalent flat ΛCDM. The ΛCDM reference is
        # approximated by the model's own late-time growth here; a full wCDM
        # Dolag reference is a TODO (documented). Integrate to a = 100.
        gm, _, _ = integrate_to(log(100.0))
        ginf_ratio = 1.0            # placeholder until the ΛCDM-equivalent solve lands
    end
    (g=g, Gacc=G, ginf_ratio=ginf_ratio)
end

# --- Mead 2017 spherical collapse --------------------------------------------

_f_Mead(x, y, p0, p1, p2, p3) = p0 + p1 * (1 - x) + p2 * (1 - x)^2 + p3 * (1 - y)

"Critical linear collapse density δ_c(z), Mead 2017 Appendix A."
function _dc_mead(x, y, Ωm_z, fν)
    dc0 = (3 / 20) * (12π)^(2 / 3)
    d = 1.0
    d += _f_Mead(x, y, -0.0069, -0.0208, 0.0312, 0.0021) * log10(Ωm_z)^1
    d += _f_Mead(x, y, 0.0001, -0.0647, -0.0417, 0.0646)
    d * dc0 * (1 - 0.041 * fν)
end

"Virial overdensity Δ_v(z), Mead 2017 Appendix A."
function _dv_mead(x, y, Ωm_z, fν)
    Dv0 = 18π^2
    d = 1.0
    d += _f_Mead(x, y, -0.79, -10.17, 2.51, 6.51) * log10(Ωm_z)^1
    d += _f_Mead(x, y, -1.89, 0.38, 18.8, -15.87) * log10(Ωm_z)^2
    d * Dv0 * (1 + 0.763 * fν)
end

# --- linear-spectrum shape quantities ----------------------------------------

"σ(R) on the linear spectrum, R in Mpc (top-hat). Thin wrapper over `σ_R`."
_sigma(P::MatterPowerSpectrum, R) = σ_R(P, R)

"Displacement variance σ_v² = (1/3)∫ Δ²_lin(k)/k² dln k, in Mpc²."
function _sigma_v2(P::MatterPowerSpectrum)
    lo, hi = first(P.logk), last(P.logk)
    quadgk(lk -> (k = exp(lk); dimensionless_power(P, k) / k^2), lo, hi; rtol=1e-6)[1] / 3
end

"Effective index n_eff = -3 - dlnσ²/dlnR at the collapse scale R (top-hat)."
function _neff_at(P::MatterPowerSpectrum, R)
    h = 1e-3
    lnσ2(lr) = 2 * log(_sigma(P, exp(lr)))
    -3 - (lnσ2(log(R) + h) - lnσ2(log(R) - h)) / (2h)
end

# --- Eisenstein-Hu no-wiggle transfer, for BAO de-wiggling -------------------

"""
Eisenstein & Hu (1998) no-wiggle transfer function T_nw(k), k in 1/Mpc. Used
only to define the smooth broadband against which the BAO wiggle is isolated.
"""
function _eh_nowiggle(c::Cosmology, k)
    h = c.h
    Ωm = (Ω_b(c) + Ω_c(c))
    Ωb = Ω_b(c)
    ωm = Ωm * h^2
    fb = Ωb / Ωm
    θ = c.Tcmb / 2.7
    # sound horizon (fitting form) and shape parameter
    s = 44.5 * log(9.83 / ωm) / sqrt(1 + 10 * (Ωb * h^2)^0.75)   # Mpc
    αΓ = 1 - 0.328 * log(431 * ωm) * fb + 0.38 * log(22.3 * ωm) * fb^2
    kh = k                                   # already 1/Mpc
    Γeff = Ωm * h * (αΓ + (1 - αΓ) / (1 + (0.43 * kh * s)^4))
    q = k * θ^2 / (Γeff * h)                  # note Γeff carries an h
    L = log(2 * ℯ + 1.8 * q)
    Cq = 14.2 + 731 / (1 + 62.5 * q)
    L / (L + Cq * q^2)
end

"""
De-wiggled linear P(k): split into an EH-no-wiggle smooth part and the residual
wiggle (isolated by Gaussian-smoothing ln[P/P_nw] in ln k with σ=0.25), then
damp the wiggle by e^{-k²σ_v²}. Returns a function k ↦ P_dwl(k).
"""
function _dewiggle(P::MatterPowerSpectrum, c::Cosmology, σv2)
    lo, hi = first(P.logk), last(P.logk)
    # ratio r(lnk) = ln[ P_lin / (k^ns T_nw²) ]  (broadband divided out)
    ns = c.n_s
    ratio(lk) = begin
        k = exp(lk)
        log(power(P, k) / (k^ns * _eh_nowiggle(c, k)^2))
    end
    # Gaussian smoothing of the ratio in ln k, σ = 0.25 (Mead §4.1 footnote).
    σsm = 0.25
    smooth(lk) = begin
        num = quadgk(u -> ratio(u) * exp(-(u - lk)^2 / (2σsm^2)), max(lo, lk - 4σsm), min(hi, lk + 4σsm); rtol=1e-4)[1]
        den = quadgk(u -> exp(-(u - lk)^2 / (2σsm^2)), max(lo, lk - 4σsm), min(hi, lk + 4σsm); rtol=1e-4)[1]
        num / den
    end
    function (k)
        lk = log(k)
        Psm = exp(smooth(lk)) * k^ns * _eh_nowiggle(c, k)^2   # smooth broadband
        Pw = power(P, k) - Psm                                # the wiggle
        Psm + exp(-k^2 * σv2) * Pw
    end
end

# --- NFW profile --------------------------------------------------------------

"""
Normalised Fourier transform of the truncated NFW profile, W(M,k), with
W(M,0)=1. `c` is the concentration, `rv` the virial radius (Mpc), `k` in 1/Mpc.
"""
function _nfw_window(k, rv, conc)
    rs = rv / conc
    kr = k * rs
    kc = (1 + conc) * kr
    si_kc, ci_kc = sinint(kc), cosint(kc)
    si_k, ci_k = sinint(kr), cosint(kr)
    norm = log(1 + conc) - conc / (1 + conc)
    (sin(kr) * (si_kc - si_k) - sin(conc * kr) / kc +
     cos(kr) * (ci_kc - ci_k)) / norm
end

# --- the model ---------------------------------------------------------------

struct HMcodeSpectrum
    lin::MatterPowerSpectrum
    δc::Float64
    Δv::Float64
    σ8c::Float64
    neff::Float64
    ρ̄m::Float64          # comoving mean matter density, M⊙/Mpc³
    par::NamedTuple
    dewiggled            # k ↦ P_dwl(k)
    Rgrid::Vector{Float64}
    σgrid::Vector{Float64}
    growth::NamedTuple
end

"""
    hmcode_power(P)

Build the HMcode-2020 non-linear spectrum on top of a linear
[`MatterPowerSpectrum`](@ref) `P` (which must reach kmax ≳ 10/Mpc). Evaluate
with [`power`](@ref): `power(hm, k)` returns P_nl(k) in Mpc³.
"""
function hmcode_power(P::MatterPowerSpectrum)
    c = P.cosmo
    z = P.z
    a = scale_factor(z)
    E2 = E(c, a)^2
    Ωm_z = sum(ρ_over_ρc0(s, a) for s in c.species
               if (s isa Baryons || s isa ColdDarkMatter || s isa MassiveNeutrinos); init=0.0) / E2
    Ωm0 = Ω_b(c) + Ω_c(c) + _Ω_or_zero(c, MassiveNeutrinos)
    fν = Ωm0 > 0 ? _Ω_or_zero(c, MassiveNeutrinos) / Ωm0 : 0.0

    gr = _hmcode_growth(c, z)
    x = gr.g / a
    y = gr.Gacc / a
    δc = _dc_mead(x, y, Ωm_z, fν)
    Δv = _dv_mead(x, y, Ωm_z, fν)

    σ8c = _sigma(P, 8.0 / c.h)                       # cold σ8 (=total, massless ν)
    # collapse scale: σ(R_nl) = δc, then n_eff there
    Rnl = find_zero(lr -> _sigma(P, exp(lr)) - δc, (log(1e-3), log(50.0)), Bisection(); xatol=1e-5)
    Rnl = exp(Rnl)
    neff = _neff_at(P, Rnl)

    # fitted parameters (k in h/Mpc → convert with h). Clamps per CAMB.
    h = c.h
    kd = 0.05699 * σ8c^(-1.089) * h
    fdamp = clamp(0.2696 * σ8c^0.9403, 1e-3, 0.99)
    nd = 2.853
    kstar = 0.05618 * σ8c^(-1.013) * h
    η = 0.1281 * σ8c^(-0.3644)
    B = 5.196
    α = clamp(1.875 * 1.603^neff, 0.5, 2.0)
    par = (; kd, fdamp, nd, kstar, η, B, α)

    ρ̄m = Ωm0 * _RHO_CRIT_MSUN_MPC3 * h^2             # M⊙/Mpc³ comoving
    σv2 = _sigma_v2(P)
    dewig = _dewiggle(P, c, σv2)

    # σ(R) table for the mass integral, R ∈ [1e-3, 50] Mpc.
    Rgrid = exp.(range(log(1e-3), log(50.0); length=256))
    σgrid = [_sigma(P, R) for R in Rgrid]

    HMcodeSpectrum(P, δc, Δv, σ8c, neff, ρ̄m, par, dewig, Rgrid, σgrid, gr)
end

"σ(R) from the cached table (log-log linear interp), R in Mpc."
function _sigma_tab(hm::HMcodeSpectrum, R)
    lr = log(R)
    lo, hi = log(hm.Rgrid[1]), log(hm.Rgrid[end])
    lr = clamp(lr, lo, hi)
    i = searchsortedlast(hm.Rgrid, exp(lr))
    i = clamp(i, 1, length(hm.Rgrid) - 1)
    t = (lr - log(hm.Rgrid[i])) / (log(hm.Rgrid[i+1]) - log(hm.Rgrid[i]))
    exp((1 - t) * log(hm.σgrid[i]) + t * log(hm.σgrid[i+1]))
end

"Concentration c(M,z) via Bullock+2001 formation redshift, with Dolag factor."
function _concentration(hm::HMcodeSpectrum, M)
    # z_f from g(z_f)/g(z) · σ(γM) = δc, γ = 0.01. In practice σ(γM) ≥ δc means
    # z_f ≥ z; use the growth ratio at formation. We solve for the growth value
    # g_f = g(z) · δc / σ(γM_L), then map back through the (monotone) growth.
    γ = 0.01
    R_of_M(m) = (3m / (4π * hm.ρ̄m))^(1 / 3)
    σγ = _sigma_tab(hm, R_of_M(γ * M))
    # Bullock: g(z_f) = g(z)·δc/σ(γM) (capped so z_f ≥ z, i.e. g_f ≤ g_z)
    g_z = hm.growth.g
    g_f = g_z * hm.δc / σγ
    zf_factor = min(g_f / g_z, 1.0)      # = (1+z)/(1+z_f) proxy via growth ~ a
    # c = B (1+z_f)/(1+z). With g ∝ a in MD, (1+z_f)/(1+z) ≈ g_z/g_f.
    conc = hm.par.B / zf_factor
    conc * hm.growth.ginf_ratio
end

"""
    power(hm, k)

HMcode-2020 non-linear matter power spectrum P(k) in Mpc³, k in 1/Mpc.
"""
function power(hm::HMcodeSpectrum, k)
    P = hm.lin
    # --- two-halo: de-wiggled linear, perturbatively damped (eq 16) ----------
    xkd = (k / hm.par.kd)^hm.par.nd
    P2H = hm.dewiggled(k) * (1 - hm.par.fdamp * xkd / (1 + xkd))

    # --- one-halo: mass integral of W², with k⁴ large-scale damping ----------
    # integrate over ln M; ν = δc/σ(M), f_ST(ν)dν = f_ST(ν)(dν/dlnM)dlnM.
    p, q, A = 0.3, 0.707, 0.21616
    R_of_M(M) = (3M / (4π * hm.ρ̄m))^(1 / 3)
    integrand(lnM) = begin
        M = exp(lnM)
        R = R_of_M(M)
        σ = _sigma_tab(hm, R)
        ν = hm.δc / σ
        # dν/dlnM = -ν dlnσ/dlnM
        dlnσ_dlnM = begin
            hM = 1e-3
            (log(_sigma_tab(hm, R_of_M(M * exp(hM)))) -
             log(_sigma_tab(hm, R_of_M(M * exp(-hM))))) / (2hM)
        end
        dν_dlnM = -ν * dlnσ_dlnM
        fST = A * (1 + (q * ν^2)^(-p)) * exp(-q * ν^2 / 2)
        rv = (3M / (4π * hm.Δv * hm.ρ̄m))^(1 / 3)
        conc = _concentration(hm, M)
        W = _nfw_window(k * ν^hm.par.η, rv, conc)     # η-bloating
        (M / hm.ρ̄m) * W^2 * fST * dν_dlnM
    end
    P1H_raw = quadgk(integrand, log(1e2), log(1e18); rtol=1e-4)[1]
    xks = (k / hm.par.kstar)^4
    P1H = P1H_raw * xks / (1 + xks)

    # --- α-smoothed sum (eq 23) ----------------------------------------------
    Δ2_2H = P2H * k^3 / (2π^2)
    Δ2_1H = P1H * k^3 / (2π^2)
    Δ2 = (Δ2_2H^hm.par.α + Δ2_1H^hm.par.α)^(1 / hm.par.α)
    Δ2 * 2π^2 / k^3
end

dimensionless_power(hm::HMcodeSpectrum, k) = power(hm, k) * k^3 / (2π^2)
