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

Cold vs total matter: the whole halo model — σ, the mass function, the one- and
two-halo terms — runs on the *cold* (cdm+baryon) spectrum, as HMcode-2020
prescribes (neutrinos do not cluster into haloes). The linear input `P` is
Cosmic's `matter_power_spectrum`, already cold (`δ_matter_comoving` sums
cdm+baryon only). For massive ν the total-matter nonlinear spectrum is then
reconstructed as `P_cb^NL · (P_mm^lin/P_cb^lin)` — the neutrinos stay linear and
carry that suppression through unchanged — by passing `P_total =
matter_power_spectrum(...; total = true)` to `hmcode_power`. Validated against
CAMB `mead2020`: cold prescription <0.5% for k ≲ 1 h/Mpc, and the total-matter
reconstruction to <0.1% across 0.05 < k[h] < 8 for Σmν = 0.15 eV. This replaced
the earlier crash (`hmcode_power` previously threw for every massive-ν cosmology
via a `.Ω`-field mismatch; fixed).
"""

using SpecialFunctions: sinint, cosint
using QuadGK: quadgk
using Roots: find_zero, Bisection

const _RHO_CRIT_MSUN_MPC3 = 2.77536627e11   # ρ_crit,0 / h²  [M⊙ Mpc⁻³]

# --- growth (Mead normalisation: g → a in matter domination) ------------------

# Scale-independent linear growth for a background given by E²(a) and the
# matter fraction Ω_m(a). RK4 in x = ln a for
#     d²g/dlna² + (2 + dlnH/dlna) dg/dlna - (3/2)Ω_m(a) g = 0,
# seeded g = a (growing mode) at a = 1e-2, deep in matter domination. Returns
# g(a_target), dg/dlna, and G = ∫₀ᵃ (g/a')da'.
function _growth_solve(E2, Ωm, a_target; n=512)
    dlnH(a) = (log(E2(a * exp(1e-5))) - log(E2(a * exp(-1e-5)))) / (4e-5)
    rhs(g, gp, a) = (gp, -(2 + dlnH(a)) * gp + 1.5 * Ωm(a) * g)
    xi = log(1e-2)
    xf = log(a_target)
    g, gp, G = exp(xi), exp(xi), exp(xi)
    dx = (xf - xi) / n
    for i in 1:n
        x = xi + (i - 1) * dx
        k1g, k1p = rhs(g, gp, exp(x))
        k2g, k2p = rhs(g + dx / 2 * k1g, gp + dx / 2 * k1p, exp(x + dx / 2))
        k3g, k3p = rhs(g + dx / 2 * k2g, gp + dx / 2 * k2p, exp(x + dx / 2))
        k4g, k4p = rhs(g + dx * k3g, gp + dx * k3p, exp(x + dx))
        gmid = g + dx / 2 * k1g
        g += dx / 6 * (k1g + 2k2g + 2k3g + k4g)
        gp += dx / 6 * (k1p + 2k2p + 2k3p + k4p)
        G += dx * gmid
    end
    (g, gp, G)
end

# Same growth ODE, returning the full monotone trajectory (a, g) so the
# formation-redshift concentration can invert g(a). g is Mead-normalised
# (g → a in matter domination), not g(0) = 1.
function _growth_traj(E2, Ωm, a_target; n=512)
    dlnH(a) = (log(E2(a * exp(1e-5))) - log(E2(a * exp(-1e-5)))) / (4e-5)
    rhs(g, gp, a) = (gp, -(2 + dlnH(a)) * gp + 1.5 * Ωm(a) * g)
    xi = log(1e-2)
    xf = log(a_target)
    dx = (xf - xi) / n
    as = Vector{Float64}(undef, n + 1)
    gs = Vector{Float64}(undef, n + 1)
    g, gp = exp(xi), exp(xi)
    as[1], gs[1] = exp(xi), g
    for i in 1:n
        x = xi + (i - 1) * dx
        k1g, k1p = rhs(g, gp, exp(x))
        k2g, k2p = rhs(g + dx / 2 * k1g, gp + dx / 2 * k1p, exp(x + dx / 2))
        k3g, k3p = rhs(g + dx / 2 * k2g, gp + dx / 2 * k2p, exp(x + dx / 2))
        k4g, k4p = rhs(g + dx * k3g, gp + dx * k3p, exp(x + dx))
        g += dx / 6 * (k1g + 2k2g + 2k3g + k4g)
        gp += dx / 6 * (k1p + 2k2p + 2k3p + k4p)
        as[i+1] = exp(x + dx)
        gs[i+1] = g
    end
    (as, gs)
end

"""
    _hmcode_growth(c, z)

Returns `(g, Gacc, ginf_ratio)`: the growth `g = g(z)` (g → a in matter
domination), the accumulated growth `Gacc = ∫₀ᵃ g/a' da'`, and the Dolag
dark-energy concentration factor `ginf_ratio`.

The Dolag factor is the HMcode-2020 (CAMB `halofit.f90`, imead==3) form:
`(ĝ_de(z_c)/ĝ_ΛCDM(z_c)) · (ĝ_ΛCDM(z)/ĝ_de(z))` at collapse redshift
`z_c = 10`, with each growth `ĝ` normalised to unity today in its own
cosmology. The ΛCDM reference is the flat, `w = -1` model with the same
matter and radiation content (Λ set for flatness), solved with the same
growth ODE — not the placeholder `1.0` that stood here before. It reduces
to 1 identically for flat ΛCDM (the two backgrounds coincide) and for
`z ≥ z_c`.
"""
function _hmcode_growth(c::Cosmology, z)
    a_target = scale_factor(z)
    Ωm_a(a) = begin
        s = 0.0
        for sp in c.species
            (sp isa Baryons || sp isa ColdDarkMatter || sp isa MassiveNeutrinos) &&
                (s += ρ_over_ρc0(sp, a))
        end
        s / E(c, a)^2
    end
    E2m(a) = E(c, a)^2
    as, gs = _growth_traj(E2m, Ωm_a, a_target)
    g, G = gs[end], 0.0
    G = _growth_solve(E2m, Ωm_a, a_target)[3]
    # invert g(a) for the formation-redshift concentration; gs is monotone
    a_of_g = linear_interpolation(gs, as, extrapolation_bc=Line())

    # Dolag correction only for models whose growth history departs from ΛCDM.
    has_de = any(s -> (s isa AbstractDarkEnergy && !(s isa CosmologicalConstant)) ||
                      s isa Curvature, c.species)
    ginf_ratio = 1.0
    z_c = 10.0
    if has_de && z < z_c
        # matter + radiation densities today (ν counted with matter, as in the
        # HMcode reference) and the flat ΛCDM-equivalent background.
        Ωm0 = sum(ρ_over_ρc0(s, 1.0) for s in c.species
                  if s isa Baryons || s isa ColdDarkMatter || s isa MassiveNeutrinos)
        Ωr0 = sum(ρ_over_ρc0(s, 1.0) for s in c.species
                  if s isa Photons || s isa MasslessNeutrinos)
        ΩΛ0 = 1 - Ωm0 - Ωr0
        E2L(a) = Ωr0 / a^4 + Ωm0 / a^3 + ΩΛ0
        ΩmL(a) = (Ωm0 / a^3) / E2L(a)
        a_c = 1 / (1 + z_c)
        gde_z = g
        gde_zc, _, _ = _growth_solve(E2m, Ωm_a, a_c)
        gde_0, _, _ = _growth_solve(E2m, Ωm_a, 1.0)
        gL_z, _, _ = _growth_solve(E2L, ΩmL, a_target)
        gL_zc, _, _ = _growth_solve(E2L, ΩmL, a_c)
        gL_0, _, _ = _growth_solve(E2L, ΩmL, 1.0)
        # each ĝ normalised to today, per the CAMB convention
        ĝde_zc = gde_zc / gde_0; ĝL_zc = gL_zc / gL_0
        ĝde_z = gde_z / gde_0; ĝL_z = gL_z / gL_0
        ginf_ratio = (ĝde_zc / ĝL_zc) * (ĝL_z / ĝde_z)
    end
    (g=g, Gacc=G, ginf_ratio=ginf_ratio, a_of_g=a_of_g, a_z=a_target)
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
    lin::MatterPowerSpectrum       # cold (cdm+baryon) linear P — drives the halo model
    lin_total::MatterPowerSpectrum # total (cdm+baryon+ν) linear P — for the output ratio
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
    hmcode_power(P; P_total = P)

Build the HMcode-2020 non-linear spectrum on top of a linear *cold* (cdm+baryon)
[`MatterPowerSpectrum`](@ref) `P` (which must reach kmax ≳ 10/Mpc). Evaluate
with [`power`](@ref): `power(hm, k)` returns P_nl(k) in Mpc³.

The entire halo model — σ(R), the mass function, the one- and two-halo terms —
runs on the cold spectrum `P`, exactly as HMcode-2020 prescribes (neutrinos do
not cluster into haloes). For a massive-ν cosmology the *total*-matter nonlinear
spectrum is then reconstructed from the cold halo-model result by the linear
total/cold ratio,

    P_mm(k) = P_cb^NL(k) · [P_mm^lin(k) / P_cb^lin(k)],

because the neutrinos stay linear and carry that scale-dependent suppression
(→ f_cb² at high k) through the nonlinear evolution unchanged. Pass `P_total`
from `matter_power_spectrum(...; total = true)` to activate it; for massless ν
`P_total = P` makes the ratio 1 and the result is the pure cold spectrum.
Validated against CAMB `mead2020` (Σmν = 0.15 eV) to <0.1% across 0.05 < k[h] < 8.
"""
function hmcode_power(P::MatterPowerSpectrum; P_total::MatterPowerSpectrum=P)
    c = P.cosmo
    z = P.z
    a = scale_factor(z)
    E2 = E(c, a)^2
    Ωm_z = sum(ρ_over_ρc0(s, a) for s in c.species
               if (s isa Baryons || s isa ColdDarkMatter || s isa MassiveNeutrinos); init=0.0) / E2
    # massive-ν rest-mass density today (ρ/ρ_c0 at a=1) — NOT `_Ω_or_zero`,
    # which reads a `.Ω` field the MassiveNeutrinos struct does not carry
    # (it stores Ω_rel, the would-be-massless density); that mismatch made
    # hmcode_power throw for every massive-ν cosmology.
    Ω_mν = sum(ρ_over_ρc0(s, 1.0) for s in c.species
               if s isa MassiveNeutrinos; init=0.0)
    Ωm0 = Ω_b(c) + Ω_c(c) + Ω_mν
    fν = Ωm0 > 0 ? Ω_mν / Ωm0 : 0.0

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

    HMcodeSpectrum(P, P_total, δc, Δv, σ8c, neff, ρ̄m, par, dewig, Rgrid, σgrid, gr)
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
    # Formation redshift z_f from the collapse condition g(a_f) = g(z)·δc/σ(γM),
    # γ = 0.01 (CAMB halofit.f90 `zcoll_bull`, imead==3). Then
    # c = B·(1+z_f)/(1+z) = B·a_z/a_f. Inverting the growth for a_f — rather
    # than the old g_f/g_z proxy, which assumes g ∝ a and so silently drops the
    # g(z)≈0.78 growth-suppression factor — matters because the concentration
    # sets the NFW window only at HIGH k (W→1 at low k for any c). The proxy
    # under-concentrated small haloes ~1/g(z) ≈ 28%, and that was exactly the
    # high-k P(k) deficit vs CAMB (~9% at k≈8).
    γ = 0.01
    R_of_M(m) = (3m / (4π * hm.ρ̄m))^(1 / 3)
    σγ = _sigma_tab(hm, R_of_M(γ * M))
    g_z = hm.growth.g
    g_f = g_z * hm.δc / σγ
    # g_f ≥ g_z (σγ ≤ δc) ⇒ no earlier collapse, z_f = z; else invert g(a_f)=g_f
    zf_ratio = g_f >= g_z ? 1.0 : hm.growth.a_z / hm.growth.a_of_g(g_f)
    hm.par.B * zf_ratio * hm.growth.ginf_ratio
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
    P_cb = Δ2 * 2π^2 / k^3
    # cold → total: multiply by the linear total/cold power ratio. δ(k) here is
    # the ℛ-normalised transfer, so the primordial factor cancels and the ratio
    # is (δ_total/δ_cold)². Identically 1 for massless ν (lin_total === lin).
    hm.lin_total === hm.lin && return P_cb
    r = transfer(hm.lin_total, k) / transfer(hm.lin, k)
    P_cb * r^2
end

dimensionless_power(hm::HMcodeSpectrum, k) = power(hm, k) * k^3 / (2π^2)
