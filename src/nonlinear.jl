"""
Non-linear matter power spectrum.

The linear P(k) from the Boltzmann solve is only the starting point: below
k ~ 0.1/Mpc gravitational collapse moves power around in a way linear theory
cannot follow. This module adds the standard fitting-function corrections on
top of the linear spectrum.

Two prescriptions are provided, sharing one piece of machinery -- the
non-linear scale. Both need the Gaussian-smoothed variance

    σ²(R) = ∫ Δ²_lin(k) e^{-k²R²} dln k

the radius R_nl where σ = 1 (with k_σ = 1/R_nl the non-linear wavenumber), and
the shape of the spectrum there, captured by the effective slope and curvature

    n_eff = -3 - dlnσ²/dlnR,   C = -d²lnσ²/dlnR²    at R = R_nl.

`halofit_power` is Takahashi et al. (2012) with the Bird et al. (2012) massive-ν
terms -- the coefficients are transcribed from CLASS 3.3.4's `halofit.c`, so it
reproduces CLASS's own `pk_nl` and serves as the validation baseline.
`hmcode_power` (in `hmcode.jl`) is the more accurate halo-model calculation and
is the default for serious use.
"""

using QuadGK: quadgk
using Roots: find_zero, Bisection

# --- the non-linear scale (shared) -------------------------------------------

"""
    _gauss_moments(P, R)

The three CLASS halofit integrals at smoothing radius `R`:

    s1 = ∫ Δ²(k) e^{-k²R²} dln k                 (= σ²_Gauss(R))
    s2 = ∫ Δ²(k) e^{-k²R²} · 2x²    dln k         x² = k²R²
    s3 = ∫ Δ²(k) e^{-k²R²} · 4x²(1-x²) dln k

from which dlnσ²/dlnR = -s2/s1 and d²lnσ²/dlnR² = -(s2/s1)² - s3/s1.
"""
function _gauss_moments(P::MatterPowerSpectrum, R; rtol=1e-6)
    lo, hi = first(P.logk), last(P.logk)
    f(lk) = begin
        k = exp(lk)
        x2 = k^2 * R^2
        d = dimensionless_power(P, k) * exp(-x2)
        (d, d * 2x2, d * 4x2 * (1 - x2))
    end
    s1 = quadgk(lk -> f(lk)[1], lo, hi; rtol)[1]
    s2 = quadgk(lk -> f(lk)[2], lo, hi; rtol)[1]
    s3 = quadgk(lk -> f(lk)[3], lo, hi; rtol)[1]
    s1, s2, s3
end

"""
    nonlinear_scale(P)

Return `(k_σ, n_eff, C)` — the non-linear wavenumber (1/Mpc), the effective
slope, and the curvature of the linear spectrum at the non-linear scale. `k_σ`
is found from σ_Gauss(1/k_σ) = 1 by bisection in log R, exactly as CLASS does.
Throws if the linear spectrum does not reach far enough into the non-linear
regime (σ never reaches 1 over the available k-range).
"""
function nonlinear_scale(P::MatterPowerSpectrum)
    σ2(R) = _gauss_moments(P, R)[1]
    # σ decreases with R; bracket the root σ(R)=1 in log10 R.
    Rmin = 1 / exp(last(P.logk))       # smallest resolvable smoothing radius
    Rmax = 1 / exp(first(P.logk))
    σ2(Rmin) >= 1 || error("linear P(k) does not reach the non-linear scale; raise kmax")
    σ2(Rmax) <= 1 || error("linear P(k) already non-linear at the largest scale; lower kmin")
    logR = find_zero(lr -> σ2(10.0^lr) - 1, (log10(Rmin), log10(Rmax)), Bisection(); xatol=1e-6)
    R = 10.0^logR
    s1, s2, s3 = _gauss_moments(P, R)
    d1 = -s2 / s1
    d2 = -(s2 / s1)^2 - s3 / s1
    (k_σ=1 / R, n_eff=-3 - d1, C=-d2)
end

"Matter/radiation/DE density parameters and w, fν at the spectrum's redshift."
function _halofit_background(c::Cosmology, z)
    a = scale_factor(z)
    E2 = E(c, a)^2
    ismatter(s) = s isa Baryons || s isa ColdDarkMatter || s isa MassiveNeutrinos
    israd(s) = s isa Photons || s isa MasslessNeutrinos
    Ω_m_z = sum(ρ_over_ρc0(s, a) for s in c.species if ismatter(s); init=0.0) / E2
    Ω_r_z = sum(ρ_over_ρc0(s, a) for s in c.species if israd(s); init=0.0) / E2
    Ω_v_z = 1 - Ω_m_z - Ω_r_z                       # CLASS convention (flat: = Ω_Λ(z))
    Ω_m0 = Ω_m(c) + _Ω_or_zero(c, MassiveNeutrinos)
    fν = Ω_m0 > 0 ? _Ω_or_zero(c, MassiveNeutrinos) / Ω_m0 : 0.0
    w_de = -1.0
    for s in c.species
        s isa AbstractDarkEnergy && (w_de = w(s, a))
    end
    (Ω_m=Ω_m_z, Ω_v=Ω_v_z, w=w_de, fν=fν)
end

# --- Takahashi 2012 (+ Bird 2012) halofit ------------------------------------

"""
    HalofitSpectrum

Non-linear P(k) from Takahashi et al. (2012). Callable via [`power`](@ref):
`power(hf, k)` returns the non-linear P(k) in Mpc³. Also carries the linear
spectrum it was built from and the non-linear-scale diagnostics.
"""
struct HalofitSpectrum
    lin::MatterPowerSpectrum
    k_σ::Float64
    n_eff::Float64
    C::Float64
    coef::NamedTuple
end

"""
    halofit_power(P)

Build the Takahashi 2012 + Bird 2012 non-linear spectrum on top of a linear
[`MatterPowerSpectrum`](@ref) `P`. The linear spectrum must reach well past the
non-linear scale (kmax ≳ 10/Mpc) so the Gaussian variance integral converges.
"""
function halofit_power(P::MatterPowerSpectrum)
    ns = nonlinear_scale(P)
    n, Cc = ns.n_eff, ns.C
    bg = _halofit_background(P.cosmo, P.z)
    Ωm, Ωv, w, fν = bg.Ω_m, bg.Ω_v, bg.w, bg.fν

    # Takahashi 2012 eqs A6–A13, with the Bird 2012 β term and neutrino factors,
    # transcribed from CLASS 3.3.4 external/Halofit/halofit.c.
    gam = 0.1971 - 0.0843n + 0.8460Cc
    a = 10.0^(1.5222 + 2.8553n + 2.3706n^2 + 0.9903n^3 + 0.2250n^4 - 0.6038Cc + 0.1749Ωv * (1 + w))
    b = 10.0^(-0.5642 + 0.5864n + 0.5716n^2 - 1.5474Cc + 0.2279Ωv * (1 + w))
    cc = 10.0^(0.3698 + 2.0404n + 0.8161n^2 + 0.5869Cc)
    xμ = 0.0
    xν = 10.0^(5.2105 + 3.6902n)
    α = abs(6.0835 + 1.3373n - 0.1959n^2 - 5.5274Cc)
    β = 2.0379 - 0.7354n + 0.3157n^2 + 1.2490n^3 + 0.3980n^4 - 0.1682Cc + fν * (1.081 + 0.395n^2)

    if abs(1 - Ωm) > 0.01
        frac = Ωv / (1 - Ωm)
        f1 = frac * Ωm^(-0.0307) + (1 - frac) * Ωm^(-0.0732)
        f2 = frac * Ωm^(-0.0585) + (1 - frac) * Ωm^(-0.1423)
        f3 = frac * Ωm^0.0743 + (1 - frac) * Ωm^0.0725
    else
        f1 = f2 = f3 = 1.0
    end

    coef = (; gam, a, b, c=cc, xμ, xν, α, β, f1, f2, f3, fν, h=P.cosmo.h)
    HalofitSpectrum(P, ns.k_σ, n, Cc, coef)
end

"""
    power(hf, k)

Non-linear matter power spectrum P(k) in Mpc³ at wavenumber `k` (1/Mpc).
Below the non-linear scale it returns the linear spectrum; above it, the
Takahashi halo (one-halo) plus quasi-linear (two-halo) sum.
"""
function power(hf::HalofitSpectrum, k)
    P = hf.lin
    Δ2_lin = dimensionless_power(P, k)     # Δ²(k) = k³P(k)/2π²
    y = k / hf.k_σ
    C = hf.coef

    # One-halo term (Takahashi A2–A5), with the Bird 0.977 neutrino factor.
    pk_halo = C.a * y^(3 * C.f1) / (1 + C.b * y^C.f2 + (C.f3 * C.c * y)^(3 - C.gam))
    pk_halo = pk_halo / (1 + C.xμ / y + C.xν / y^2) * (1 + C.fν * 0.977)

    # Two-halo (quasi-linear) term, with the Bird neutrino boost of the linear
    # power that feeds the α/β suppression. kh = k in h/Mpc for the 47.48/1.5.
    kh = k / C.h
    Δ2_aa = Δ2_lin * (1 + C.fν * 47.48 * kh^2 / (1 + 1.5 * kh^2))
    pk_quasi = Δ2_lin * (1 + Δ2_aa)^C.β / (1 + Δ2_aa * C.α) * exp(-y / 4 - y^2 / 8)

    Δ2_nl = pk_halo + pk_quasi
    Δ2_nl * 2π^2 / k^3
end

"Dimensionless non-linear power Δ²(k) = k³P(k)/2π²."
dimensionless_power(hf::HalofitSpectrum, k) = power(hf, k) * k^3 / (2π^2)
