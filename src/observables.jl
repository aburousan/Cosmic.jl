"""
Observables built from the perturbation hierarchy: transfer functions, the
matter power spectrum, σ₈, and the growth factor.

The Boltzmann solver is run with the primordial curvature perturbation set to
ℛ = 1 for every mode, so what it returns *is* the transfer function. The system
is linear, so the primordial spectrum

    P_ℛ(k) = A_s (k/k_pivot)^(n_s - 1) · 2π²/k³

is applied afterwards rather than being carried through the integration. That is
worth doing deliberately: it means changing A_s or n_s costs nothing, and only a
change to the *background* forces a re-solve.
"""

using QuadGK: quadgk
using Interpolations: linear_interpolation, Line

"""
    primordial_power(c, k)

Dimensionless primordial curvature power spectrum Δ²_ℛ(k). By default the power
law Δ²_ℛ(k) = A_s (k/k_pivot)^(n_s − 1), but if the cosmology carries a
`primordial` override (any callable `k ↦ Δ²_ℛ(k)`, e.g. an inflaton solution or a
tabulated broken-scale-invariance spectrum) that is used instead. The override
lets features, steps, running, or a computed inflationary spectrum enter without
touching the transfer-function machinery, since the spectrum is applied linearly
on top of the ℛ = 1 transfer functions.
"""
primordial_power(c::Cosmology, k) = _primordial_power(c.primordial, c, k)
_primordial_power(::Nothing, c::Cosmology, k) = c.A_s * (k / c.k_pivot)^(c.n_s - 1)
_primordial_power(p, c::Cosmology, k) = p(k)

"""
    MatterPowerSpectrum

Tabulated transfer function δ_m(k) at one redshift, ready to be turned into P(k).
Built by [`matter_power_spectrum`](@ref).
"""
struct MatterPowerSpectrum{C,I}
    cosmo::C
    z::Float64
    logk::Vector{Float64}
    δm::Vector{Float64}      # transfer function, ℛ = 1 normalisation
    interp::I
end

"""
    bao_k_grid(rec, kmin, kmax; per_decade = 10, per_decade_bao = 70,
               bao_center = 3.0, bao_width = 4.0)

Wavenumber grid with the sampling density boosted across the BAO region.

A uniform grid in ln k is the obvious choice and it is wrong, quietly. The baryon
acoustic oscillations in P(k) are periodic in *linear* k with period
Δk = 2π/r_s ≈ 0.04 /Mpc, so a log grid that comfortably resolves the turnover at
k_eq gives only two or three points per wiggle by k ~ 0.2, and the interpolated
spectrum comes out visibly kinked -- the wiggles get aliased, not smoothed.

CLASS solves this by varying the density (`perturbations.c`, and the defaults in
`precisions.h`): 10 points per decade in general, rising to 70 per decade inside a
window centred on 3·k_rec with logarithmic width 4, where k_rec = 2π/r_s(z_*). That
window is exactly where the oscillations live, so the cost is paid only where it
buys something.

    dln k = ln10 / [ n_pk + (n_bao - n_pk)·(1 - tanh( ((ln k - ln(3 k_rec))/ln 4)^4 )) ]
"""
function bao_k_grid(rec::RecombinationSolution, kmin, kmax;
    per_decade=10.0, per_decade_bao=70.0, bao_center=3.0, bao_width=4.0)

    k_rec = 2π / sound_horizon(rec.cosmo, z_star(rec))
    lnc = log(bao_center * k_rec)
    lnw = log(bao_width)

    ks = [float(kmin)]
    k = float(kmin)
    while k < kmax
        u = (log(k) - lnc) / lnw
        dens = per_decade + (per_decade_bao - per_decade) * (1 - tanh(u^4))
        k *= 10^(1 / dens)
        push!(ks, k)
    end
    ks
end

"""
    matter_power_spectrum(c, rec; z = 0, kmin = 1e-4, kmax = 10, nk = nothing)

Solve the hierarchy on a grid of wavenumbers and tabulate the matter transfer
function at redshift `z`. `k` is in 1/Mpc.

With `nk = nothing` (the default) the grid comes from [`bao_k_grid`](@ref), which
refines itself across the acoustic oscillations. Passing an integer `nk` forces a
plain log-spaced grid of that size instead -- simpler, faster, and liable to alias
the BAO wiggles, so it is for quick looks rather than for results.

One Boltzmann solve per mode: this is the expensive step, and it is embarrassingly
parallel, so the modes are threaded.
"""
function matter_power_spectrum(c::Cosmology, rec::RecombinationSolution;
    z=0.0, kmin=1e-4, kmax=10.0, nk=nothing, lmax_γ=25, lmax_ν=32)

    bg = BackgroundCache(c, rec)
    logk = nk === nothing ? log.(bao_k_grid(rec, kmin, kmax)) :
           collect(range(log(kmin), log(kmax); length=nk))
    δm = zeros(length(logk))
    a = scale_factor(z)

    Threads.@threads for i in eachindex(logk)
        k = exp(logk[i])
        p = solve_perturbations(c, bg, k; lmax_γ, lmax_ν)
        δm[i] = δ_matter_comoving(p, a)
    end

    MatterPowerSpectrum(c, float(z), logk, δm,
        linear_interpolation(logk, δm; extrapolation_bc=Line()))
end

"""
    transfer(P, k)

Matter transfer function δ_m(k) at the spectrum's redshift, ℛ = 1 normalised.
"""
transfer(P::MatterPowerSpectrum, k) = P.interp(log(k))

"""
    power(P, k)

Matter power spectrum P(k) in Mpc³:

    P(k) = δ_m(k)² · A_s (k/k_pivot)^(n_s-1) · 2π²/k³
"""
function power(P::MatterPowerSpectrum, k)
    δ = transfer(P, k)
    δ^2 * primordial_power(P.cosmo, k) * 2π^2 / k^3
end

"""
    dimensionless_power(P, k)

Δ²(k) = k³P(k)/(2π²) = δ_m(k)² Δ²_ℛ(k). The variance of the matter field per
logarithmic interval in k, which is the form σ₈ is an integral of.
"""
function dimensionless_power(P::MatterPowerSpectrum, k)
    transfer(P, k)^2 * primordial_power(P.cosmo, k)
end

"""
    top_hat_window(x)

Fourier transform of a real-space spherical top hat,
W(x) = 3(sin x - x cos x)/x³, with W(0) = 1.
"""
function top_hat_window(x)
    x < 1e-4 && return 1 - x^2 / 10          # the naive form loses all precision here
    3 * (sin(x) - x * cos(x)) / x^3
end

"""
    σ_R(P, R)

RMS linear matter fluctuation in spheres of comoving radius `R` (Mpc):

    σ_R² = ∫ dln k  Δ²(k)  W(kR)²
"""
function σ_R(P::MatterPowerSpectrum, R; rtol=1e-6)
    integrand(logk) = begin
        k = exp(logk)
        dimensionless_power(P, k) * top_hat_window(k * R)^2
    end
    sqrt(quadgk(integrand, first(P.logk), last(P.logk); rtol)[1])
end

"""
    σ8(P)

σ₈: the RMS matter fluctuation in spheres of radius 8 h⁻¹ Mpc. Note the radius is
in h⁻¹ Mpc while the code works in Mpc, so R = 8/h.
"""
σ8(P::MatterPowerSpectrum) = σ_R(P, 8.0 / P.cosmo.h)

# --- Growth -----------------------------------------------------------------

"""
    growth_factor(c, rec, z; k = 1e-3)

Linear growth factor D(z), normalised to D(0) = 1.

Measured from an actual mode rather than from the usual scale-independent growth
ODE. That matters: with massive neutrinos (or any component that free-streams)
growth *is* scale dependent, and the textbook ODE quietly assumes it is not. The
default k = 1e-3 /Mpc is large-scale enough to recover the scale-independent
answer in plain ΛCDM.
"""
function growth_factor(c::Cosmology, rec::RecombinationSolution, z; k=1e-3)
    bg = BackgroundCache(c, rec)
    p = solve_perturbations(c, bg, k)
    δ_matter(p, scale_factor(z)) / δ_matter(p, 1.0)
end

"""
    growth_rate(c, rec, z; k = 1e-3)

f = dln D / dln a, the quantity redshift-space distortions measure. In ΛCDM it is
close to Ω_m(a)^0.55.
"""
function growth_rate(c::Cosmology, rec::RecombinationSolution, z; k=1e-3, h=1e-3)
    bg = BackgroundCache(c, rec)
    p = solve_perturbations(c, bg, k)
    a = scale_factor(z)
    x = log(a)
    d(xx) = log(abs(δ_matter(p, exp(xx))))
    (d(x + h) - d(x - h)) / (2h)
end
