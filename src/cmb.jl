"""
CMB angular power spectra by line-of-sight integration.

The direct approach -- evolve the photon hierarchy up to ℓ = 2000 and read off the
multipoles today -- is hopeless: it needs thousands of coupled equations per mode.
Seljak & Zaldarriaga (1996) observed that the ℓ-dependence can be factored out of
the dynamics entirely. Everything that *sources* anisotropy happens at
recombination and is confined to ℓ ≤ 2; the projection of that source onto today's
sky is pure geometry, carried by a spherical Bessel function:

    Θ_ℓ(k) = ∫₀^{η₀} dη  S(k, η)  j_ℓ(k(η₀ - η))

So we evolve only a short hierarchy, build the source S once, and project. The cost
drops by orders of magnitude, which is what made CMB parameter estimation possible.

The source function has four pieces, each a distinct physical effect:

    S = g·(Θ₀ + ψ + Π/4)              intrinsic temperature + Sachs-Wolfe
      + e^{-τ}·(φ' + ψ')              integrated Sachs-Wolfe
      - (1/k)·d(g·v_b)/dη             Doppler, from the baryon velocity
      + (3/4k²)·d²(g·Π)/dη²           polarization quadrupole

with the visibility g = κ̇e^{-τ} confining the first, third and fourth to the last
scattering surface. The second is not: e^{-τ} is order unity all the way to us, so
any late-time evolution of the potentials contributes -- which in a Λ universe they
do, and that is the late ISW effect.
"""

using SpecialFunctions: sphericalbesselj
using Interpolations: cubic_spline_interpolation, linear_interpolation, Line
using QuadGK: quadgk

"""
    SourceFunction

The line-of-sight sources for one mode k, tabulated on a conformal-time grid.
`S_T` is the temperature source, `S_E` the polarization source (which is just the
visibility-weighted quadrupole, g·Π).
"""
struct SourceFunction
    k::Float64
    η::Vector{Float64}
    S_T::Vector{Float64}
    S_E::Vector{Float64}
    S_L::Vector{Float64}   # lensing-potential source -W(χ)(φ+ψ); empty if unused
end

SourceFunction(k, η, S_T, S_E) = SourceFunction(k, η, S_T, S_E, Float64[])

"""
    source_function(c, bg, p; nη = 4000)

Build S_T and S_E for the solved mode `p`.

The derivatives d(g v_b)/dη and d²(g Π)/dη² are taken numerically on the grid.
That is safe here because g v_b and g Π are smooth and *compactly supported* --
the visibility function kills them away from last scattering -- so there is no
cancellation to amplify and no tail to lose.
"""
function source_function(c::Cosmology, bg::BackgroundCache, p::PerturbationSolution;
    a_min=1e-5, nη=4000, sw=1.0, isw=1.0, dop=1.0, pol=1.0, Wlens=nothing)

    k = p.k
    L = _layout(p)
    rec = p.rec

    # Grid uniform in ln a. The visibility is ~20 Mpc wide in η and sits at
    # η ≈ 280 Mpc, and this grid resolves it with tens of points.
    xs = range(log(a_min), 0.0; length=nη)
    ηs = [bg.η(x) for x in xs]
    η0 = bg.η(0.0)

    gs = zeros(nη)       # visibility g = κ̇ e^{-τ}
    e_τ = zeros(nη)      # e^{-τ}
    gvb = zeros(nη)      # g · v_b
    gΠ = zeros(nη)      # g · Π
    ST0 = zeros(nη)      # the non-derivative part of S_T

    for (i, x) in enumerate(xs)
        a = exp(x)
        z = redshift(a)
        u = p.sol(x)

        κ̇ = bg.κ̇(x)
        τ = optical_depth(rec, z)
        eτ = exp(-τ)
        g = κ̇ * eτ

        Θ0 = u[L.iγ] / 4
        Θ2 = u[L.iγ+2] / 4
        ΘP0 = u[L.iG] / 4
        ΘP2 = u[L.iG+2] / 4
        Π = Θ2 + ΘP0 + ΘP2

        v_b = u[4] / k                      # v_b = θ_b / k

        φ, ψ, φ̇ = _metric_at(c, p, u, L, a)
        # ψ' follows from φ' and the (slowly varying) anisotropic stress; taking
        # it numerically below is simpler and accurate enough for the ISW term,
        # which is a smooth, small contribution.

        gs[i] = g
        e_τ[i] = eτ
        gvb[i] = g * v_b
        gΠ[i] = g * Π
        ST0[i] = g * (Θ0 + ψ + Π / 4)
    end

    # φ + ψ on the grid, for the ISW term (its η-derivative is taken numerically).
    φψ = zeros(nη)
    for (i, x) in enumerate(xs)
        u = p.sol(x)
        φ, ψ, _ = _metric_at(c, p, u, L, exp(x))
        φψ[i] = φ + ψ
    end

    d_gvb = _ddη(gvb, ηs)
    d2_gΠ = _ddη(_ddη(gΠ, ηs), ηs)
    d_φψ = _ddη(φψ, ηs)

    S_T = similar(ηs)
    for i in eachindex(ηs)
        # Each term is a distinct physical effect, and each carries a weight so
        # it can be switched off and its contribution seen in isolation:
        #
        #   sw   g·(Θ₀ + ψ + Π/4)      intrinsic temperature + Sachs-Wolfe
        #   isw  e^{-τ}·(φ' + ψ')      integrated Sachs-Wolfe
        #   dop  (1/k)·d(g·v_b)/dη     Doppler, from the baryon velocity
        #   pol  (3/4k²)·d²(g·Π)/dη²   polarization quadrupole
        #
        # The Doppler sign is +, given v_b = θ_b/k. Callin (2006) writes a minus,
        # but with the opposite velocity convention. It matters enormously and
        # fails quietly: with the sign flipped the acoustic *peaks* still land in
        # the right places (they are set by the geometry, not the interference),
        # while the troughs -- which the Doppler term fills, being π/2 out of
        # phase with the density -- come out up to 4x too deep. The spectrum looks
        # plausible and is wrong.
        S_T[i] = sw * ST0[i] + isw * e_τ[i] * d_φψ[i] +
                 dop * d_gvb[i] / k + pol * 3 * d2_gΠ[i] / (4k^2)
    end
    # The E-mode source is just g·Π. Its ℓ-dependent geometric prefactor belongs
    # at projection time, not here.
    S_E = gΠ

    S_L = Wlens === nothing ? Float64[] :
          [-Wlens(ηs[i]) * φψ[i] for i in 1:nη]

    SourceFunction(k, collect(ηs), S_T, S_E, S_L)
end

"Central-difference derivative of `y` with respect to a non-uniform grid `x`."
function _ddη(y, x)
    n = length(y)
    d = similar(y)
    @inbounds for i in 2:(n-1)
        h1 = x[i] - x[i-1]
        h2 = x[i+1] - x[i]
        d[i] = (y[i+1] - y[i-1]) / (h1 + h2)
    end
    d[1] = (y[2] - y[1]) / (x[2] - x[1])
    d[n] = (y[n] - y[n-1]) / (x[n] - x[n-1])
    d
end

"Metric potentials and φ̇ at one grid point, reusing the perturbation state."
function _metric_at(c::Cosmology, p::PerturbationSolution, u, L, a)
    k = p.k
    H0 = Constants.H0_in_invMpc(c.h)
    ℋ_ = ℋ(c, a)
    ργ = Ω_γ(c) / a^4
    ρν = _Ω_or_zero(c, MasslessNeutrinos) / a^4
    ρb = Ω_b(c) / a^3
    ρc = Ω_c(c) / a^3
    G = MassiveNuGrid(L.nq == 0 ? 1 : L.nq)
    mν = Tuple(get_all_species(c, MassiveNeutrinos))
    dc = _dcdm_params(c)
    ρdc = dc === nothing ? 0.0 : exp(dc.lnρd(log(a)))
    ρdr = dc === nothing ? 0.0 : exp(dc.lnρr(log(a)))
    de = _fld_species(c)
    ρfld = de === nothing ? 0.0 : ρ_over_ρc0(de, a)
    wfld = de === nothing ? -1.0 : w(de, a)
    mg = nothing
    if p.mg !== nothing && log(a) >= p.mg.x_on
        x = log(a)
        αB = p.mg.α_B(x)
        αK = p.mg.α_K(x)
        begin
            M2 = p.mg.M_star2(x)
            rpm = (4 / 3) * (ργ + ρν) + ρb + ρc     # mg forbids dcdm/fld/mν
            E2 = (ℋ_ / (a * H0))^2
            mg = (αB, αK, p.mg.α_M(x), M2, -1.5 * rpm / E2,
                3 * H0^2 * a^2 * rpm / (ℋ_^2 * M2))
        end
    end
    _metric(u, L, k, spatial_curvature_K(c), a, ℋ_, H0, ργ, ρν, ρb, ρc, G, mν,
        ρdc, ρdr, ρfld, wfld, false, 1.0, 0.0, mg)[1:3]
end

# --- Projection --------------------------------------------------------------

"""
    BesselTable

Precomputed j_ℓ(x) on a uniform x grid, one row per ℓ.

The projection integral evaluates j_ℓ at every (ℓ, k, η) triple, which for a real
spectrum is order 1e9 calls -- far too many to compute from scratch. j_ℓ oscillates
with period 2π in x, so a uniform grid with Δx ≈ 0.15 resolves it to well below the
accuracy of everything else here, and linear interpolation off the table costs a
few nanoseconds. This is what makes the whole calculation tractable, and is exactly
what CAMB's `bessels.f90` exists to do.
"""
struct BesselTable
    ℓs::Vector{Int}
    xmax::Float64
    dx::Float64
    tab::Matrix{Float64}    # [iℓ, ix]
end

function BesselTable(ℓs::Vector{Int}, xmax; dx=0.15)
    nx = ceil(Int, xmax / dx) + 2
    tab = zeros(length(ℓs), nx)
    Threads.@threads for i in eachindex(ℓs)
        ℓ = ℓs[i]
        for j in 1:nx
            tab[i, j] = sphericalbesselj(ℓ, (j - 1) * dx)
        end
    end
    BesselTable(ℓs, xmax, dx, tab)
end

@inline function (B::BesselTable)(iℓ::Int, x::Float64)
    x <= 0 && return 0.0
    t = x / B.dx
    j = floor(Int, t) + 1
    j >= size(B.tab, 2) && return 0.0
    f = t - (j - 1)
    @inbounds (1 - f) * B.tab[iℓ, j] + f * B.tab[iℓ, j+1]
end

"""
    Θ_ℓ_of_k(S, iℓ, B, η0)

Project the source onto multipole ℓ:  Θ_ℓ(k) = ∫ dη S_T(η) j_ℓ(k(η₀-η)).
"""
function Θ_ℓ_of_k(S::SourceFunction, iℓ::Int, B::BesselTable, η0)
    acc = 0.0
    η, ST = S.η, S.S_T
    @inbounds for i in 1:(length(η)-1)
        dη = η[i+1] - η[i]
        f1 = ST[i] * B(iℓ, S.k * (η0 - η[i]))
        f2 = ST[i+1] * B(iℓ, S.k * (η0 - η[i+1]))
        acc += 0.5 * (f1 + f2) * dη
    end
    acc
end

"""
    ΘE_ℓ_of_k(S, ℓ, η0)

The E-mode transfer function. Polarization is generated only by the quadrupole Π
at scattering, and its projection carries an extra geometric factor:

    Θ^E_ℓ(k) = (3/4)·√[(ℓ+2)!/(ℓ-2)!] · ∫ dη g(η)Π(η) · j_ℓ(x)/x²,   x = k(η₀-η)

The 1/x² is why E-modes are sourced only by the *quadrupole* and vanish for a pure
monopole -- Thomson scattering polarizes only when the incident radiation is
anisotropic at ℓ = 2.
"""
function ΘE_ℓ_of_k(S::SourceFunction, iℓ::Int, ℓ::Int, B::BesselTable, η0)
    ℓ < 2 && return 0.0
    # √[(ℓ+2)!/(ℓ-2)!] = √[(ℓ+2)(ℓ+1)ℓ(ℓ-1)]
    pref = 0.75 * sqrt((ℓ + 2) * (ℓ + 1) * ℓ * (ℓ - 1))
    acc = 0.0
    η, SE = S.η, S.S_E
    @inbounds for i in 1:(length(η)-1)
        dη = η[i+1] - η[i]
        f = 0.0
        for (j, w) in ((i, 0.5), (i + 1, 0.5))
            x = S.k * (η0 - η[j])
            x < 1e-6 && continue
            f += w * SE[j] * B(iℓ, x) / x^2
        end
        acc += f * dη
    end
    pref * acc
end

"Unintegrated scalar sources required by the non-flat radial operators."
function _curved_scalar_source(c, bg, p; a_min=1e-5, nη=4000, Wlens=nothing)
    k = p.k
    L = _layout(p)
    xs = range(log(a_min), 0.0; length=nη)
    etas = [bg.η(x) for x in xs]
    S0 = zeros(nη); S1 = similar(S0); S2 = similar(S0); SP = similar(S0)
    SL = Wlens === nothing ? Float64[] : zeros(nη)
    gtheta_b = zeros(nη)
    for (i, x) in enumerate(xs)
        a = exp(x)
        u = p.sol(x)
        z = redshift(a)
        etau = exp(-optical_depth(p.rec, z))
        g = bg.κ̇(x) * etau
        phi, psi, phidot = _metric_at(c, p, u, L, a)
        Wlens === nothing || (SL[i] = -Wlens(etas[i]) * (phi + psi))
        P = (u[L.iγ+2] + u[L.iG] + u[L.iG+2]) / 8
        # Exact integration-by-parts form of the Newtonian-gauge temperature
        # source.  The direct form
        #   (e^-τ phi' + g Theta0, e^-τ k psi + g theta_b/k, g P)
        # is algebraically equivalent, but its high-k T0/T1 cancellation is
        # numerically ill-conditioned.  This is the stable form used in CLASS
        # perturbations.c; the boundary term vanishes because g theta_b is zero
        # both before visibility turns on and at the observer.
        S0[i] = g * (u[L.iγ] / 4 + phi) + 2etau * phidot
        S1[i] = etau * k * (psi - phi)
        S2[i] = g * P
        SP[i] = sqrt(6) * g * P
        gtheta_b[i] = g * u[4]
    end
    S0 .+= _ddη(gtheta_b, etas) ./ k^2
    (k=k, eta=collect(etas), S0=S0, S1=S1, S2=S2, SP=SP, SL=SL)
end

function _curved_scalar_radials(K, k, q, elllist, distance)
    R0 = zeros(length(elllist)); R1 = zeros(length(elllist))
    R2 = zeros(length(elllist)); RE = zeros(length(elllist))
    if distance <= 0 || sqrt(abs(K)) * distance < 1e-8
        R0[findall(==(0), elllist)] .= 1
        return R0, R1, R2, RE
    end
    rootK = sqrt(abs(K))
    sgnK = K > 0 ? 1 : -1
    chi = rootK * distance
    beta = q / rootK
    lmax = maximum(elllist)
    Phi, dPhi, d2Phi = _hyperspherical_bessels(sgnK, beta, lmax, chi)
    sinchi = _sinK(sgnK, chi)
    cscgen = rootK / (k * sinchi)
    rk = rootK / k
    s2 = _curvature_streaming(K, k, 2)
    for (i, ell) in enumerate(elllist)
        if sgnK == 1 && ell >= round(Int, beta)
            continue
        end
        phi = Phi[ell+1]
        R0[i] = phi
        R1[i] = rk * dPhi[ell+1]
        R2[i] = (3rk^2 * d2Phi[ell+1] + phi) / (2s2)
        if ell >= 2
            RE[i] = sqrt(3 / 8 * (ell + 2) * (ell + 1) * ell * (ell - 1)) /
                s2 * cscgen^2 * phi
        end
    end
    R0, R1, R2, RE
end

function _project_scalar_curved(S, elllist, K, k, q, eta0)
    # All three arrays are integration accumulators.  `similar(TT)` is
    # deliberately not used here: it returns uninitialised storage, so adding
    # the first quadrature panel would retain arbitrary heap contents.
    TT = zeros(length(elllist)); EE = zeros(length(elllist))
    LL = zeros(length(elllist))
    lens = hasproperty(S, :SL) && !isempty(S.SL)
    r00, r10, r20, re0 = _curved_scalar_radials(K, k, q, elllist, eta0 - S.eta[1])
    for i in 1:(length(S.eta)-1)
        r01, r11, r21, re1 = _curved_scalar_radials(K, k, q, elllist, eta0 - S.eta[i+1])
        deta = S.eta[i+1] - S.eta[i]
        @inbounds for j in eachindex(elllist)
            f0 = S.S0[i] * r00[j] + S.S1[i] * r10[j] + S.S2[i] * r20[j]
            f1 = S.S0[i+1] * r01[j] + S.S1[i+1] * r11[j] + S.S2[i+1] * r21[j]
            TT[j] += deta * (f0 + f1) / 2
            EE[j] += deta * (S.SP[i] * re0[j] + S.SP[i+1] * re1[j]) / 2
            # The lensing potential is a scalar on the sky: its radial
            # function is the plain hyperspherical Φ^ν_ℓ, same as the monopole.
            lens && (LL[j] += deta * (S.SL[i] * r00[j] + S.SL[i+1] * r01[j]) / 2)
        end
        r00, r10, r20, re0 = r01, r11, r21, re1
    end
    TT, EE, LL
end

function _scalar_mode_grid(K, nk, kmax)
    nk >= 12 || throw(ArgumentError("scalar spectra require nk >= 12"))
    if K < 0
        qmax2 = kmax^2 + K
        qmax2 > (5e-5)^2 || throw(ArgumentError("kmax does not reach the open scalar continuum"))
        qmax = sqrt(qmax2)
        qsplit = min(0.02, qmax / 2)
        nlo = min(max(8, nk ÷ 10), nk - 2)
        nhi = nk - nlo
        qs = vcat(exp.(range(log(5e-5), log(qsplit); length=nlo)),
            range(qsplit, qmax; length=nhi + 1)[2:end])
        return sqrt.(qs.^2 .- K), qs, false
    else
        rootK = sqrt(K)
        numax = max(3, floor(Int, sqrt(kmax^2 / K + 1)))
        nus = collect(3:numax)
        qs = rootK .* nus
        return sqrt.(qs.^2 .- K), qs, true
    end
end

# --- C_ℓ ---------------------------------------------------------------------

"""
    CMBSpectra

Angular power spectra. `TT`, `EE`, `TE` are `Vector`s indexed by ℓ from `ℓ_min`,
in μK² (raw C_ℓ, not ℓ(ℓ+1)C_ℓ/2π -- use [`D_ℓ`](@ref) for that).
"""
struct CMBSpectra
    ℓ::Vector{Int}
    TT::Vector{Float64}
    EE::Vector{Float64}
    TE::Vector{Float64}
    φφ::Vector{Float64}   # lensing potential C_ℓ^φφ (dimensionless); empty unless lensing=true
    Tφ::Vector{Float64}   # ISW-lensing cross C_ℓ^Tφ in μK; empty unless lensing=true
end

CMBSpectra(ℓ, TT, EE, TE) = CMBSpectra(ℓ, TT, EE, TE, Float64[], Float64[])

"""
    D_ℓ(spec, which = :TT)

ℓ(ℓ+1)C_ℓ/2π in μK² -- the form the CMB is conventionally plotted in, chosen
because a scale-invariant Sachs-Wolfe plateau is flat in it.
"""
function D_ℓ(s::CMBSpectra, which::Symbol=:TT)
    C = which === :TT ? s.TT : which === :EE ? s.EE : which === :TE ? s.TE :
        which === :φφ ? s.φφ : which === :Tφ ? s.Tφ :
        throw(ArgumentError("which must be :TT, :EE, :TE, :φφ or :Tφ"))
    [ℓ * (ℓ + 1) * C[i] / (2π) for (i, ℓ) in enumerate(s.ℓ)]
end

"""
    cmb_spectra(c, rec; lmax = 1500, nk = 3000, kmax = nothing, nk_solve = nothing)

Compute C_ℓ^TT, C_ℓ^EE and C_ℓ^TE.

    C_ℓ = 4π ∫ dln k  Δ²_ℛ(k)  Θ_ℓ(k)²

Two grids matter and both are easy to under-resolve:

* **k**. Θ_ℓ(k) oscillates with period ≈ 2π/η₀ ≈ 4.4e-4 /Mpc, so the k-grid must
  resolve *that*, not the smooth envelope. Under-sampling does not blur the
  spectrum, it aliases -- producing spurious wiggles that look like physics. The
  grid is therefore linear at high k, not logarithmic.
* **k_max**. A multipole ℓ is sourced by k ≈ ℓ/η₀, so k_max must reach beyond
  ℓ_max/η₀ or the highest multipoles are simply missing power.

Cost is one Boltzmann solve per k, threaded.

`nk_solve` decouples the two grids. The Bessel-scale oscillation above lives in
the *projection* Θ_ℓ(k) = ∫ dη S(k,η) j_ℓ(k(η₀−η)); the line-of-sight source
S(k,η) itself varies in k only on the sound-horizon scale Δk ≈ 2π/r_s ≈
0.04 /Mpc. Passing `nk_solve = n` therefore solves the Boltzmann hierarchy on
only `n` k-points, interpolates the source tables linearly onto the full
`nk`-point grid, and projects densely — the exact integral, evaluated on
interpolated sources. `nk_solve = 600` carries ~80 solve points per source
oscillation and reproduces the all-k result to ≲0.1% while cutting the
hierarchy cost ~8×. The default `nothing` solves at every k; tighten by raising
`nk_solve` and confirming the spectra stop changing, same as any grid knob.
"""
function cmb_spectra(c::Cosmology, rec::RecombinationSolution;
    lmax=1500, nk=3000, kmax=nothing, lmax_γ=25, lmax_ν=32,
    ℓs=nothing, nη=4000, sw=1.0, isw=1.0, dop=1.0, pol=1.0,
    lensing=false, lensing_kernel::Symbol=:lastscatter, ic::Symbol=:adiabatic,
    nk_solve=nothing, mg::Union{Nothing,HorndeskiFunctions}=nothing)

    bg = BackgroundCache(c, rec)
    η0 = bg.η(0.0)
    kmax === nothing && (kmax = 2.2 * lmax / η0)
    K = spatial_curvature_K(c)
    Wlens = lensing ? _lensing_kernel(c, bg, rec; kernel=lensing_kernel) : nothing

    if !iszero(K)
        (sw, isw, dop, pol) == (1.0, 1.0, 1.0, 1.0) || throw(ArgumentError(
            "effect switches are defined for the integrated flat source; curved spectra require all switches = 1"))
        ks, qs, discrete = _scalar_mode_grid(K, nk, kmax)
        ell_list = ℓs === nothing ? _default_ℓs(lmax) : collect(Int.(ℓs))
        nells = length(ell_list)
        ThetaT = zeros(nells, length(ks)); ThetaE = similar(ThetaT)
        ThetaL = zeros(nells, length(ks))
        Threads.@threads for j in eachindex(ks)
            p = solve_perturbations(c, bg, ks[j]; lmax_γ, lmax_ν, ic)
            S = _curved_scalar_source(c, bg, p; nη, Wlens)
            t, e, lphi = _project_scalar_curved(S, ell_list, K, ks[j], qs[j], η0)
            ThetaT[:, j] .= t; ThetaE[:, j] .= e; ThetaL[:, j] .= lphi
        end
        temp2 = (c.Tcmb * 1e6)^2
        TT = zeros(nells); EE = similar(TT); TE = similar(TT)
        φφ = lensing ? zeros(nells) : Float64[]
        Tφ = lensing ? zeros(nells) : Float64[]
        for i in eachindex(ell_list)
            if discrete
                weight(j) = qs[j] * sqrt(K) / ks[j]^2
                TT[i] = 4pi * sum(j -> primordial_power(c, ks[j]) * ThetaT[i,j]^2 * weight(j), eachindex(ks)) * temp2
                EE[i] = 4pi * sum(j -> primordial_power(c, ks[j]) * ThetaE[i,j]^2 * weight(j), eachindex(ks)) * temp2
                TE[i] = 4pi * sum(j -> primordial_power(c, ks[j]) * ThetaT[i,j] * ThetaE[i,j] * weight(j), eachindex(ks)) * temp2
                if lensing
                    φφ[i] = 4pi * sum(j -> primordial_power(c, ks[j]) * ThetaL[i,j]^2 * weight(j), eachindex(ks))
                    Tφ[i] = 4pi * sum(j -> primordial_power(c, ks[j]) * ThetaT[i,j] * ThetaL[i,j] * weight(j), eachindex(ks)) * sqrt(temp2)
                end
            else
                TT[i] = 4pi * _trapz(j -> primordial_power(c, ks[j]) * ThetaT[i,j]^2 / ks[j], ks) * temp2
                EE[i] = 4pi * _trapz(j -> primordial_power(c, ks[j]) * ThetaE[i,j]^2 / ks[j], ks) * temp2
                TE[i] = 4pi * _trapz(j -> primordial_power(c, ks[j]) * ThetaT[i,j] * ThetaE[i,j] / ks[j], ks) * temp2
                if lensing
                    φφ[i] = 4pi * _trapz(j -> primordial_power(c, ks[j]) * ThetaL[i,j]^2 / ks[j], ks)
                    Tφ[i] = 4pi * _trapz(j -> primordial_power(c, ks[j]) * ThetaT[i,j] * ThetaL[i,j] / ks[j], ks) * sqrt(temp2)
                end
            end
        end
        return CMBSpectra(ell_list, TT, EE, TE, φφ, Tφ)
    end

    # Log-spaced at low k (where the spectrum is smooth and few modes matter),
    # linear at high k (where the acoustic oscillations must be resolved).
    #
    # The linear spacing must be a small fraction of 2π/η₀ ≈ 4.4e-4. It is
    # tempting to economise here and it does not fail loudly: at Δk = 4.8e-5 the
    # first acoustic peak lands at ℓ = 180 instead of 221 and is 27% too high,
    # which looks like a plausible spectrum rather than an artefact. Only at
    # Δk ≲ 1.2e-5 does the peak converge onto CAMB's.
    # For deliberately low-l runs kmax can lie below the usual 0.02/Mpc
    # transition.  Keeping a fixed split in that case makes the second segment
    # run backwards and invalidates the quadrature.
    k_split = min(0.02, kmax / 2)
    n_lo = min(max(80, nk ÷ 10), nk - 2)
    n_hi = nk - n_lo
    ks = vcat(exp.(range(log(5e-5), log(k_split); length=n_lo)),
        range(k_split, kmax; length=n_hi + 1)[2:end])

    # The acoustic peaks converge only for Δk ≲ 1.2e-5 in the linear segment
    # (see the resolution note above). An under-resolved run does not fail
    # loudly — it produces a plausible-looking spectrum with the first peak
    # shifted and misweighted — and the error does NOT cancel between two
    # cosmologies (different η₀ ⇒ different k-grids).
    Δk_hi = (kmax - k_split) / n_hi
    if Δk_hi > 1.5e-5
        nk_need = ceil(Int, (kmax - k_split) / 1.2e-5 + n_lo)
        @warn "cmb_spectra: k-grid under-resolved for converged acoustic peaks" Δk = Δk_hi Δk_converged = 1.2e-5 nk = nk nk_recommended = nk_need maxlog = 1
    end

    ℓ_list = ℓs === nothing ? _default_ℓs(lmax) : ℓs
    B = BesselTable(ℓ_list, kmax * η0 * 1.02)

    nℓ = length(ℓ_list)
    ΘT = zeros(nℓ, length(ks))
    ΘE = zeros(nℓ, length(ks))
    ΘL = zeros(nℓ, length(ks))

    if nk_solve === nothing
        Threads.@threads for j in eachindex(ks)
            k = ks[j]
            p = solve_perturbations(c, bg, k; lmax_γ, lmax_ν, ic, mg)
            S = source_function(c, bg, p; nη, sw, isw, dop, pol, Wlens)
            for (i, ℓ) in enumerate(ℓ_list)
                ΘT[i, j] = Θ_ℓ_of_k(S, i, B, η0)
                ΘE[i, j] = ΘE_ℓ_of_k(S, i, ℓ, B, η0)
                lensing && (ΘL[i, j] = Θφ_ℓ_of_k(S, i, B, η0))
            end
        end
    else
        # Source-interpolation fast path. The line-of-sight SOURCES S(k,η) are
        # smooth in k — at fixed η they oscillate with the sound-horizon period
        # Δk ≈ 2π/r_s ≈ 0.04 Mpc⁻¹ — while the transfer functions Θ_ℓ(k)
        # oscillate on the far finer Bessel scale Δk ≈ π/η₀. So the expensive
        # Boltzmann hierarchy needs solving only on a grid that resolves the
        # sources; the dense projection grid reads them by interpolation.
        #
        # The solve grid must be *linear* across the whole acoustic range: the
        # source oscillation has a fixed k-period, so a log segment reaching to
        # k_split = 0.02 starves the first-peak region (k ≈ 0.015) of points —
        # errors there sat at 0.7% while ℓ ≥ 400 converged to 0.03%. Log
        # spacing is right only super-horizon, where the sources vary on the
        # scale of k itself; hand over to linear at 1e-3, safely inside the
        # horizon at recombination. With nk_solve = 600 the linear segment
        # carries ~140 points per source oscillation.
        k_log = min(1e-3, k_split)
        n_log = max(40, nk_solve ÷ 15)
        ks_solve = vcat(exp.(range(log(5e-5), log(k_log); length=n_log)),
            range(k_log, kmax; length=nk_solve - n_log + 1)[2:end])
        nks = length(ks_solve)
        ST_s = Matrix{Float64}(undef, nks, nη)
        SE_s = Matrix{Float64}(undef, nks, nη)
        SL_s = lensing ? Matrix{Float64}(undef, nks, nη) : Matrix{Float64}(undef, 0, 0)
        ηgrid = Ref{Vector{Float64}}()
        Threads.@threads for j in 1:nks
            p = solve_perturbations(c, bg, ks_solve[j]; lmax_γ, lmax_ν, ic, mg)
            S = source_function(c, bg, p; nη, sw, isw, dop, pol, Wlens)
            ST_s[j, :] = S.S_T
            SE_s[j, :] = S.S_E
            lensing && (SL_s[j, :] = S.S_L)
            j == 1 && (ηgrid[] = S.η)
        end
        ηs_common = ηgrid[]
        Threads.@threads for j in eachindex(ks)
            k = ks[j]
            # bracketing weights on the solve grid
            jr = searchsortedfirst(ks_solve, k)
            jl = clamp(jr - 1, 1, nks - 1)
            jr = jl + 1
            w = (k - ks_solve[jl]) / (ks_solve[jr] - ks_solve[jl])
            w = clamp(w, 0.0, 1.0)
            ST = @views (1 - w) .* ST_s[jl, :] .+ w .* ST_s[jr, :]
            SE = @views (1 - w) .* SE_s[jl, :] .+ w .* SE_s[jr, :]
            SL = lensing ? (@views (1 - w) .* SL_s[jl, :] .+ w .* SL_s[jr, :]) :
                 Float64[]
            S = SourceFunction(k, ηs_common, ST, SE, SL)
            for (i, ℓ) in enumerate(ℓ_list)
                ΘT[i, j] = Θ_ℓ_of_k(S, i, B, η0)
                ΘE[i, j] = ΘE_ℓ_of_k(S, i, ℓ, B, η0)
                lensing && (ΘL[i, j] = Θφ_ℓ_of_k(S, i, B, η0))
            end
        end
    end

    Tcmb_μK = c.Tcmb * 1e6
    TT = zeros(nℓ)
    EE = zeros(nℓ)
    TE = zeros(nℓ)
    φφ = lensing ? zeros(nℓ) : Float64[]
    Tφ = lensing ? zeros(nℓ) : Float64[]
    for i in 1:nℓ
        fT(j) = primordial_power(c, ks[j]) * ΘT[i, j]^2 / ks[j]
        fE(j) = primordial_power(c, ks[j]) * ΘE[i, j]^2 / ks[j]
        fX(j) = primordial_power(c, ks[j]) * ΘT[i, j] * ΘE[i, j] / ks[j]
        TT[i] = 4π * _trapz(fT, ks) * Tcmb_μK^2
        EE[i] = 4π * _trapz(fE, ks) * Tcmb_μK^2
        TE[i] = 4π * _trapz(fX, ks) * Tcmb_μK^2
        if lensing
            fL(j) = primordial_power(c, ks[j]) * ΘL[i, j]^2 / ks[j]
            fTL(j) = primordial_power(c, ks[j]) * ΘT[i, j] * ΘL[i, j] / ks[j]
            φφ[i] = 4π * _trapz(fL, ks)
            Tφ[i] = 4π * _trapz(fTL, ks) * Tcmb_μK
        end
    end

    CMBSpectra(ℓ_list, TT, EE, TE, φφ, Tφ)
end

"Trapezoid over the k-grid, with f given as a function of index."
function _trapz(f, ks)
    acc = 0.0
    @inbounds for j in 1:(length(ks)-1)
        acc += 0.5 * (f(j) + f(j + 1)) * (ks[j+1] - ks[j])
    end
    acc
end

# --- Correlated adiabatic + isocurvature spectra ----------------------------

"""
    PrimordialMatrix(modes; amp, tilt, k_pivot = 0.05)

Primordial cross-power specification for a set of scalar `modes` (any of
`:adiabatic`, `:cdi`, `:bi`, `:nid`, `:niv`). `amp` and `tilt` are symmetric
`n×n` matrices over `modes`; the dimensionless cross power of modes `i,j` is

    𝒫_ij(k) = amp[i,j] · (k/k_pivot)^(tilt[i,j] − 1),

so the diagonal holds each mode's auto-amplitude and tilt (the adiabatic entry is
just `A_s`, `n_s`) and the off-diagonal holds the cross-correlations. The total
harmonic spectrum sums every ordered pair,

    C_ℓ^XY = 4π Σ_ij ∫ dk/k 𝒫_ij(k) Θ_ℓ^{X,i}(k) Θ_ℓ^{Y,j}(k).

`amp` should be positive semi-definite for a physical spectrum; this is not
enforced (an unphysical matrix simply yields C_ℓ that can go negative).
"""
struct PrimordialMatrix
    modes::Vector{Symbol}
    amp::Matrix{Float64}
    tilt::Matrix{Float64}
    k_pivot::Float64
end

function PrimordialMatrix(modes::AbstractVector{Symbol};
    amp::AbstractMatrix, tilt::AbstractMatrix, k_pivot=0.05)
    n = length(modes)
    size(amp) == (n, n) || throw(ArgumentError("amp must be $(n)×$(n)"))
    size(tilt) == (n, n) || throw(ArgumentError("tilt must be $(n)×$(n)"))
    PrimordialMatrix(collect(Symbol, modes), Matrix{Float64}(amp),
        Matrix{Float64}(tilt), Float64(k_pivot))
end

@inline _prim_cross(pm::PrimordialMatrix, a, b, k) =
    pm.amp[a, b] * (k / pm.k_pivot)^(pm.tilt[a, b] - 1)

"""
    cmb_spectra(c, rec, prim::PrimordialMatrix; kwargs...)

Correlated scalar CMB spectra for a mixture of adiabatic and isocurvature modes,
combined through the primordial matrix `prim` (see [`PrimordialMatrix`](@ref)).
Each mode is evolved with its own initial condition, projected to transfer
functions Θ_ℓ^{T,i}, Θ_ℓ^{E,i}, and summed over all mode pairs. Flat geometry
only. Returns a [`CMBSpectra`](@ref) (φφ/Tφ empty).

Passing `PrimordialMatrix([:adiabatic]; amp=[A_s;;], tilt=[n_s;;])` reproduces the
ordinary adiabatic `cmb_spectra`.
"""
function cmb_spectra(c::Cosmology, rec::RecombinationSolution, prim::PrimordialMatrix;
    lmax=2500, nk=3000, kmax=nothing, lmax_γ=25, lmax_ν=32, ℓs=nothing, nη=4000)

    bg = BackgroundCache(c, rec)
    η0 = bg.η(0.0)
    K = spatial_curvature_K(c)
    iszero(K) || throw(ArgumentError(
        "correlated isocurvature spectra are implemented for flat geometry only"))
    kmax === nothing && (kmax = 2.2 * lmax / η0)

    k_split = min(0.02, kmax / 2)
    n_lo = min(max(80, nk ÷ 10), nk - 2)
    n_hi = nk - n_lo
    ks = vcat(exp.(range(log(5e-5), log(k_split); length=n_lo)),
        range(k_split, kmax; length=n_hi + 1)[2:end])

    # The acoustic peaks converge only for Δk ≲ 1.2e-5 in the linear segment
    # (see the resolution note above). An under-resolved run does not fail
    # loudly — it produces a plausible-looking spectrum with the first peak
    # shifted and misweighted — and the error does NOT cancel between two
    # cosmologies (different η₀ ⇒ different k-grids).
    Δk_hi = (kmax - k_split) / n_hi
    if Δk_hi > 1.5e-5
        nk_need = ceil(Int, (kmax - k_split) / 1.2e-5 + n_lo)
        @warn "cmb_spectra: k-grid under-resolved for converged acoustic peaks" Δk = Δk_hi Δk_converged = 1.2e-5 nk = nk nk_recommended = nk_need maxlog = 1
    end

    ℓ_list = ℓs === nothing ? _default_ℓs(lmax) : ℓs
    B = BesselTable(ℓ_list, kmax * η0 * 1.02)
    nℓ = length(ℓ_list)
    nm = length(prim.modes)

    # Per-mode transfer functions on the shared (ℓ, k) grid.
    ΘT = [zeros(nℓ, length(ks)) for _ in 1:nm]
    ΘE = [zeros(nℓ, length(ks)) for _ in 1:nm]
    for (m, mode) in enumerate(prim.modes)
        Threads.@threads for j in eachindex(ks)
            p = solve_perturbations(c, bg, ks[j]; lmax_γ, lmax_ν, ic=mode)
            S = source_function(c, bg, p; nη)
            for (i, ℓ) in enumerate(ℓ_list)
                ΘT[m][i, j] = Θ_ℓ_of_k(S, i, B, η0)
                ΘE[m][i, j] = ΘE_ℓ_of_k(S, i, ℓ, B, η0)
            end
        end
    end

    Tcmb_μK = c.Tcmb * 1e6
    TT = zeros(nℓ)
    EE = zeros(nℓ)
    TE = zeros(nℓ)
    for i in 1:nℓ
        for a in 1:nm, b in 1:nm
            fT(j) = _prim_cross(prim, a, b, ks[j]) * ΘT[a][i, j] * ΘT[b][i, j] / ks[j]
            fE(j) = _prim_cross(prim, a, b, ks[j]) * ΘE[a][i, j] * ΘE[b][i, j] / ks[j]
            fX(j) = _prim_cross(prim, a, b, ks[j]) * ΘT[a][i, j] * ΘE[b][i, j] / ks[j]
            TT[i] += 4π * _trapz(fT, ks) * Tcmb_μK^2
            EE[i] += 4π * _trapz(fE, ks) * Tcmb_μK^2
            TE[i] += 4π * _trapz(fX, ks) * Tcmb_μK^2
        end
    end
    CMBSpectra(ℓ_list, TT, EE, TE)
end

"""
Multipoles to actually compute. C_ℓ is smooth in ℓ, so sampling it and
interpolating is standard practice and saves a large factor -- but the sampling
must stay dense through the acoustic peaks, where it is *not* slowly varying.
"""
function _default_ℓs(lmax)
    ℓs = Int[]
    ℓ = 2
    while ℓ <= lmax
        push!(ℓs, ℓ)
        step = ℓ < 50 ? 1 : ℓ < 200 ? 5 : ℓ < 1000 ? 10 : 25
        ℓ += step
    end
    ℓs
end
