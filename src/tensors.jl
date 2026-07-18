"""
Tensor perturbations and their CMB temperature/E/B-mode spectra.

The evolution equations use the total-angular-momentum hierarchy of Hu & White,
Phys. Rev. D 56, 596 (1997), in the same flat-space normalization used by CLASS
3.3.4.  The gravitational wave is a dynamical metric degree of freedom,

    h'' + 2 (a'/a) h' + k²h = 8πGa² πᵀ,

and the source πᵀ is computed from the evolving photon, massless-neutrino, and
momentum-resolved massive-neutrino distributions.  Thus particle anisotropic
stress feeds back into the metric rather than being appended after the solve.

No tight-coupling or massless-neutrino approximation is used.  The early photon
system is stiff and is integrated implicitly, as in the scalar solver.
"""

struct TensorLayout
    lmax_γ::Int
    lmax_ν::Int
    lmax_m::Int
    nq::Int
    nmν::Int
    iF::Int
    iG::Int
    iN::Int
    iM::Int
    ih::Int
    ihdot::Int
    n::Int
end

function TensorLayout(lmax_γ, lmax_ν; lmax_m=8, nq=0, nmν=0)
    lmax_γ >= 4 || throw(ArgumentError("tensor photon hierarchy needs lmax_γ ≥ 4"))
    lmax_ν >= 4 || throw(ArgumentError("tensor neutrino hierarchy needs lmax_ν ≥ 4"))
    lmax_m >= 4 || throw(ArgumentError("tensor massive-neutrino hierarchy needs lmax_m ≥ 4"))
    iF = 1
    iG = iF + lmax_γ + 1
    iN = iG + lmax_γ + 1
    iM = iN + lmax_ν + 1
    ih = iM + nmν * nq * (lmax_m + 1)
    ihdot = ih + 1
    TensorLayout(lmax_γ, lmax_ν, lmax_m, nq, nmν,
        iF, iG, iN, iM, ih, ihdot, ihdot)
end

@inline _itν(L::TensorLayout, s, iq, ℓ) =
    L.iM + ((s - 1) * L.nq + (iq - 1)) * (L.lmax_m + 1) + ℓ

struct TensorPerturbationSolution{C,R,S}
    cosmo::C
    rec::R
    k::Float64
    sol::S
    x_start::Float64
    lmax_γ::Int
    lmax_ν::Int
    lmax_m::Int
    nq::Int
    nmν::Int
end

_tensor_layout(p::TensorPerturbationSolution) = TensorLayout(
    p.lmax_γ, p.lmax_ν; lmax_m=p.lmax_m, nq=p.nq, nmν=p.nmν)

@inline function _tensor_P(u, L)
    F0 = u[L.iF]
    F2 = u[L.iF+2]
    F4 = u[L.iF+4]
    G0 = u[L.iG]
    G2 = u[L.iG+2]
    G4 = u[L.iG+4]
    -(F0 / 10 + F2 / 7 + 3F4 / 70 - 3G0 / 5 + 6G2 / 7 - 3G4 / 70) / sqrt(6)
end

@inline _tensor_stress_combo(u, i0) =
    u[i0] / 15 + 2u[i0+2] / 21 + u[i0+4] / 35

function _tensor_massive_source(u, L, G, mν, a, H0)
    source = 0.0
    for (s, m) in enumerate(mν)
        am = m.y1 * a
        integral = 0.0
        @inbounds for iq in 1:L.nq
            q = G.q[iq]
            ε = sqrt(q^2 + am^2)
            i0 = _itν(L, s, iq, 0)
            combo = u[i0] / 15 + 2u[i0+2] / 21 + u[i0+4] / 35
            integral += G.w[iq] * G.f0[iq] * q^4 / ε * combo
        end
        # m.Ω_rel/F_massless is the same phase-space normalization used by the
        # massive-neutrino background and scalar hierarchy.
        norm = m.Ω_rel / F_massless / a^4
        source += -4sqrt(6) * H0^2 * a^2 * norm * integral
    end
    source
end

function _tensor_rhs!(du, u, p, x)
    bg, k, K, L, G, mν = p
    mgf = length(p) >= 7 ? p[7] : nothing
    c = bg.cosmo
    a = exp(x)
    Hconf = getproperty(bg, Symbol("\u210b"))(x)
    η = bg.η(x)
    κ = bg.κ̇(x)
    H0 = Constants.H0_in_invMpc(c.h)
    hdot = u[L.ihdot]
    P = _tensor_P(u, L)
    s(ell) = _curvature_streaming(K, k, ell)
    cotend = _curvature_cot_closure(K, k, η)

    # Photon temperature hierarchy F_l^(2), including the exact Thomson source.
    @inbounds begin
        du[L.iF] = -k * u[L.iF+1] + sqrt(6) * hdot -
                    κ * (u[L.iF] + sqrt(6) * P)
        for ℓ in 1:(L.lmax_γ-1)
            du[L.iF+ℓ] = k / (2ℓ + 1) *
                (ℓ * s(ℓ) * u[L.iF+ℓ-1] -
                 (ℓ + 1) * s(ℓ + 1) * u[L.iF+ℓ+1]) - κ * u[L.iF+ℓ]
        end
        ℓ = L.lmax_γ
        du[L.iF+ℓ] = k * (s(ℓ) * u[L.iF+ℓ-1] - (ℓ + 1) * cotend * u[L.iF+ℓ]) -
                         κ * u[L.iF+ℓ]

        # Linear polarization hierarchy G_l^(2).
        du[L.iG] = -k * u[L.iG+1] - κ * (u[L.iG] - sqrt(6) * P)
        for ℓ in 1:(L.lmax_γ-1)
            du[L.iG+ℓ] = k / (2ℓ + 1) *
                (ℓ * s(ℓ) * u[L.iG+ℓ-1] -
                 (ℓ + 1) * s(ℓ + 1) * u[L.iG+ℓ+1]) - κ * u[L.iG+ℓ]
        end
        ℓ = L.lmax_γ
        du[L.iG+ℓ] = k * (s(ℓ) * u[L.iG+ℓ-1] - (ℓ + 1) * cotend * u[L.iG+ℓ]) -
                         κ * u[L.iG+ℓ]

        # Massless, collisionless neutrinos.
        du[L.iN] = -k * u[L.iN+1] + sqrt(6) * hdot
        for ℓ in 1:(L.lmax_ν-1)
            du[L.iN+ℓ] = k / (2ℓ + 1) *
                (ℓ * s(ℓ) * u[L.iN+ℓ-1] -
                 (ℓ + 1) * s(ℓ + 1) * u[L.iN+ℓ+1])
        end
        ℓ = L.lmax_ν
        du[L.iN+ℓ] = k * (s(ℓ) * u[L.iN+ℓ-1] - (ℓ + 1) * cotend * u[L.iN+ℓ])
    end

    # Exact momentum-dependent massive-neutrino tensor hierarchy.
    for (s, m) in enumerate(mν), iq in 1:L.nq
        q = G.q[iq]
        ε = sqrt(q^2 + (m.y1 * a)^2)
        qkε = k * q / ε
        i0 = _itν(L, s, iq, 0)
        dlnf = G.dlnf0[iq]
        du[i0] = -qkε * u[i0+1] - sqrt(6) * hdot * dlnf / 4
        @inbounds for ℓ in 1:(L.lmax_m-1)
            du[i0+ℓ] = qkε / (2ℓ + 1) *
                (ℓ * s(ℓ) * u[i0+ℓ-1] -
                 (ℓ + 1) * s(ℓ + 1) * u[i0+ℓ+1])
        end
        ℓ = L.lmax_m
        du[i0+ℓ] = qkε * s(ℓ) * u[i0+ℓ-1] - (ℓ + 1) * k * cotend * u[i0+ℓ]
    end

    ργ = Ω_γ(c) / a^4
    ρν = _Ω_or_zero(c, MasslessNeutrinos) / a^4
    gw_source = -4sqrt(6) * H0^2 * a^2 *
        (ργ * _tensor_stress_combo(u, L.iF) + ρν * _tensor_stress_combo(u, L.iN))
    L.nmν > 0 && (gw_source += _tensor_massive_source(u, L, G, mν, a, H0))

    du[L.ih] = hdot
    # Horndeski (alpha_T = 0): the run of the Planck mass adds Hubble friction
    # and rescales the matter source, h'' + (2+α_M)ℋh' + k²h = source/M*²
    # (B&S 1404.3713 tensor sector; the GW speed is untouched by assumption).
    αM_t, M2_t = mgf === nothing ? (0.0, 1.0) : (mgf.α_M(x), mgf.M_star2(x))
    du[L.ihdot] = -(2 + αM_t) * Hconf * hdot - (k^2 + 2K) * u[L.ih] +
                  gw_source / M2_t

    # The integration variable is ln(a), not conformal time.
    du ./= Hconf
    nothing
end

function _tensor_initial_conditions(c, bg, k, K, L, G, mν, x0)
    η = bg.η(x0)
    Ωγ = Ω_γ(c)
    Ωfs = _Ω_or_zero(c, MasslessNeutrinos) + sum(m.Ω_rel for m in mν; init=0.0)
    Ωr = Ωγ + Ωfs

    # CLASS/Hu-White normalization: h0=1/sqrt(6), with the tensor primordial
    # spectrum A_t=r A_s applied only when forming C_l.  The O((kη)^2) term is
    # the regular solution including free-streaming anisotropic stress.
    h0 = 1 / sqrt(6)
    if !iszero(K)
        q2 = k^2 + 3K
        q2 > 0 || throw(ArgumentError(
            "tensor super-curvature modes q²=k²+3K≤0 require separate initial data"))
        h0 *= sqrt(k^2 * (k^2 - K) / ((k^2 + 3K) * (k^2 + 2K)))
        K < 0 && (h0 *= sqrt(tanh(pi * sqrt(q2 / -K) / 2)))
    end
    h2 = -h0 * (k^2 + 2K) * η^2 / (6 + (8 / 5) * Ωfs / Ωr)
    u = zeros(L.n)
    u[L.ih] = h0 + h2
    u[L.ihdot] = 2h2 / η
    u[L.iN] = sqrt(6) * h2

    for (s, _) in enumerate(mν), iq in 1:L.nq
        u[_itν(L, s, iq, 0)] = sqrt(6) * h2 * (-G.dlnf0[iq] / 4)
    end
    u
end

"""
    solve_tensor_perturbations(c, bg, k; ...)

Evolve one primordial gravitational-wave mode and the particle hierarchies that
source and are sourced by it.  For nonzero curvature this uses the exact curved
streaming coefficients, curved hierarchy closure, `(k²+2K)h` propagation term,
and curvature-normalised regular initial mode.  No flat radial projection is
performed by this function.
"""
function solve_tensor_perturbations(c::Cosmology, bg::BackgroundCache, k::Real;
    lmax_γ=12, lmax_ν=16, lmax_m=12, nq=20,
    a_end=1.0, reltol=1e-8, abstol=1e-11, mg::Union{Nothing,HorndeskiFunctions}=nothing)

    K = spatial_curvature_K(c)
    mν = Tuple(get_all_species(c, MassiveNeutrinos))
    G = MassiveNuGrid(nq)
    L = TensorLayout(lmax_γ, lmax_ν; lmax_m, nq, nmν=length(mν))
    k = float(k)

    x_horizon = log(1e-9)
    for x in range(log(1e-9), log(1e-2); length=600)
        k * bg.η(x) > 1e-2 && break
        x_horizon = x
    end
    x0 = clamp(min(x_horizon, log(1e-5)), log(1e-10), log(a_end))
    u0 = _tensor_initial_conditions(c, bg, k, K, L, G, mν, x0)
    prob = ODEProblem(_tensor_rhs!, u0, (x0, log(a_end)), (bg, k, K, L, G, mν, mg))
    sol = solve(prob, KenCarp4(); reltol, abstol)
    TensorPerturbationSolution(c, bg.rec, k, sol, x0,
        lmax_γ, lmax_ν, lmax_m, nq, length(mν))
end

solve_tensor_perturbations(c::Cosmology, rec::RecombinationSolution, k::Real; kws...) =
    solve_tensor_perturbations(c, BackgroundCache(c, rec), k; kws...)

tensor_h(p::TensorPerturbationSolution, a) = p.sol(log(a))[_tensor_layout(p).ih]
tensor_hdot(p::TensorPerturbationSolution, a) = p.sol(log(a))[_tensor_layout(p).ihdot]

struct TensorSourceFunction
    k::Float64
    η::Vector{Float64}
    S_T::Vector{Float64}
    S_P::Vector{Float64}
end

function tensor_source_function(bg::BackgroundCache, p::TensorPerturbationSolution; nη=3000)
    L = _tensor_layout(p)
    xs = range(p.x_start, 0.0; length=nη)
    ηs = [bg.η(x) for x in xs]
    ST = zeros(nη)
    SP = zeros(nη)
    for (i, x) in enumerate(xs)
        u = p.sol(x)
        z = redshift(exp(x))
        eτ = exp(-optical_depth(p.rec, z))
        g = bg.κ̇(x) * eτ
        P = _tensor_P(u, L)
        ST[i] = -u[L.ihdot] * eτ + g * P
        # Historical CMBFAST/CAMB sign convention, also used by CLASS.
        SP[i] = sqrt(6) * g * P
    end
    TensorSourceFunction(p.k, collect(ηs), ST, SP)
end

struct TensorBesselTable
    ℓs::Vector{Int}
    dx::Float64
    j::Matrix{Float64}
    jp::Matrix{Float64}
    jpp::Matrix{Float64}
end

function TensorBesselTable(ℓs::Vector{Int}, xmax; dx=0.08)
    nx = ceil(Int, xmax / dx) + 2
    jt = zeros(length(ℓs), nx)
    jpt = similar(jt)
    jppt = similar(jt)
    Threads.@threads for i in eachindex(ℓs)
        ℓ = ℓs[i]
        for n in 2:nx
            x = (n - 1) * dx
            j = sphericalbesselj(ℓ, x)
            jp = ℓ * j / x - sphericalbesselj(ℓ + 1, x)
            jt[i, n] = j
            jpt[i, n] = jp
            jppt[i, n] = -2jp / x - (1 - ℓ * (ℓ + 1) / x^2) * j
        end
    end
    TensorBesselTable(ℓs, dx, jt, jpt, jppt)
end

@inline function _tensor_bessel(B::TensorBesselTable, iℓ, x)
    if x <= 0
        return (0.0, 0.0, 0.0)
    end
    t = x / B.dx
    n = floor(Int, t) + 1
    n >= size(B.j, 2) && return (0.0, 0.0, 0.0)
    f = t - (n - 1)
    @inbounds begin
        j = (1 - f) * B.j[iℓ, n] + f * B.j[iℓ, n+1]
        jp = (1 - f) * B.jp[iℓ, n] + f * B.jp[iℓ, n+1]
        jpp = (1 - f) * B.jpp[iℓ, n] + f * B.jpp[iℓ, n+1]
    end
    (j, jp, jpp)
end

@inline function _tensor_radials(B, iℓ, ℓ, x)
    if x == 0
        # j_2(x)=x²/15+O(x⁴).  The tensor temperature and E operators
        # therefore both tend to 1/5; all ℓ>2 functions vanish at the origin.
        return ℓ == 2 ? (0.2, 0.2, 0.0) : (0.0, 0.0, 0.0)
    end
    if x < 4B.dx
        # Linear interpolation cannot preserve j_l∝x^l through the first cell;
        # evaluate the exact special function and analytic recurrences there.
        j = sphericalbesselj(ℓ, x)
        jp = ℓ * j / x - sphericalbesselj(ℓ + 1, x)
        jpp = -2jp / x - (1 - ℓ * (ℓ + 1) / x^2) * j
    else
        j, jp, jpp = _tensor_bessel(B, iℓ, x)
    end
    T = sqrt(3 / 8 * (ℓ + 2) * (ℓ + 1) * ℓ * (ℓ - 1)) * j / x^2
    E = (jpp + 4jp / x - (1 - 2 / x^2) * j) / 4
    Bm = (jp + 2j / x) / 2
    (T, E, Bm)
end

function _project_tensor(S, iℓ, ℓ, B, η0)
    tT = tE = tB = 0.0
    for i in 1:(length(S.η)-1)
        dη = S.η[i+1] - S.η[i]
        aT = aE = aB = 0.0
        for j in (i, i + 1)
            x = S.k * (η0 - S.η[j])
            rT, rE, rB = _tensor_radials(B, iℓ, ℓ, x)
            aT += S.S_T[j] * rT
            aE += S.S_P[j] * rE
            aB += S.S_P[j] * rB
        end
        tT += dη * aT / 2
        tE += dη * aE / 2
        tB += dη * aB / 2
    end
    (tT, tE, tB)
end

"Exact curved-FRW tensor temperature/E/B radial functions for all requested l."
function _curved_tensor_radials(K, k, q, elllist, distance)
    T = zeros(length(elllist)); E = zeros(length(elllist)); Bm = zeros(length(elllist))
    if distance <= 0 || sqrt(abs(K)) * distance < 1e-8
        for (i, ell) in enumerate(elllist)
            if ell == 2
                T[i], E[i], Bm[i] = 0.2, 0.2, 0.0
            end
        end
        return T, E, Bm
    end

    rootK = sqrt(abs(K))
    sgnK = K > 0 ? 1 : -1
    chi = rootK * distance
    beta = q / rootK
    lmax = maximum(elllist)
    Phi, dPhi, d2Phi = _hyperspherical_bessels(sgnK, beta, lmax, chi)
    sinchi = _sinK(sgnK, chi)
    coschi = _cosK(sgnK, chi)
    cscgen = rootK / (k * sinchi)
    cotgen = cscgen * coschi
    rk = rootK / k
    k2 = k^2
    sminus = sqrt(1 - K / k2)
    splus = sqrt(1 + 2K / k2)
    sb = sqrt(1 + 3K / k2)

    for (i, ell) in enumerate(elllist)
        # A closed mode with beta<=ell has no such angular multipole.
        if sgnK == 1 && ell >= round(Int, beta)
            continue
        end
        phi = Phi[ell+1]; dphi = dPhi[ell+1]; d2phi = d2Phi[ell+1]
        T[i] = sqrt(3 / 8 * (ell + 2) * (ell + 1) * ell * (ell - 1)) /
            (splus * sminus) * cscgen^2 * phi
        E[i] = (rk^2 * d2phi + 4cotgen * rk * dphi -
            (1 + 4K / k2 - 2cotgen^2) * phi) / (4splus * sminus)
        Bm[i] = sb / (2sminus * splus) * (rk * dphi + 2cotgen * phi)
    end
    T, E, Bm
end

"Project one curved tensor source without any flat-Bessel substitution."
function _project_tensor_curved(S, elllist, K, k, q, eta0)
    # These are quadrature accumulators, not output buffers.  They must all
    # start at the additive identity; `similar` would leave E/B undefined.
    tT = zeros(length(elllist)); tE = zeros(length(elllist));
    tB = zeros(length(elllist))
    rT0, rE0, rB0 = _curved_tensor_radials(K, k, q, elllist, eta0 - S.η[1])
    for i in 1:(length(S.η)-1)
        rT1, rE1, rB1 = _curved_tensor_radials(K, k, q, elllist, eta0 - S.η[i+1])
        dη = S.η[i+1] - S.η[i]
        @inbounds for j in eachindex(elllist)
            tT[j] += dη * (S.S_T[i] * rT0[j] + S.S_T[i+1] * rT1[j]) / 2
            tE[j] += dη * (S.S_P[i] * rE0[j] + S.S_P[i+1] * rE1[j]) / 2
            tB[j] += dη * (S.S_P[i] * rB0[j] + S.S_P[i+1] * rB1[j]) / 2
        end
        rT0, rE0, rB0 = rT1, rE1, rB1
    end
    tT, tE, tB
end

function _tensor_mode_grid(K, eta0, lmax, nk, kmax)
    nk >= 12 || throw(ArgumentError("tensor spectra require nk >= 12"))
    if iszero(K)
        ksplit = min(0.015, kmax / 2)
        nlo = min(max(8, nk ÷ 4), nk - 2)
        nhi = nk - nlo
        ks = vcat(exp.(range(log(2e-6), log(ksplit); length=nlo)),
            range(ksplit, kmax; length=nhi + 1)[2:end])
        return ks, copy(ks), false
    elseif K < 0
        # q is continuous for the sub-curvature spectrum; k²=q²-3K.
        qmax2 = kmax^2 + 3K
        qmax2 > (2e-6)^2 || throw(ArgumentError("kmax does not reach the open tensor continuum"))
        qmax = sqrt(qmax2)
        qsplit = min(0.015, qmax / 2)
        nlo = min(max(8, nk ÷ 4), nk - 2)
        nhi = nk - nlo
        qs = vcat(exp.(range(log(2e-6), log(qsplit); length=nlo)),
            range(qsplit, qmax; length=nhi + 1)[2:end])
        return sqrt.(qs.^2 .- 3K), qs, false
    else
        # On S^3, beta=q/sqrt(K)=3,4,... is discrete for tensors.  Retain every
        # eigenmode: subsampling this list would be a geometric approximation.
        rootK = sqrt(K)
        numax = max(3, floor(Int, sqrt(kmax^2 / K + 3)))
        nus = collect(3:numax)
        qs = rootK .* nus
        return sqrt.(qs.^2 .- 3K), qs, true
    end
end

struct TensorCMBSpectra
    ℓ::Vector{Int}
    TT::Vector{Float64}
    EE::Vector{Float64}
    BB::Vector{Float64}
    TE::Vector{Float64}
end

function D_ℓ(s::TensorCMBSpectra, which::Symbol=:BB)
    C = which === :TT ? s.TT : which === :EE ? s.EE : which === :BB ? s.BB :
        which === :TE ? s.TE : throw(ArgumentError("tensor spectrum must be :TT, :EE, :BB, or :TE"))
    [ℓ * (ℓ + 1) * C[i] / (2π) for (i, ℓ) in enumerate(s.ℓ)]
end

"""
    tensor_cmb_spectra(c, rec; r=0.01, n_t=-r/8, ...)

Primordial tensor CMB spectra.  `r` is the tensor-to-scalar ratio at `k_pivot`;
the default tilt is the single-field consistency relation but can be supplied
independently.  The returned spectra are raw C_l in μK² and include primordial
tensor TT, EE, BB, and TE (not lensing-generated scalar BB).
"""
function tensor_cmb_spectra(c::Cosmology, rec::RecombinationSolution;
    r=0.01, n_t=-r / 8, lmax=300, nk=nothing, kmax=nothing, ℓs=nothing,
    lmax_γ=12, lmax_ν=16, lmax_m=12, nq=20, nη=3000)

    r >= 0 || throw(ArgumentError("tensor-to-scalar ratio r must be nonnegative"))
    bg = BackgroundCache(c, rec)
    η0 = bg.η(0.0)
    # The tensor E/B radial operators do not decay at large x (E → -j/2,
    # B → j'/2), so EE and BB keep accumulating from k far above lmax/η0 --
    # truncating at the usual ~2·lmax/η0 leaves them up to an order of
    # magnitude low at low lmax while TT/TE (whose radial falls as j/x²) look
    # converged. The grid must reach the tensor source cutoff (~0.06/Mpc,
    # where h' at recombination has died) and resolve the transfer oscillation
    # period 2π/η0 with ~40 points per cycle, or EE/BB alias by several
    # percent. Validated against CLASS (flat, open) and CAMB (closed) at
    # ℓ ≤ 20 to ≤1%.
    kmax === nothing && (kmax = max(2.0 * lmax / η0, 0.06))
    nk === nothing && (nk = ceil(Int, 40 * kmax * η0 / (2π)))
    K = spatial_curvature_K(c)
    ks, qs, discrete = _tensor_mode_grid(K, η0, lmax, nk, kmax)
    ℓlist = ℓs === nothing ? _default_ℓs(lmax) : collect(Int.(ℓs))
    B = iszero(K) ? TensorBesselTable(ℓlist, 1.02 * kmax * η0) : nothing
    ΔT = zeros(length(ℓlist), length(ks))
    ΔE = similar(ΔT)
    ΔB = similar(ΔT)

    Threads.@threads for ik in eachindex(ks)
        p = solve_tensor_perturbations(c, bg, ks[ik]; lmax_γ, lmax_ν, lmax_m, nq)
        S = tensor_source_function(bg, p; nη)
        if iszero(K)
            for (i, ℓ) in enumerate(ℓlist)
                ΔT[i, ik], ΔE[i, ik], ΔB[i, ik] = _project_tensor(S, i, ℓ, B, η0)
            end
        else
            tT, tE, tB = _project_tensor_curved(S, ℓlist, K, ks[ik], qs[ik], η0)
            ΔT[:, ik] .= tT; ΔE[:, ik] .= tE; ΔB[:, ik] .= tB
        end
    end

    Pt(k) = r * c.A_s * (k / c.k_pivot)^n_t
    TμK2 = (c.Tcmb * 1e6)^2
    TT = zeros(length(ℓlist)); EE = similar(TT); BB = similar(TT); TE = similar(TT)
    for i in eachindex(ℓlist)
        if discrete
            # Integral dk is replaced by the exact S^3 mode sum with
            # dk=(q/k)dq and dq=sqrt(K).
            weight(j) = qs[j] * sqrt(K) / ks[j]^2
            TT[i] = 4π * sum(j -> Pt(ks[j]) * ΔT[i, j]^2 * weight(j), eachindex(ks)) * TμK2
            EE[i] = 4π * sum(j -> Pt(ks[j]) * ΔE[i, j]^2 * weight(j), eachindex(ks)) * TμK2
            BB[i] = 4π * sum(j -> Pt(ks[j]) * ΔB[i, j]^2 * weight(j), eachindex(ks)) * TμK2
            TE[i] = 4π * sum(j -> Pt(ks[j]) * ΔT[i, j] * ΔE[i, j] * weight(j), eachindex(ks)) * TμK2
        else
            TT[i] = 4π * _trapz(j -> Pt(ks[j]) * ΔT[i, j]^2 / ks[j], ks) * TμK2
            EE[i] = 4π * _trapz(j -> Pt(ks[j]) * ΔE[i, j]^2 / ks[j], ks) * TμK2
            BB[i] = 4π * _trapz(j -> Pt(ks[j]) * ΔB[i, j]^2 / ks[j], ks) * TμK2
            TE[i] = 4π * _trapz(j -> Pt(ks[j]) * ΔT[i, j] * ΔE[i, j] / ks[j], ks) * TμK2
        end
    end
    TensorCMBSpectra(ℓlist, TT, EE, BB, TE)
end
