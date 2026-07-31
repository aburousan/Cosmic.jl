"""
Finite-temperature (plasma) QED corrections to the n <-> p weak rates.

Brown & Sawyer (2001), as laid out in Pitrou et al. Phys. Rep. 754 (2018) eq. 107 and
§5.15-5.16. This is the last correction of the four, and the one that closes the gap to
CLASS/PRIMAT.

Everything above the Born rate so far has been a *zero-temperature* correction: the
Coulomb attraction of the outgoing electron, the Sirlin radiative factor, the recoil of
a finite-mass nucleon. But the decay does not happen in vacuum. It happens inside a
plasma at T ~ 1 MeV, dense with real photons and real electrons and positrons, and that
changes the rate three ways:

  * **Real photon exchange with the bath.** The electron can absorb or emit a *thermal*
    photon during the transition. Weighted by the Bose occupation of the bath, this is
    the `_L_photon` term.
  * **Bremsstrahlung, differenced.** Emission of a photon into (or out of) an already
    occupied mode. Only the *difference* from the vacuum bremsstrahlung already counted
    in the Sirlin factor belongs here, which is what `_L_brems` computes -- adding the
    whole thing would double-count.
  * **The thermal mass shift and e-e scattering.** The electron's dispersion relation is
    not the vacuum one inside a plasma: it picks up a thermal mass δm² ~ e²T²/6. That
    shifts the available phase space (`_L_mass`, a 1D integral) and adds electron-electron
    and electron-positron scattering (`_L_ee`, a 2D one).

None of it survives to T = 0, which is why -- unlike the other three -- this correction
does *not* enter the neutron-lifetime normalisation F_n. It changes the rates without
changing the decay they are normalised against, so it does not cancel, and it is worth
about -0.2% on Y_p.

**On tables.** PRIMAT and PRyMordial both ship this as a precomputed file, and PRyMordial
evaluates it with `vegas` (adaptive Monte Carlo). Neither is a statement that the physics
is uncomputable -- those tables were produced by doing exactly these integrals, once. They
are a cache, because at fixed neutrino-decoupling history and chemical potential
the result is a function of temperature alone. So the standard ξ=0 table is
computed here from scratch and stored in `data/bbn/thermal_qed.dat`; deleting it
recomputes the production quadrature, which takes minutes on a multi-core server.
Nonzero ξ values build a separate in-memory interpolation and never reuse the ξ=0
table.

**Numerics.** The true-photon and differential-bremsstrahlung pieces contain a severe
infrared cancellation. They are therefore combined pointwise and integrated on a
logarithmic photon-energy grid; convergence of their small residual must be checked, not
inferred from convergence of either large term separately. The remaining threshold
structure gives a leading `1/n` Gauss-grid error, so production tables use the
Richardson limit `2L(2n)-L(n)` (with final order 256). Below `10^8.2 K` the
electron/positron/photon bath is dilute and the well-behaved terms are negligible, while
the finite photon-momentum cutoff can leave a spurious bremsstrahlung residual. Following
current PRIMAT, both directions are set exactly to zero there.
"""

using QuadGK: gauss, quadgk

# --- shared kernels ----------------------------------------------------------

@inline _fd2(en, x, ξ=0.0) =
    (a = x * en - ξ; a > 300 ? 0.0 : a < -300 ? 1.0 : 1 / (exp(a) + 1))

"Derivative of the Fermi-Dirac occupation with respect to energy."
@inline function _dfd2(en, x)
    a = x * en
    abs(a) > 300 && return 0.0
    e = exp(a)
    -x * e / (1 + e)^2
end

"Bose occupation of the photon bath."
@inline _be(y) = (y <= 0 || y > 300) ? 0.0 : 1 / expm1(y)

"χ̃(E) = (E - s q)² · f_ν(E - s q). The neutrino side of the matrix element."
@inline function _chitilde(en, znu, s, ξ_ν=0.0)
    d = en - s * _q_np
    d^2 * _fd2(d, znu, s * ξ_ν)
end

"χ(E) = χ̃(E) · [1 - f_e(E)], the full Born integrand kernel."
@inline _chi_born(en, x, znu, s, ξ_ν=0.0) =
    _chitilde(en, znu, s, ξ_ν) * _fd2(-en, x)

# The two photon-emission form factors, Brown & Sawyer.
@inline function _A_ff(E, k)
    p = sqrt(max(E^2 - 1, 0.0))
    p <= 0 && return 0.0
    (2E^2 + k^2) * log((E + p) / (E - p)) - 4p * E
end
@inline function _B_ff(E)
    p = sqrt(max(E^2 - 1, 0.0))
    p <= 0 && return 0.0
    2E * log((E + p) / (E - p)) - 4p
end

"Coulomb factor, applied only when the produced charged lepton is the electron."
@inline _fstat(s, sgnE, b) = (s * sgnE > 0) ? _fermi_coulomb(b) : 1.0

# --- (1) real thermal-photon exchange ---------------------------------------

function _int_photon(E, k, x, znu, s, ξ_ν=0.0)
    E <= 1 && return 0.0
    p = sqrt(E^2 - 1)
    b = p / E
    A = _A_ff(E, k)
    B = _B_ff(E)
    fp = _fd2(-E, x) * _fstat(s, 1, b)
    fm = _fd2(E, x) * _fstat(s, -1, b)

    t1 = A * (fp * (_chitilde(E - k, znu, s, ξ_ν) + _chitilde(E + k, znu, s, ξ_ν) - 2 * _chitilde(E, znu, s, ξ_ν)) +
              fm * (_chitilde(-E + k, znu, s, ξ_ν) + _chitilde(-E - k, znu, s, ξ_ν) - 2 * _chitilde(-E, znu, s, ξ_ν)))
    t2 = k * B * (fp * (_chitilde(E - k, znu, s, ξ_ν) - _chitilde(E + k, znu, s, ξ_ν)) +
                  fm * (_chitilde(-E + k, znu, s, ξ_ν) - _chitilde(-E - k, znu, s, ξ_ν)))

    _ALPHA_FS / (2π) * (_be(x * k) / k) * (t1 - t2)
end

# --- (2) bremsstrahlung, differenced against the vacuum piece ----------------

function _int_brems(E, k, x, znu, s, ξ_ν=0.0)
    E <= 1 && return 0.0
    p = sqrt(E^2 - 1)
    b = p / E
    q = _q_np
    A = _A_ff(E, k)
    B = _B_ff(E)
    Fp = A + k * B
    Fm = A - k * B

    r1 = Fp * _chitilde(E + k, znu, s, ξ_ν)
    if abs(k) < abs(E - s * q)
        r1 -= Fp * _fd2(E - s * q, znu, s * ξ_ν) * (abs(E - s * q) - k)^2
    end
    r1 *= _fd2(-E, x) * _fstat(s, 1, b)

    r2 = Fm * _chitilde(-E + k, znu, s, ξ_ν)
    if abs(k) < abs(E + s * q)
        r2 -= Fp * _fd2(-E - s * q, znu, s * ξ_ν) * (abs(E + s * q) - k)^2
    end
    r2 *= _fd2(E, x) * _fstat(s, -1, b)

    _ALPHA_FS / (2π * k) * (r1 + r2)
end

# --- (3) thermal mass shift (1D) --------------------------------------------

function _int_mass(E, x, znu, s, ξ_ν=0.0)
    E <= 1 && return 0.0
    p = sqrt(E^2 - 1)
    -(_ALPHA_FS * E / (2π * p)) * (2π^2) / (3 * x^2) *
    (_chi_born(E, x, znu, s, ξ_ν) + _chi_born(-E, x, znu, s, ξ_ν))
end

# --- (4) electron-electron / electron-positron scattering (2D) ---------------

function _int_ee(e1, e2, x, znu, s, ξ_ν=0.0)
    (e1 <= 1 || e2 <= 1) && return 0.0
    p1 = sqrt(e1^2 - 1)
    p2 = sqrt(e2^2 - 1)
    (p1 <= 0 || p2 <= 0) && return 0.0
    abs(p1 - p2) < 1e-12 && return 0.0            # integrable log singularity at p1 = p2

    L = log((e1 * e2 + p1 * p2 + 1) / (e1 * e2 - p1 * p2 + 1))
    G = log(((p1 + p2) / (p1 - p2))^2)
    F2 = _fd2(e2, x)
    dF2 = _dfd2(e2, x)

    br = -(1 / 4) * G^2 * (dF2 * (p2 / p1) * (e1^2 / e2) * (e1 + e2) +
                           F2 * e1^2 / (p1 * p2) * (e2 + e1 / e2^2)) +
         G * (dF2 * (p2^2 * e1 / e2 * (1 / p1^2 + 2) - e1^2 * p2 / p1 * L) +
              F2 * (e1 / (p1^2 * e2^2) * (e2^2 + 2p1^2 + 1) - (e1^2 + e2^2) / (e1 + e2) -
                    (e1^2 * e2) / (p1 * p2) * L)) -
         F2 * (4 * e1 * p2 / p1 + 2 * e2 * L)

    _ALPHA_FS / (2π) *
    (_chi_born(e1, x, znu, s, ξ_ν) + _chi_born(-e1, x, znu, s, ξ_ν)) * br
end

# --- assembly ----------------------------------------------------------------

"2D tensor-product Gauss-Legendre over [a1,b1] x [a2,b2]."
function _gl2(f, a1, b1, a2, b2; n=96)
    nodes, weights = _gauss_rule(n)
    h1 = (b1 - a1) / 2
    h2 = (b2 - a2) / 2
    m1 = (a1 + b1) / 2
    m2 = (a2 + b2) / 2
    acc = 0.0
    @inbounds for i in eachindex(nodes), j in eachindex(nodes)
        u = m1 + h1 * nodes[i]
        v = m2 + h2 * nodes[j]
        acc += weights[i] * weights[j] * f(u, v)
    end
    h1 * h2 * acc
end

"""
    thermal_weak_correction(T_γ, T_ν; direction, n)

The finite-temperature QED correction to the n→p (`direction = +1`) or p→n (`-1`)
rate, in the same units as the Born integral Λ (so it is added to Λ, then divided by
τ_n·F_n).

Both corrections are switched off below T = 10^8.2 K, following current PRIMAT: the
thermal plasma contribution is negligible there and the finite photon-momentum cutoff
otherwise leaves a spurious infrared residual.
"""
function _thermal_grid_correction(x, znu, s, n, ξ_ν)
    M = max(10.0, 20 / x)

    # The photon-exchange and bremsstrahlung terms are individually large and of
    # opposite sign -- at T ~ 1 MeV they are +8.2 and -5.7, and their sum is +0.54. Six
    # percent of either. Integrating them separately and subtracting is catastrophic:
    # a 1% quadrature error in each becomes a 16% error in the answer, and the result
    # creeps toward the truth so slowly that n = 768 is still 13% off. This is the IR
    # cancellation between real and virtual emission, and it is exact -- so cancel it
    # *pointwise, inside the integrand*, and let the quadrature see only the small
    # residual it actually has to resolve.
    #
    # The k-grid is logarithmic for the same reason: both integrands carry a 1/k (and
    # the photon term a Bose 1/k²) whose divergence is killed by a numerator vanishing
    # to matching order, so all the structure sits at small k. A uniform grid on
    # [0.001, 40] puts almost no points there. Substituting u = ln k (dk = k du) both
    # concentrates the nodes and makes the transformed integrand vanish at the endpoint.
    lnk_lo, lnk_hi = log(1e-4), log(M)
    Lphbr = _gl2((E, u) -> begin
            k = exp(u)
            k * (_int_photon(E, k, x, znu, s, ξ_ν) +
                 _int_brems(E, k, x, znu, s, ξ_ν))
        end, 1.001, M, lnk_lo, lnk_hi; n)

    # e-e scattering, in the (e1+e2, e1-e2) variables Brown & Sawyer use. The Jacobian
    # is 1/2. The domain splits at e1 = e2, where the integrand has an integrable
    # logarithmic singularity, so the two halves are done separately and the singular
    # line is never a quadrature node.
    W = max(10.0, 15 / x)
    lo, hi = 2.002, 2 + W
    g(sp, dm) = 0.5 * _int_ee(
        (sp + dm) / 2, (sp - dm) / 2, x, znu, s, ξ_ν)
    Lee = _gl2(g, lo, hi, -W, -0.001; n) + _gl2(g, lo, hi, 0.001, W; n)

    Lphbr + Lee
end

_thermal_mass_correction(x, znu, s, ξ_ν) =
    quadgk(E -> _int_mass(E, x, znu, s, ξ_ν),
        1.0, max(25.0, 150 / x); rtol=1e-8)[1]

function thermal_weak_correction(T_γ, T_ν; direction=1, n=256, ξ_ν=0.0)
    # T in Kelvin, for the dilute-plasma/infrared cutoff.
    T_K = T_γ / 8.617333262e-11
    T_K < 10^8.2 && return 0.0
    s = direction
    x = _ME_MEV / T_γ
    znu = _ME_MEV / T_ν
    _thermal_grid_correction(x, znu, s, n, ξ_ν) +
        _thermal_mass_correction(x, znu, s, ξ_ν)
end

"""
Richardson-extrapolated thermal correction.

The kinematic step functions leave a leading error proportional to the inverse
tensor-grid order. Evaluating at `n/2` and `n` therefore removes that term as
`2L(n)-L(n/2)`. This is the production value used in caches; the unextrapolated
function remains public for explicit convergence studies.
"""
function _thermal_weak_converged(T_γ, T_ν; direction=1, n=256, ξ_ν=0.0)
    iseven(n) || throw(ArgumentError("thermal Richardson order n must be even"))
    n >= 4 || throw(ArgumentError("thermal Richardson order n must be at least 4"))
    T_γ / 8.617333262e-11 < 10^8.2 && return 0.0
    x = _ME_MEV / T_γ
    znu = _ME_MEV / T_ν
    coarse = _thermal_grid_correction(x, znu, direction, n ÷ 2, ξ_ν)
    fine = _thermal_grid_correction(x, znu, direction, n, ξ_ν)
    _thermal_mass_correction(x, znu, direction, ξ_ν) + 2fine - coarse
end

# --- tabulation and cache ----------------------------------------------------
#
# For fixed ξ and the solved Standard-Model decoupling trajectory, the correction
# depends on temperature rather than late-time cosmological parameters. It is
# therefore computed once on a temperature grid. The shipped disk cache is only
# the standard ξ=0 case; arbitrary ξ values use distinct in-memory tables.
#
# The difference here is that the cache is *ours*. `data/bbn/thermal_qed.dat` is generated
# by the integrals above, and deleting it simply causes them to run again. Nothing in it is
# taken on faith from another code.

const _THERMAL_CACHE = joinpath(_BBN_DIR, "thermal_qed.dat")

"Compute the thermal correction on a temperature grid and write the cache."
function generate_thermal_cache(; n_T=120, T_lo=4e-4, T_hi=9.0, n=256)
    Ts = exp.(range(log(T_lo), log(T_hi); length=n_T))       # MeV
    np = zeros(n_T)
    pn = zeros(n_T)
    Threads.@threads for i in eachindex(Ts)
        Tγ = Ts[i]
        Tν, _ = nu_temperatures(Tγ)
        np[i] = _thermal_weak_converged(Tγ, Tν; direction=1, n)
        pn[i] = _thermal_weak_converged(Tγ, Tν; direction=-1, n)
    end
    open(_THERMAL_CACHE, "w") do io
        println(io, "# Finite-temperature QED corrections to the n<->p weak rates.")
        println(io, "# Brown & Sawyer / Pitrou et al. Phys. Rep. 754 eq. 107, computed by")
        println(io, "# Cosmic.generate_thermal_cache(); Richardson 2L($n)-L($(n ÷ 2)).")
        println(io, "# T_gamma [MeV]        Lambda_th(n->p)        Lambda_th(p->n)")
        for i in eachindex(Ts)
            println(io, "$(Ts[i])  $(np[i])  $(pn[i])")
        end
    end
    Ts, np, pn
end

function _thermal_interpolants_for_xi(ξ_ν; n_T=120, T_lo=4e-4, T_hi=9.0, n=256)
    Ts = exp.(range(log(T_lo), log(T_hi); length=n_T))
    np = zeros(n_T)
    pn = zeros(n_T)
    Threads.@threads for i in eachindex(Ts)
        Tγ = Ts[i]
        Tν, _ = nu_temperatures(Tγ)
        np[i] = _thermal_weak_converged(Tγ, Tν; direction=1, n, ξ_ν)
        pn[i] = _thermal_weak_converged(Tγ, Tν; direction=-1, n, ξ_ν)
    end
    lnT = range(log(first(Ts)), log(last(Ts)); length=length(Ts))
    (cubic_spline_interpolation(lnT, np; extrapolation_bc=Line()),
        cubic_spline_interpolation(lnT, pn; extrapolation_bc=Line()))
end

function _thermal_interpolants()
    if !isfile(_THERMAL_CACHE)
        @info "Computing finite-temperature QED corrections to the weak rates (once; cached afterwards)"
        generate_thermal_cache()
    end
    d = readdlm(_THERMAL_CACHE, comments=true, comment_char='#')
    Ts = Float64.(d[:, 1])
    lnT = range(log(first(Ts)), log(last(Ts)); length=length(Ts))
    (cubic_spline_interpolation(lnT, Float64.(d[:, 2]); extrapolation_bc=Line()),
        cubic_spline_interpolation(lnT, Float64.(d[:, 3]); extrapolation_bc=Line()))
end

const _TH_NP = Ref{Any}(nothing)
const _TH_PN = Ref{Any}(nothing)
const _TH_XI = Dict{Float64,Tuple{Any,Any}}()
const _TH_XI_LOCK = ReentrantLock()

"Interpolated thermal correction Λ_th(T_γ), in the same units as the Born Λ."
function _thermal_lookup(T_γ, s, ξ_ν=0.0)
    T_γ / 8.617333262e-11 < 10^8.2 && return 0.0
    if iszero(ξ_ν)
        if _TH_NP[] === nothing
            a, b = _thermal_interpolants()
            _TH_NP[] = a
            _TH_PN[] = b
        end
        itp = s > 0 ? _TH_NP[] : _TH_PN[]
        return itp(log(T_γ))
    end
    key = Float64(ξ_ν)
    pair = lock(_TH_XI_LOCK) do
        get!(_TH_XI, key) do
            _thermal_interpolants_for_xi(key)
        end
    end
    itp = s > 0 ? pair[1] : pair[2]
    itp(log(T_γ))
end
