"""
Second-order Einstein--Boltzmann angular machinery.

This file contains pieces that are exact and independently testable before the
full second-order solver is exposed: the K=0 Fourier-triangle formulation,
rotations of first-order scalar multipoles, the complete (ell,m) hierarchy
layout, angular couplings, quadratic source groups, the linear-in-second-order
radiation operator, and flat line-of-sight projection kernels. The production
driver is added only after every Einstein, Liouville, collision, matter and
recombination source has a reference gate. This prevents a partial hierarchy
or a flat kernel from being mistaken for a general physical calculation.

Conventions and equations are those of Pitrou (arXiv:0809.3245) and Pettinari's
thesis (arXiv:1405.2280), chapters 5--6 and appendices A--B.  SONG is used only
as an independent numerical reference at matched flat, adiabatic settings.
"""

using SpecialFunctions: loggamma, sphericalbesselj
using OrdinaryDiffEq: ODEProblem, solve, Rodas5P

# --- Fourier triangle and rotations -----------------------------------------

"""
    SecondOrderTriangle(k1, k2, k)

The convolution geometry `kvec = k1vec + k2vec`, with `kvec` along z,
`phi1 = 0` and `phi2 = pi`.  The constructor enforces the exact triangle
condition; degenerate triangles are retained because squeezed/collinear limits
are physical and are required convergence gates.
"""
struct SecondOrderTriangle{T<:AbstractFloat}
    k1::T
    k2::T
    k::T
    cos1::T
    sin1::T
    cos2::T
    sin2::T
end

function SecondOrderTriangle(k1::Real, k2::Real, k::Real)
    T = promote_type(typeof(float(k1)), typeof(float(k2)), typeof(float(k)))
    a, b, c = T(k1), T(k2), T(k)
    (a > 0 && b > 0 && c > 0) || throw(ArgumentError("k1, k2 and k must be positive"))
    tol = 64eps(T) * max(a, b, c)
    c <= a + b + tol || throw(DomainError((a, b, c), "k > k1+k2 violates triangle closure"))
    c + tol >= abs(a - b) || throw(DomainError((a, b, c), "k < |k1-k2| violates triangle closure"))
    c1 = clamp((c*c + a*a - b*b) / (2c*a), -one(T), one(T))
    c2 = clamp((c*c - a*a + b*b) / (2c*b), -one(T), one(T))
    s1 = sqrt(max(zero(T), one(T) - c1*c1))
    s2 = sqrt(max(zero(T), one(T) - c2*c2))
    # This is a stronger closure check than the cosine formula alone near a
    # degenerate triangle: transverse momenta must cancel.
    abs(a*s1 - b*s2) <= 256eps(T) * max(a, b, c) ||
        throw(ArgumentError("loss of precision while constructing Fourier triangle"))
    SecondOrderTriangle{T}(a, b, c, c1, s1, c2, s2)
end

"Spherical component xi_[m].k_i in the Pettinari/Pitrou helicity convention."
function _k_spherical(t::SecondOrderTriangle, which::Int, m::Int)
    m in -1:1 || throw(ArgumentError("a vector has only m=-1,0,1 components"))
    if which == 1
        kval, co, si, azsign = t.k1, t.cos1, t.sin1, one(t.k1)
    elseif which == 2
        kval, co, si, azsign = t.k2, t.cos2, t.sin2, -one(t.k2) # phi2=pi
    else
        throw(ArgumentError("which must be 1 or 2"))
    end
    m == 0 ? kval * co : (m == -1 ? one(kval) : -one(kval)) * azsign * kval * si / sqrt(2)
end

@inline _logfact(n::Int) = loggamma(n + 1.0)

"Normalised scalar rotation sqrt(4pi/(2l+1)) Y_lm(theta,phi)."
function _scalar_rotation(l::Int, m::Int, costheta::Real, phi_is_pi::Bool=false)
    l >= 0 || return 0.0
    abs(m) <= l || return 0.0
    m < 0 && return (isodd(-m) ? -1.0 : 1.0) *
        _scalar_rotation(l, -m, costheta, phi_is_pi)
    x = clamp(float(costheta), -1.0, 1.0)
    # Associated Legendre P_m^m including the Condon--Shortley phase.
    pmm = 1.0
    if m > 0
        somx2 = sqrt(max(0.0, (1 - x) * (1 + x)))
        fact = 1.0
        for _ in 1:m
            pmm *= -fact * somx2
            fact += 2
        end
    end
    plm = if l == m
        pmm
    else
        pmmp1 = x * (2m + 1) * pmm
        if l == m + 1
            pmmp1
        else
            pprev, pcur = pmm, pmmp1
            for ll in (m + 2):l
                pnext = ((2ll - 1) * x * pcur - (ll + m - 1) * pprev) / (ll - m)
                pprev, pcur = pcur, pnext
            end
            pcur
        end
    end
    norm = exp(0.5 * (_logfact(l - m) - _logfact(l + m)))
    phase = phi_is_pi && isodd(m) ? -1.0 : 1.0
    phase * norm * plm
end

_rotate_scalar(t::SecondOrderTriangle, which::Int, l::Int, m::Int) =
    _scalar_rotation(l, m, which == 1 ? t.cos1 : t.cos2, which == 2)

"""
    _tensor_product(t, a, b, m)

Helicity-`m` component of the symmetric trace-free product
`k_a^{<i} k_b^{j>}` for `m=-2:2`. This is the tensor called
`k1_ten_k2[m]` in SONG and the tensor projection in Pettinari Appendix A.
"""
function _tensor_product(t::SecondOrderTriangle, a::Int, b::Int, m::Int)
    m in -2:2 || throw(ArgumentError("a rank-two STF tensor has only m=-2:2"))
    a in 1:2 && b in 1:2 || throw(ArgumentError("a and b must be 1 or 2"))
    am, a0, ap = (_k_spherical(t, a, mm) for mm in (-1, 0, 1))
    bm, b0, bp = (_k_spherical(t, b, mm) for mm in (-1, 0, 1))
    dotab = a == b ? (a == 1 ? t.k1^2 : t.k2^2) :
        (t.k^2 - t.k1^2 - t.k2^2) / 2
    m == -2 && return sqrt(2 / 3) * am * bm
    m == -1 && return (am * b0 + bm * a0) / sqrt(3)
    m == 0 && return a0 * b0 - dotab / 3
    m == 1 && return (ap * b0 + bp * a0) / sqrt(3)
    sqrt(2 / 3) * ap * bp
end

# --- Quadratic Einstein sources --------------------------------------------

"First-order Newtonian-gauge metric values at one wave number and time."
struct FirstOrderMetricSnapshot{T<:Number}
    phi::T
    psi::T
    phi_prime::T
    psi_prime::T
    phi_prime_prime::T
end

"""
    _einstein_quadratic_sources(t, q1, q2, Hc, Hc_prime;
        a=1, rho_dipole1=0, rho_dipole2=0)

Exact quadratic Newtonian-gauge Einstein kernels for one Fourier triangle.
The result uses Cosmic/SONG's convention `X=Xbar+X1+X2/2`, so the
symmetrised kernels in Pettinari thesis eqs. (3.96)--(3.100) receive their
required overall factor two.

`rho_dipole1/2` must be the complete sums `sum_s rho_s dipole_s` over every
active matter and radiation species. They enter the tetrad contribution to
the longitudinal constraint.
"""
function _einstein_quadratic_sources(t::SecondOrderTriangle,
    q1::FirstOrderMetricSnapshot, q2::FirstOrderMetricSnapshot,
    Hc::Real, Hc_prime::Real; a::Real=1,
    rho_dipole1::Number=0, rho_dipole2::Number=0)
    k1, k2, k = t.k1, t.k2, t.k
    k1sq, k2sq, ksq = k1^2, k2^2, k^2
    k1k2 = (ksq - k1sq - k2sq) / 2
    k10, k20 = _k_spherical(t, 1, 0), _k_spherical(t, 2, 0)
    p1, p2 = q1.phi, q2.phi
    s1, s2 = q1.psi, q2.psi
    pd1, pd2 = q1.phi_prime, q2.phi_prime
    sd1, sd2 = q1.psi_prime, q2.psi_prime
    pdd1, pdd2 = q1.phi_prime_prime, q2.phi_prime_prime

    QTT = -12s1*s2*Hc^2 + p1*p2*(3k1k2 + 4k1sq + 4k2sq) +
          6pd2*(p1-s1)*Hc + 6pd1*(p2-s2)*Hc - 3pd1*pd2
    QST = 2s1*s2*Hc*(k10+k20) - 2s1*p2*Hc*k10 - 2s2*p1*Hc*k20 -
          2p1*pd2*(k10+2k20) - 2p2*pd1*(2k10+k20) +
          s1*pd2*k10 + s2*pd1*k20 +
          (a^2/2) * ((k10/k1)*rho_dipole1*(s2+p2) +
                     (k20/k2)*rho_dipole2*(s1+p1))
    QTR = -12s1*s2*(Hc^2+2Hc_prime) + s1*s2*(ksq+k1sq+k2sq) +
          p1*p2*(3k1k2+4k1sq+4k2sq) + p1*s2*(k1k2-2k2sq) +
          p2*s1*(k1k2-2k1sq) + 6(pdd2+2Hc*pd2)*(p1-s1) +
          6(pdd1+2Hc*pd1)*(p2-s2) - 3sd2*(4Hc*s1+pd1) -
          3sd1*(4Hc*s2+pd2) + 3pd1*pd2

    QT = promote_type(typeof(QTT), typeof(QST), typeof(QTR))
    QSS = Dict{Int,QT}()
    QSS_prime = Dict{Int,QT}()
    for m in -2:2
        t12 = _tensor_product(t, 1, 2, m)
        t11 = _tensor_product(t, 1, 1, m)
        t22 = _tensor_product(t, 2, 2, m)
        QSS[m] = p1*p2*(-3t12-2t11-2t22) +
                 s1*s2*(-t12-t11-t22) + s1*p2*(t12+t11) +
                 s2*p1*(t12+t22)
        QSS_prime[m] = (pd1*p2+p1*pd2)*(-3t12-2t11-2t22) +
                       (sd1*s2+s1*sd2)*(-t12-t11-t22) +
                       (sd1*p2+s1*pd2)*(t12+t11) +
                       (sd2*p1+s2*pd1)*(t12+t22)
    end

    QTT *= 2
    QST *= 2
    QTR *= 2
    for m in -2:2
        QSS[m] *= 2
        QSS_prime[m] *= 2
    end
    iszero(Hc) && throw(DomainError(Hc, "Poisson form of phi' is singular at Hc=0"))
    scalar = (
        phi_prime_poisson=-QTT/(6Hc),
        phi_prime_longitudinal=QST/(2k),
        psi=-3QSS[0]/(2ksq),
        psi_prime=-3QSS_prime[0]/(2ksq)-3Hc*QSS[0]/ksq,
        phi_prime_prime=-(QTR+QTT)/6,
    )
    vector = Dict(m => sqrt(3)*QSS[m]/k for m in (-1, 1))
    tensor = Dict(m => -QSS[m] for m in (-2, 2))
    (QTT=QTT, QST=QST, QTR=QTR, QSS=QSS, QSS_prime=QSS_prime,
     scalar=scalar, omega_prime=vector, gamma_prime_prime=tensor)
end

# --- Angular state layout ----------------------------------------------------

"All (l,m) pairs with `lmin <= l <= lmax`, including both helicities."
function _lm_pairs(lmin::Int, lmax::Int)
    lmax < lmin && return Tuple{Int,Int}[]
    [(l, m) for l in lmin:lmax for m in -l:l]
end

"""
    SecondOrderLayout(lmax_I, lmax_E, lmax_N)

Dense index map for photon intensity, photon E/B polarization and massless
neutrino brightness.  Both signs of m are stored: this costs memory, but makes
reality, parity and exchange symmetries explicit gates rather than assumptions.
"""
struct SecondOrderLayout
    lmax_I::Int
    lmax_E::Int
    lmax_N::Int
    I::Dict{Tuple{Int,Int},Int}
    E::Dict{Tuple{Int,Int},Int}
    B::Dict{Tuple{Int,Int},Int}
    N::Dict{Tuple{Int,Int},Int}
    n::Int
end

function SecondOrderLayout(lmax_I::Int=12, lmax_E::Int=lmax_I, lmax_N::Int=12)
    lmax_I >= 2 || throw(ArgumentError("lmax_I must be at least 2"))
    lmax_E >= 2 || throw(ArgumentError("lmax_E must be at least 2"))
    lmax_N >= 2 || throw(ArgumentError("lmax_N must be at least 2"))
    n = 0
    make(ps) = Dict(p => (n += 1) for p in ps)
    I = make(_lm_pairs(0, lmax_I))
    E = make(_lm_pairs(2, lmax_E))
    B = make(_lm_pairs(2, lmax_E))
    N = make(_lm_pairs(0, lmax_N))
    SecondOrderLayout(lmax_I, lmax_E, lmax_N, I, E, B, N, n)
end

@inline _so_get(y, d, l, m) = get(d, (l, m), 0) == 0 ? zero(eltype(y)) : y[d[(l, m)]]

# --- C/D/R/K angular couplings ----------------------------------------------

function _Cplus(l::Int, m1::Int, m::Int)
    abs(m1 - m) <= 1 || return 0.0
    if m1 == m
        return sqrt(max(0, (l + 1)^2 - m^2)) / (2l + 3)
    end
    s = m1 - m # +1 or -1
    -sqrt(max(0, (l + 1 + s*m) * (l + 2 + s*m))) / (sqrt(2) * (2l + 3))
end

function _Cminus(l::Int, m1::Int, m::Int)
    l == 0 && return 0.0
    abs(m1 - m) <= 1 || return 0.0
    if m1 == m
        return sqrt(max(0, l^2 - m^2)) / (2l - 1)
    end
    s = m1 - m
    sqrt(max(0, (l - 1 - s*m) * (l - s*m))) / (sqrt(2) * (2l - 1))
end

_Dplus(l, m1, m) = l < 2 ? 0.0 : sqrt((l - 1) * (l + 3)) / (l + 1) * _Cplus(l, m1, m)
_Dminus(l, m1, m) = l < 2 ? 0.0 : sqrt(max(0, l*l - 4)) / l * _Cminus(l, m1, m)

function _Dzero(l::Int, m1::Int, m::Int)
    l < 2 && return 0.0
    if m1 == m
        return -2m / (l * (l + 1))
    elseif abs(m1 - m) == 1
        s = m1 - m
        return -s * sqrt(max(0, 2 * (l + 1 + s*m) * (l - s*m))) / (l * (l + 1))
    end
    0.0
end

_Rplus(l, m1, m) = -(l + 2) * _Cplus(l, m1, m)
_Rminus(l, m1, m) = (l - 1) * _Cminus(l, m1, m)
_Kplus(l, m1, m) = -(l + 2) * _Dplus(l, m1, m)
_Kminus(l, m1, m) = (l - 1) * _Dminus(l, m1, m)
_Kzero(l, m1, m) = -_Dzero(l, m1, m)

@inline _raw(v::AbstractVector, l::Int) = l < 0 || l >= length(v) ? zero(eltype(v)) : v[l + 1]

"First-order scalar photon/polarization/massless-relic multipoles in their k-aligned frame."
struct FirstOrderRadiationSnapshot{VI<:AbstractVector,VE<:AbstractVector,VN<:AbstractVector}
    I::VI
    E::VE
    N::VN
end

"""
    _scalar_E_from_G(G)

Recover the first-order scalar helicity multipoles `E_l` from Cosmic's
already-evolved polarization hierarchy `G_l`; no second polarization hierarchy
is integrated.  This is Tram--Lesgourgues (arXiv:1305.3261), eq. (2.40b),
multiplied by four so that `E` has the same brightness normalization as `I`:

    E_l = (2l+1)/sqrt((l-1)l(l+1)(l+2)) *
          [-l(l-1)G_l + sum_{j=0}^{l-2} i^(l-j)(1+(-1)^(l+j))(2j+1)G_j].

Only equal-parity terms survive, so the returned multipoles are real whenever
`G` is real.  This correspondence is used only for the flat second-order
Fourier solver.  It is not an exact curved-space shortcut; curved second-order
polarization requires the unreduced spin hierarchy and curved product harmonics.
"""
function _scalar_E_from_G(G::AbstractVector)
    length(G) >= 3 || throw(ArgumentError("scalar polarization needs G_0 through G_2"))
    T = promote_type(eltype(G), Float64)
    E = zeros(T, length(G))
    for l in 2:(length(G)-1)
        bracket = -l * (l - 1) * G[l+1]
        for j in (l % 2):2:(l-2)
            phase = isodd((l-j) ÷ 2) ? -one(T) : one(T)
            bracket += 2 * phase * (2j + 1) * G[j+1]
        end
        E[l+1] = (2l+1) * bracket / sqrt((l-1)*l*(l+1)*(l+2))
    end
    E
end

"Build K=0 quadratic-kernel radiation input from Cosmic's existing scalar solve."
function _flat_first_order_radiation_snapshot(p::PerturbationSolution, a::Real)
    iszero(spatial_curvature_K(p.cosmo)) || throw(ArgumentError(
        "curved second-order polarization requires the full curved spin hierarchy"))
    L = _layout(p)
    u = p.sol(log(float(a)))
    I = [(2l+1)*u[L.iγ+l] for l in 0:L.lmax_γ]
    G = [u[L.iG+l] for l in 0:L.lmax_γ]
    N = [(2l+1)*u[L.iν+l] for l in 0:L.lmax_ν]
    FirstOrderRadiationSnapshot(I, _scalar_E_from_G(G), N)
end

@inline function _rotated_raw(v::AbstractVector, t::SecondOrderTriangle,
    which::Int, l::Int, m::Int)
    _raw(v, l) * _rotate_scalar(t, which, l, m)
end

"Summed harmonic product used by the quadratic Liouville/collision operators."
function _angular_product(t::SecondOrderTriangle, l::Int, m::Int,
    vector_leg::Int, multipole_leg::Int, shift::Int, coupling::Function)
    ll = l + shift
    ll < 0 && return 0.0
    sum(_rotate_scalar(t, vector_leg, 1, m2) *
        _rotate_scalar(t, multipole_leg, ll, m - m2) *
        coupling(l, m - m2, m) for m2 in -1:1)
end

"""
    _quadratic_radiation_liouville(t, q1, q2, r1, r2, L)

Quadratic Newtonian-gauge Liouville sources for photon I/E/B and massless
relic brightness. This is the full time-delay, redshift and lensing block of
Pettinari thesis eqs. (4.143), (4.145), (4.146) / Pitrou Appendix A. It returns
the right-hand-side convention used by `_second_order_radiation_linear!`.
"""
function _quadratic_radiation_liouville(t::SecondOrderTriangle,
    q1::FirstOrderMetricSnapshot, q2::FirstOrderMetricSnapshot,
    r1::FirstOrderRadiationSnapshot, r2::FirstOrderRadiationSnapshot,
    L::SecondOrderLayout)
    T = promote_type(eltype(r1.I), eltype(r2.I), typeof(q1.phi), typeof(q2.phi))
    SI, SE, SB, SN = (Dict{Tuple{Int,Int},T}() for _ in 1:4)
    product(l, m, a, b, shift, f) = _angular_product(t, l, m, a, b, shift, f)
    for (l, m) in keys(L.I)
        s = l == 0 ? -8(q1.phi*q2.phi_prime + q2.phi*q1.phi_prime) : zero(T)
        if l == 1 && abs(m) <= 1
            s += 4(_k_spherical(t, 1, m)*q1.psi*(q2.psi-q2.phi) +
                   _k_spherical(t, 2, m)*q2.psi*(q1.psi-q1.phi))
        end
        # Time delay: vector and multipole from the same leg (22/11).
        s += _raw(r2.I, l+1)*product(l,m,2,2, 1,_Cplus)*t.k2*(q1.phi+q1.psi) -
             _raw(r2.I, l-1)*product(l,m,2,2,-1,_Cminus)*t.k2*(q1.phi+q1.psi) +
             _raw(r1.I, l+1)*product(l,m,1,1, 1,_Cplus)*t.k1*(q2.phi+q2.psi) -
             _raw(r1.I, l-1)*product(l,m,1,1,-1,_Cminus)*t.k1*(q2.phi+q2.psi)
        # Redshift.
        s += -4q1.phi_prime*_rotated_raw(r2.I,t,2,l,m) +
             4q1.psi*t.k1*(_raw(r2.I,l+1)*product(l,m,1,2, 1,_Cplus) -
                            _raw(r2.I,l-1)*product(l,m,1,2,-1,_Cminus)) -
             4q2.phi_prime*_rotated_raw(r1.I,t,1,l,m) +
             4q2.psi*t.k2*(_raw(r1.I,l+1)*product(l,m,2,1, 1,_Cplus) -
                            _raw(r1.I,l-1)*product(l,m,2,1,-1,_Cminus))
        # Lensing.
        s += t.k1*(q1.phi+q1.psi)*(_raw(r2.I,l+1)*product(l,m,1,2, 1,_Rplus) -
                                    _raw(r2.I,l-1)*product(l,m,1,2,-1,_Rminus)) +
             t.k2*(q2.phi+q2.psi)*(_raw(r1.I,l+1)*product(l,m,2,1, 1,_Rplus) -
                                    _raw(r1.I,l-1)*product(l,m,2,1,-1,_Rminus))
        SI[(l,m)] = -s # projected Liouville operator was on the left-hand side
    end

    for (l, m) in keys(L.E)
        s = -4q1.phi_prime*_rotated_raw(r2.E,t,2,l,m) +
            _raw(r2.E,l+1)*(t.k1*(q1.phi+q1.psi)*product(l,m,1,2, 1,_Kplus) +
                             t.k2*(q1.phi+q1.psi)*product(l,m,2,2, 1,_Dplus) +
                             4t.k1*q1.psi*product(l,m,1,2, 1,_Dplus)) -
            _raw(r2.E,l-1)*(t.k1*(q1.phi+q1.psi)*product(l,m,1,2,-1,_Kminus) +
                             t.k2*(q1.phi+q1.psi)*product(l,m,2,2,-1,_Dminus) +
                             4t.k1*q1.psi*product(l,m,1,2,-1,_Dminus)) -
            4q2.phi_prime*_rotated_raw(r1.E,t,1,l,m) +
            _raw(r1.E,l+1)*(t.k2*(q2.phi+q2.psi)*product(l,m,2,1, 1,_Kplus) +
                             t.k1*(q2.phi+q2.psi)*product(l,m,1,1, 1,_Dplus) +
                             4t.k2*q2.psi*product(l,m,2,1, 1,_Dplus)) -
            _raw(r1.E,l-1)*(t.k2*(q2.phi+q2.psi)*product(l,m,2,1,-1,_Kminus) +
                             t.k1*(q2.phi+q2.psi)*product(l,m,1,1,-1,_Dminus) +
                             4t.k2*q2.psi*product(l,m,2,1,-1,_Dminus))
        SE[(l,m)] = -s
        sb = -_raw(r2.E,l)*(t.k1*(q1.phi+q1.psi)*product(l,m,1,2,0,_Kzero) +
                             t.k2*(q1.phi+q1.psi)*product(l,m,2,2,0,_Dzero) +
                             4t.k1*q1.psi*product(l,m,1,2,0,_Dzero)) -
             _raw(r1.E,l)*(t.k2*(q2.phi+q2.psi)*product(l,m,2,1,0,_Kzero) +
                             t.k1*(q2.phi+q2.psi)*product(l,m,1,1,0,_Dzero) +
                             4t.k2*q2.psi*product(l,m,2,1,0,_Dzero))
        SB[(l,m)] = -sb
    end

    # Collisionless massless relics have the intensity Liouville hierarchy
    # with N replacing I and no Thomson term.
    for (l, m) in keys(L.N)
        s = l == 0 ? -8(q1.phi*q2.phi_prime + q2.phi*q1.phi_prime) : zero(T)
        if l == 1 && abs(m) <= 1
            s += 4(_k_spherical(t,1,m)*q1.psi*(q2.psi-q2.phi) +
                   _k_spherical(t,2,m)*q2.psi*(q1.psi-q1.phi))
        end
        s += -4q1.phi_prime*_rotated_raw(r2.N,t,2,l,m) +
             _raw(r2.N,l+1)*(t.k1*(q1.phi+q1.psi)*product(l,m,1,2,1,_Rplus) +
                              t.k2*(q1.phi+q1.psi)*product(l,m,2,2,1,_Cplus) +
                              4t.k1*q1.psi*product(l,m,1,2,1,_Cplus)) -
             _raw(r2.N,l-1)*(t.k1*(q1.phi+q1.psi)*product(l,m,1,2,-1,_Rminus) +
                              t.k2*(q1.phi+q1.psi)*product(l,m,2,2,-1,_Cminus) +
                              4t.k1*q1.psi*product(l,m,1,2,-1,_Cminus)) -
             4q2.phi_prime*_rotated_raw(r1.N,t,1,l,m) +
             _raw(r1.N,l+1)*(t.k2*(q2.phi+q2.psi)*product(l,m,2,1,1,_Rplus) +
                              t.k1*(q2.phi+q2.psi)*product(l,m,1,1,1,_Cplus) +
                              4t.k2*q2.psi*product(l,m,2,1,1,_Cplus)) -
             _raw(r1.N,l-1)*(t.k2*(q2.phi+q2.psi)*product(l,m,2,1,-1,_Rminus) +
                              t.k1*(q2.phi+q2.psi)*product(l,m,1,1,-1,_Cminus) +
                              4t.k2*q2.psi*product(l,m,2,1,-1,_Cminus))
        SN[(l,m)] = -s
    end
    (I=SI, E=SE, B=SB, N=SN)
end

"First-order baryon/electron data required by the exact Thomson kernel."
struct FirstOrderBaryonSnapshot{T<:Number}
    delta_b::T
    delta_e::T
    vpot_b::T
    vpot_gamma::T
end

"First-order density and irrotational velocity potential of a matter species."
struct FirstOrderMatterSnapshot{T<:Number}
    delta::T
    vpot::T
end


"""
    _quadratic_matter_liouville(t, q1, q2, d1, d2)

Quadratic Newtonian-gauge Liouville sources for a cold matter species in beta
moments: density, all three velocity helicities, and the velocity-generated
pressure and anisotropic stress. These are thesis sec. 4.6.2 and are used for
both baryons and cold dark matter before their distinct collision forces are
added. No perfect-fluid switch discards the generated rank-two moments.
"""
function _quadratic_matter_liouville(t::SecondOrderTriangle,
    q1::FirstOrderMetricSnapshot, q2::FirstOrderMetricSnapshot,
    d1::FirstOrderMatterSnapshot, d2::FirstOrderMatterSnapshot)
    k1sq, k2sq = t.k1^2, t.k2^2
    k1k2 = (t.k^2-k1sq-k2sq)/2
    mono = 2k1k2*(q1.psi-q1.phi)*d2.vpot + k2sq*(q1.phi+q1.psi)*d2.vpot +
           3q1.phi_prime*(d2.delta+2q2.phi) +
           2k1k2*(q2.psi-q2.phi)*d1.vpot + k1sq*(q2.phi+q2.psi)*d1.vpot +
           3q2.phi_prime*(d1.delta+2q1.phi)
    dipole = Dict{Int,typeof(mono)}()
    for m in -1:1
        dipole[m] = 3*_k_spherical(t,1,m)*(d2.delta*q1.psi -
                     4d1.vpot*q2.phi_prime + q2.phi*q1.psi) +
                     3*_k_spherical(t,2,m)*(d1.delta*q2.psi -
                     4d2.vpot*q1.phi_prime + q1.phi*q2.psi) -
                     3q1.psi*q2.psi*(_k_spherical(t,1,m)+_k_spherical(t,2,m))
    end
    velocity_metric = q1.psi*d2.vpot + q2.psi*d1.vpot
    pressure = 2k1k2*velocity_metric
    quadrupole = Dict(m => -15*_tensor_product(t,1,2,m)*velocity_metric for m in -2:2)
    (monopole=mono, dipole=dipole, pressure=pressure, quadrupole=quadrupole)
end

"""
    _quadratic_thomson(t, q1, q2, r1, r2, b1, b2, L;
        rho_gamma_over_rho_b)

Complete quadratic Thomson collision kernel for I/E/B, before multiplication
by `kappa_dot`. `delta_e` is an input from the perturbed recombination solve,
not an internal approximation: `delta_e=delta_b+delta_xe`. The returned baryon
monopole/dipole kernels enforce photon--baryon energy-momentum exchange exactly.
Equations are Pettinari thesis (4.151), (4.154), (4.157).
"""
function _quadratic_thomson(t::SecondOrderTriangle,
    q1::FirstOrderMetricSnapshot, q2::FirstOrderMetricSnapshot,
    r1::FirstOrderRadiationSnapshot, r2::FirstOrderRadiationSnapshot,
    b1::FirstOrderBaryonSnapshot, b2::FirstOrderBaryonSnapshot,
    L::SecondOrderLayout; rho_gamma_over_rho_b::Real)
    T = promote_type(eltype(r1.I), eltype(r2.I), typeof(b1.delta_b), typeof(b2.delta_b))
    CI, CE, CB = (Dict{Tuple{Int,Int},T}() for _ in 1:3)
    product(l,m,a,b,shift,f) = _angular_product(t,l,m,a,b,shift,f)
    v01, v02 = -t.k1*b1.vpot_b, -t.k2*b2.vpot_b
    u1(m) = -_k_spherical(t,1,m)*b1.vpot_b
    u2(m) = -_k_spherical(t,2,m)*b2.vpot_b
    k1k2 = (t.k^2-t.k1^2-t.k2^2)/2

    CI[(0,0)] = (4/3)*k1k2*(b1.vpot_b*(b2.vpot_gamma-b2.vpot_b) +
                              b2.vpot_b*(b1.vpot_gamma-b1.vpot_b))
    for (l,m) in keys(L.I)
        l == 0 && continue
        i1 = _rotated_raw(r1.I,t,1,l,m)
        i2 = _rotated_raw(r2.I,t,2,l,m)
        e1 = _rotated_raw(r1.E,t,1,l,m)
        e2 = _rotated_raw(r2.E,t,2,l,m)
        c1 = l == 1 ? 4u1(m)-i1 :
             l == 2 ? -i1 + (i1-sqrt(6)*e1)/10 : -i1
        c2 = l == 1 ? 4u2(m)-i2 :
             l == 2 ? -i2 + (i2-sqrt(6)*e2)/10 : -i2
        gain = zero(T)
        if l == 1
            gain = 3(product(l,m,1,2,-1,_Cminus)*v01*_raw(r2.I,0) +
                     product(l,m,2,1,-1,_Cminus)*v02*_raw(r1.I,0)) -
                   4(b1.delta_b*u2(m)+b2.delta_b*u1(m))
        elseif l == 2
            gain = product(l,m,1,2,-1,_Cminus)*v01*(7v02-_raw(r2.I,1)/2) +
                   product(l,m,2,1,-1,_Cminus)*v02*(7v01-_raw(r1.I,1)/2)
        elseif l == 3
            gain = (product(l,m,1,2,-1,_Cminus)*v01*(_raw(r2.I,2)-sqrt(6)*_raw(r2.E,2)) +
                    product(l,m,2,1,-1,_Cminus)*v02*(_raw(r1.I,2)-sqrt(6)*_raw(r1.E,2)))/2
        end
        loss = (q1.psi+b1.delta_e)*c2 + (q2.psi+b2.delta_e)*c1 +
               product(l,m,1,2,-1,_Cminus)*v01*_raw(r2.I,l-1) -
               product(l,m,1,2, 1,_Cplus)*v01*_raw(r2.I,l+1) +
               product(l,m,2,1,-1,_Cminus)*v02*_raw(r1.I,l-1) -
               product(l,m,2,1, 1,_Cplus)*v02*_raw(r1.I,l+1)
        CI[(l,m)] = gain + loss
    end

    for (l,m) in keys(L.E)
        e1 = _rotated_raw(r1.E,t,1,l,m)
        e2 = _rotated_raw(r2.E,t,2,l,m)
        c1 = l == 2 ? -e1-sqrt(6)*(_rotated_raw(r1.I,t,1,2,m)-sqrt(6)*e1)/10 : -e1
        c2 = l == 2 ? -e2-sqrt(6)*(_rotated_raw(r2.I,t,2,2,m)-sqrt(6)*e2)/10 : -e2
        gain = zero(T)
        if l == 2
            gain = sqrt(6)/2 *
                (product(l,m,1,2,-1,_Cminus)*v01*(_raw(r2.I,1)-2v02) +
                 product(l,m,2,1,-1,_Cminus)*v02*(_raw(r1.I,1)-2v01))
        elseif l == 3
            gain = -sqrt(6)/2 *
                (product(l,m,1,2,-1,_Dminus)*v01*(_raw(r2.I,2)-sqrt(6)*_raw(r2.E,2)) +
                 product(l,m,2,1,-1,_Dminus)*v02*(_raw(r1.I,2)-sqrt(6)*_raw(r1.E,2)))
        end
        loss = (q1.psi+b1.delta_e)*c2 + (q2.psi+b2.delta_e)*c1 +
               product(l,m,1,2,-1,_Dminus)*v01*_raw(r2.E,l-1) -
               product(l,m,1,2, 1,_Dplus)*v01*_raw(r2.E,l+1) +
               product(l,m,2,1,-1,_Dminus)*v02*_raw(r1.E,l-1) -
               product(l,m,2,1, 1,_Dplus)*v02*_raw(r1.E,l+1)
        CE[(l,m)] = gain + loss

        bgain = l == 2 ? -sqrt(6)/5 *
            (product(l,m,1,2,0,_Dzero)*v01*(_raw(r2.I,2)-sqrt(6)*_raw(r2.E,2)) +
             product(l,m,2,1,0,_Dzero)*v02*(_raw(r1.I,2)-sqrt(6)*_raw(r1.E,2))) : zero(T)
        bloss = product(l,m,1,2,0,_Dzero)*v01*_raw(r2.E,l) +
                product(l,m,2,1,0,_Dzero)*v02*_raw(r1.E,l)
        CB[(l,m)] = bgain + bloss
    end
    baryon_monopole = -rho_gamma_over_rho_b * CI[(0,0)]
    baryon_dipole = Dict(m => -rho_gamma_over_rho_b*CI[(1,m)] for m in -1:1)
    (I=CI, E=CE, B=CB, baryon_monopole=baryon_monopole,
     baryon_dipole=baryon_dipole)
end

"""
    _second_order_radiation_linear!(dy, y, L, k, kappa_dot; velocity=zeros(3))

Exact part of the second-order photon/neutrino hierarchy that is linear in the
second-order variables. `velocity[m+2]` is the electron/baryon helicity velocity
u_[m], m=-1:1. Quadratic Liouville/collision and metric terms are *not* silently
set here; the caller must add certified source vectors separately.
"""
function _second_order_radiation_linear!(dy, y, L::SecondOrderLayout, k, kappa_dot;
    velocity=zeros(eltype(y), 3))
    length(y) == L.n == length(dy) || throw(DimensionMismatch("state does not match layout"))
    fill!(dy, zero(eltype(dy)))
    for ((l, m), i) in L.I
        stream = k * (_Cplus(l, m, m) * _so_get(y, L.I, l + 1, m) -
                      _Cminus(l, m, m) * _so_get(y, L.I, l - 1, m))
        Pi = l == 2 ? (_so_get(y, L.I, 2, m) - sqrt(6) * _so_get(y, L.E, 2, m)) / 10 : zero(eltype(y))
        gain = l == 0 ? _so_get(y, L.I, 0, 0) :
               l == 1 && abs(m) <= 1 ? 4velocity[m + 2] : Pi
        dy[i] = -stream + kappa_dot * (-y[i] + gain)
    end
    for ((l, m), i) in L.E
        stream = k * (_Dplus(l, m, m) * _so_get(y, L.E, l + 1, m) -
                      _Dminus(l, m, m) * _so_get(y, L.E, l - 1, m) +
                      _Dzero(l, m, m) * _so_get(y, L.B, l, m))
        Pi = (_so_get(y, L.I, 2, m) - sqrt(6) * _so_get(y, L.E, 2, m)) / 10
        dy[i] = -stream + kappa_dot * (-y[i] - (l == 2 ? sqrt(6) * Pi : 0))
    end
    for ((l, m), i) in L.B
        stream = k * (_Dplus(l, m, m) * _so_get(y, L.B, l + 1, m) -
                      _Dminus(l, m, m) * _so_get(y, L.B, l - 1, m) -
                      _Dzero(l, m, m) * _so_get(y, L.E, l, m))
        dy[i] = -stream - kappa_dot * y[i]
    end
    for ((l, m), i) in L.N
        stream = k * (_Cplus(l, m, m) * _so_get(y, L.N, l + 1, m) -
                      _Cminus(l, m, m) * _so_get(y, L.N, l - 1, m))
        dy[i] = -stream
    end
    dy
end

# --- Wigner 3j and exact flat line-of-sight projectors -----------------------

"Integer Wigner 3j from the Racah sum, evaluated in scaled log-factorials."
function _wigner3j(j1::Int, j2::Int, j3::Int, m1::Int, m2::Int, m3::Int)
    any(x -> x < 0, (j1, j2, j3)) && return 0.0
    m1 + m2 + m3 == 0 || return 0.0
    abs(m1) <= j1 && abs(m2) <= j2 && abs(m3) <= j3 || return 0.0
    abs(j1 - j2) <= j3 <= j1 + j2 || return 0.0
    logdelta = _logfact(j1 + j2 - j3) + _logfact(j1 - j2 + j3) +
               _logfact(-j1 + j2 + j3) - _logfact(j1 + j2 + j3 + 1)
    lognorm = 0.5 * (logdelta +
        _logfact(j1 + m1) + _logfact(j1 - m1) +
        _logfact(j2 + m2) + _logfact(j2 - m2) +
        _logfact(j3 + m3) + _logfact(j3 - m3))
    zmin = max(0, j2 - j3 - m1, j1 - j3 + m2)
    zmax = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    zmin <= zmax || return 0.0
    logs = Float64[]
    signs = Float64[]
    for z in zmin:zmax
        ld = _logfact(z) + _logfact(j1 + j2 - j3 - z) +
             _logfact(j1 - m1 - z) + _logfact(j2 + m2 - z) +
             _logfact(j3 - j2 + m1 + z) + _logfact(j3 - j1 - m2 + z)
        push!(logs, -ld)
        push!(signs, isodd(z) ? -1.0 : 1.0)
    end
    scale = maximum(logs)
    s = sum(signs[i] * exp(logs[i] - scale) for i in eachindex(logs))
    phase = isodd(j1 - j2 - m3) ? -1.0 : 1.0
    phase * exp(lognorm + scale) * s
end

"Temperature/intensity projector J_{L,l,m}(x), thesis eq. (5.97)."
function _J_T(Ls::Int, l::Int, m::Int, x::Real)
    abs(m) <= min(Ls, l) || return 0.0
    acc = 0.0
    for l1 in abs(l - Ls):(l + Ls)
        isodd(l + l1 + Ls) && continue
        phase = isodd((l - l1 - Ls) ÷ 2) ? -1.0 : 1.0
        acc += phase * (2l1 + 1) * _wigner3j(l, l1, Ls, 0, 0, 0) *
               _wigner3j(l, l1, Ls, -m, 0, m) * sphericalbesselj(l1, float(x))
    end
    (isodd(m) ? -1.0 : 1.0) * (2l + 1) * acc
end

"Polarization projectors (J_EE=J_BB, J_EB=-J_BE), thesis eqs. (5.102-5.103)."
function _J_pol(Ls::Int, l::Int, m::Int, x::Real; cross::Bool=false)
    (Ls >= 2 && l >= 2 && abs(m) <= min(Ls, l)) || return 0.0
    acc = 0.0
    for l1 in abs(l - Ls):(l + Ls)
        d = l - l1 - Ls
        (cross ? isodd(d) : iseven(d)) || continue
        # i^d is real for even d; i^(d-1) is real for odd d.
        p = cross ? (d - 1) ÷ 2 : d ÷ 2
        phase = isodd(p) ? -1.0 : 1.0
        acc += phase * (2l1 + 1) * _wigner3j(l, l1, Ls, -2, 0, 2) *
               _wigner3j(l, l1, Ls, m, 0, -m) * sphericalbesselj(l1, float(x))
    end
    (isodd(m) ? -1.0 : 1.0) * (2l + 1) * acc
end

# --- Perturbed recombination (STZ arXiv:0812.3652 §3) -------------------------
#
# δ_e = δn_e/n_e and δT_M are first-order quantities that first matter at
# second order, where they modulate every Thomson collision term. The
# derivation ledger is derivations/second_order_implementation.tex
# §"Perturbed recombination"; equation numbers below are STZ's.
#
# The collision term Q is the SAME three-level-atom net rate the background
# solver uses (RECFAST branch: β evaluated at T_M by that convention), taken
# at displaced arguments — never a re-derived expression. Its partial
# derivatives are central differences of that one function. Helium is treated
# as already recombined over the hydrogen era evolved here (STZ include it via
# Saha; that extension belongs to the ledger before it enters the code).

"Net hydrogen recombination rate Q [m⁻³ s⁻¹] of the background solver's
three-level branch: ṅ_e|coll = Q, with C_P the Peebles escape factor and
β(T_M) the RECFAST-convention photoionization rate."
function _peebles_Q(n_e, n_H, T_M, Hz, fudge)
    x_H = clamp(n_e / n_H, 0.0, 1.0)          # H⁺ fraction (He recombined)
    n_1s = max(n_H - n_e, 0.0)
    aH = α_H(T_M, fudge)
    Adb = _debroglie(T_M)
    K_ion = Constants.T_of_wavenumber(Constants.L_H_ion)
    K_lya = Constants.T_of_wavenumber(Constants.L_H_alpha)
    βH = aH * Adb * exp(-(K_ion - K_lya) / T_M)
    β_ground = aH * Adb * exp(-K_ion / T_M)
    λ_lya = 1 / Constants.L_H_alpha
    K_H = λ_lya^3 / (8π * Hz)
    C_H = (1 + K_H * Constants.Λ_2s1s_H * n_1s) /
          (1 + K_H * (Constants.Λ_2s1s_H + βH) * n_1s)
    -C_H * (aH * n_e * x_H * n_H - β_ground * n_1s)
end

"Net HeII→HeI recombination rate Q_He [m⁻³ s⁻¹] of the background solver's
singlet channel (non-Sobolev escape branch K = λ³/(8πH); the Sobolev and
H-continuum refinements affect only the escape route of a channel that is
itself subdominant in δ_e after z ≈ 1600 — flagged in the ledger). x_HeII is
the ionized-helium fraction from the background solution."
function _peebles_Q_He(n_e, n_H, fHe, T_M, Hz)
    x_e = n_e / n_H
    x_HeII = clamp((x_e - min(x_e, 1.0)) / fHe, 0.0, 1.0)
    Tw = Constants.T_of_wavenumber
    aHe = α_He(T_M)
    Adb = _debroglie(T_M)
    K_2s_He1 = Tw(Constants.L_He1_ion - Constants.L_He_2s)
    K_2s1s_He1 = Tw(Constants.L_He_2s)
    K_2p2s_He1 = Tw(Constants.L_He_2p - Constants.L_He_2s)
    βHe = 4 * aHe * Adb * exp(-K_2s_He1 / T_M)
    β_ground_He = βHe * exp(-K_2s1s_He1 / T_M)
    He_Boltz = exp(min(K_2p2s_He1 / T_M, 680.0))
    K_He = (1 / Constants.L_He_2p)^3 / (8π * Hz)
    n_He1_eff = max(fHe * n_H * (1 - x_HeII), 0.0) * He_Boltz
    C_He = (1 + K_He * Constants.Λ_2s1s_He1 * n_He1_eff) /
           (1 + K_He * (Constants.Λ_2s1s_He1 + βHe) * n_He1_eff)
    # per-volume electron rate: fHe·nH·dx_HeII/dt
    -C_He * (aHe * n_e * x_HeII * fHe * n_H -
             β_ground_He * fHe * n_H * (1 - x_HeII))
end

"Total net recombination rate over both channels, expressed in exactly the
STZ eq.-39 argument list (n_e, n_H, T_M, H): the H⁺ and HeII fractions are
derived internally (x_H = min(x_e,1), x_HeII = (x_e−x_H)/f_He), so the
displaced-argument perturbation flows through both channels consistently."
_Q_total(n_e, n_H, fHe, T_M, Hz, fudge) =
    _peebles_Q(n_e, n_H, T_M, Hz, fudge) + _peebles_Q_He(n_e, n_H, fHe, T_M, Hz)


"Explicit-argument Saha x_e(nH, T): the same three ionization balances the
background `saha_state` closes, here as a function of its physical arguments
so the equilibrium initial condition can be perturbed at displaced arguments."
function _saha_xe_of(nH, T, fHe)
    A = _debroglie(T)
    S_H = A * exp(-Constants.T_of_wavenumber(Constants.L_H_ion) / T)
    S_He1 = 4A * exp(-Constants.T_of_wavenumber(Constants.L_He1_ion) / T)
    S_He2 = 1A * exp(-Constants.T_of_wavenumber(Constants.L_He2_ion) / T)
    f(xe) = begin
        ne = max(xe, 1e-12) * nH
        xH = S_H / (S_H + ne)
        r1 = S_He1 / ne
        r2 = S_He2 / ne
        xHeII = r1 / (1 + r1 + r1 * r2)
        xHeIII = r1 * r2 / (1 + r1 + r1 * r2)
        xH + fHe * (xHeII + 2xHeIII) - xe
    end
    lo, hi = 1e-10, 1 + 2fHe
    for _ in 1:200
        mid = 0.5 * (lo + hi)
        (f(lo) * f(mid) <= 0) ? (hi = mid) : (lo = mid)
        hi - lo < 1e-12 && break
    end
    0.5 * (lo + hi)
end


"Net hydrogen recombination rate [m⁻³ s⁻¹] from the SAME effective
multi-level atom (HyRec SWIFT) the background solver integrates, as an
explicit function of its physical arguments — so the perturbation partials
come from the identical rate family as the background, with no
reconciliation factors. Valid in HyRec's SWIFT window (checked by caller);
the ionized fraction is derived from the STZ argument list."
function _Q_hyrec(n_e, n_H, T_M, T_R, Hz, dω_cb, dω_bH, dNeff)
    xe = n_e / n_H
    xH = min(xe, 1.0)
    dxlna, _ = hyrec_dxHII_dlna(xe, xH, n_H * 1e-6, Hz,
        _HY_kB * T_M, _HY_kB * T_R, dω_cb, dω_bH, dNeff)
    n_H * Hz * dxlna
end

"Kompaneets coupling Λ_C = (4σ_T a_R k_B/m_e) n_e T_R⁴ (T_R − T_M) (STZ eq. 49),
in J m⁻³ s⁻¹."
function _kompaneets_Λ(n_e, T_R, T_M)
    aR = Constants.a_rad_SI
    (4 * Constants.σ_T_SI * Constants.c_SI * aR /
     (Constants.m_e_SI * Constants.c_SI^2)) * n_e * T_R^4 * (T_R - T_M) *
        Constants.k_B_SI / Constants.k_B_SI
end

"""
    PerturbedRecombination

First-order electron-density and matter-temperature perturbations for one
mode: callables `δ_e(x)`, `δT_M(x)` of x = ln a (δT_M in kelvin, transfer
normalization ℛ = 1). Built by [`perturbed_recombination`](@ref).
"""
struct PerturbedRecombination{I1,I2}
    k::Float64
    δ_e::I1
    δT_M::I2
end

"""
    perturbed_recombination(c, rec, p; z_start = 3500, z_end = 20, fudge = 1.14)

Solve the STZ system — eq. (38) for δ_e with the escape-probability
perturbation δ_H of eq. (40), and eq. (51) for δT_M — for the solved mode `p`,
as a post-processing ODE over the stored first-order solution. Starts in tight
Saha coupling (δ_e = δ_b, δT_M = ¼δ_γ T_M).
"""
function perturbed_recombination(c::Cosmology, rec::RecombinationSolution,
    p::PerturbationSolution; z_start=3500.0, z_end=20.0, fudge=1.14)
    k = p.k
    bg = BackgroundCache(c, rec)
    L = _layout(p)
    H0 = Constants.H0_in_invMpc(c.h)

    # per-instant first-order inputs from the stored mode
    function inputs(x)
        a = exp(x)
        u = p.sol(x)
        ℋ_ = bg.ℋ(x)
        ργ = Ω_γ(c) / a^4
        ρν = _Ω_or_zero(c, MasslessNeutrinos) / a^4
        ρb = Ω_b(c) / a^3
        ρc = Ω_c(c) / a^3
        G = MassiveNuGrid(L.nq == 0 ? 1 : L.nq)
        mν = Tuple(get_all_species(c, MassiveNeutrinos))
        φ, ψ, φ̇ = _metric(u, L, k, spatial_curvature_K(c), a, ℋ_, H0,
            ργ, ρν, ρb, ρc, G, mν)
        (; a, ℋ_, δ_b=u[3], θ_b=u[4], δ_γ=u[L.iγ], ψ, φ̇)
    end

    # background thermal state at redshift z
    function thermal(z)
        nH = n_H_of_z(c, z)
        xe = x_e(rec, z)
        (; nH, n_e=xe * nH, T_M=T_matter(rec, z), T_R=c.Tcmb * (1 + z),
            Hz=H_SI(c, z))
    end

    ε = 1e-5
    fHe_ = f_He(c)
    T0fid3 = (2.7255 / c.Tcmb)^3
    dω_cb = ((Ω_b(c) + Ω_c(c)) * c.h^2 - 0.14175) * T0fid3
    dω_bH = (Ω_b(c) * c.h^2 * (1 - c.Yp) - 0.02242 * (1 - 0.246738546372)) * T0fid3
    dNeff = _N_eff_of(c) - 3.046
    function rhs!(du, y, _, x)
        a = exp(x)
        z = 1 / a - 1
        δ_e, δTM = y
        In = inputs(x)
        th = thermal(z)
        n_n = th.nH * (1 + fHe_)              # H + He nuclei
        n_t0 = n_n + th.n_e
        # STZ eq. 37 defines Q⁰ from the background itself:
        # Q⁰ = (dn_e/dt + 3H n_e), evaluated by differentiating the stored
        # x_e history. This guarantees the zeroth-order balance is exact no
        # matter which solver (HyRec, Sobolev-He) produced the background;
        # the three-level functional form below supplies only the
        # displaced-argument PARTIALS for δQ.
        # widen the stencil at late times: after freeze-out dn_e/dη is a
        # small difference of slowly varying quantities, and a too-narrow
        # stencil turns interpolation noise into a percent-level Q⁰ error
        dzz = max(2.0, 5e-3 * z)
        ne_p = x_e(rec, z + dzz) * n_H_of_z(c, z + dzz)
        ne_m = x_e(rec, z - dzz) * n_H_of_z(c, z - dzz)
        dne_dt = -(1 + z) * th.Hz * (ne_p - ne_m) / (2dzz)
        Q0 = dne_dt + 3 * th.Hz * th.n_e

        # δ̇_b from continuity; δ_H per STZ eq. 40
        δ̇_b = -In.θ_b + 3 * In.φ̇                     # conformal-time derivative
        δ_H = -In.ψ - δ̇_b / (3 * In.ℋ_)
        δTR = 0.25 * In.δ_γ * th.T_R

        # δQ: displaced-argument partials of the SAME rate family the
        # background integrated. Inside HyRec's SWIFT window (the entire
        # visibility era) that is the effective multi-level atom itself —
        # strictly more general than the Peebles partials CLASS and SONG use,
        # and with the proper δT_R dependence the three-level RECFAST
        # convention cannot represent. Outside the window (helium era, deep
        # freeze-out) the three-level family stands in, and that boundary is
        # a ledger-flagged refinement, not a hidden switch.
        δTR_rel = 0.25 * In.δ_γ
        in_swift = 0.0042 < _HY_kB * th.T_R < 0.395
        if in_swift
            Qh(ne, nH, TM, TR, Hz) = _Q_hyrec(ne, nH, TM, TR, Hz, dω_cb, dω_bH, dNeff)
            dQ_ne = (Qh(th.n_e * (1 + ε), th.nH, th.T_M, th.T_R, th.Hz) -
                     Qh(th.n_e * (1 - ε), th.nH, th.T_M, th.T_R, th.Hz)) / (2ε)
            dQ_nH = (Qh(th.n_e, th.nH * (1 + ε), th.T_M, th.T_R, th.Hz) -
                     Qh(th.n_e, th.nH * (1 - ε), th.T_M, th.T_R, th.Hz)) / (2ε)
            dQ_TM = (Qh(th.n_e, th.nH, th.T_M * (1 + ε), th.T_R, th.Hz) -
                     Qh(th.n_e, th.nH, th.T_M * (1 - ε), th.T_R, th.Hz)) / (2ε)
            dQ_TR = (Qh(th.n_e, th.nH, th.T_M, th.T_R * (1 + ε), th.Hz) -
                     Qh(th.n_e, th.nH, th.T_M, th.T_R * (1 - ε), th.Hz)) / (2ε)
            dQ_H = (Qh(th.n_e, th.nH, th.T_M, th.T_R, th.Hz * (1 + ε)) -
                    Qh(th.n_e, th.nH, th.T_M, th.T_R, th.Hz * (1 - ε))) / (2ε)
            δQ = dQ_ne * δ_e + dQ_nH * In.δ_b + dQ_TM * (δTM / th.T_M) +
                 dQ_TR * δTR_rel + dQ_H * δ_H
        else
            dQ_ne = (_Q_total(th.n_e * (1 + ε), th.nH, fHe_, th.T_M, th.Hz, fudge) -
                     _Q_total(th.n_e * (1 - ε), th.nH, fHe_, th.T_M, th.Hz, fudge)) / (2ε)
            dQ_nH = (_Q_total(th.n_e, th.nH * (1 + ε), fHe_, th.T_M, th.Hz, fudge) -
                     _Q_total(th.n_e, th.nH * (1 - ε), fHe_, th.T_M, th.Hz, fudge)) / (2ε)
            dQ_TM = (_Q_total(th.n_e, th.nH, fHe_, th.T_M * (1 + ε), th.Hz, fudge) -
                     _Q_total(th.n_e, th.nH, fHe_, th.T_M * (1 - ε), th.Hz, fudge)) / (2ε)
            dQ_H = (_Q_total(th.n_e, th.nH, fHe_, th.T_M, th.Hz * (1 + ε), fudge) -
                    _Q_total(th.n_e, th.nH, fHe_, th.T_M, th.Hz * (1 - ε), fudge)) / (2ε)
            δQ = dQ_ne * δ_e + dQ_nH * In.δ_b + dQ_TM * (δTM / th.T_M) + dQ_H * δ_H
        end

        # δΛ_C likewise (same function, displaced arguments)
        Λ0 = _kompaneets_Λ(th.n_e, th.T_R, th.T_M)
        dΛ_ne = (_kompaneets_Λ(th.n_e * (1 + ε), th.T_R, th.T_M) -
                 _kompaneets_Λ(th.n_e * (1 - ε), th.T_R, th.T_M)) / (2ε)
        dΛ_TR = (_kompaneets_Λ(th.n_e, th.T_R * (1 + ε), th.T_M) -
                 _kompaneets_Λ(th.n_e, th.T_R * (1 - ε), th.T_M)) / (2ε)
        dΛ_TM = (_kompaneets_Λ(th.n_e, th.T_R, th.T_M * (1 + ε)) -
                 _kompaneets_Λ(th.n_e, th.T_R, th.T_M * (1 - ε))) / (2ε)
        δΛ = dΛ_ne * δ_e + dΛ_TR * (δTR / th.T_R) + dΛ_TM * (δTM / th.T_M)

        # conformal-time rates; a·(1 Mpc/c in s) converts the SI collision rates
        # to conformal Mpc units used by ℋ
        s_per_Mpc = Constants.Mpc_SI / Constants.c_SI
        aQ_ne = In.a * s_per_Mpc * Q0 / th.n_e
        kB = Constants.k_B_SI

        # STZ eq. 38 (divided by n_e), then to d/dx via 1/ℋ
        dδe_dη = δ̇_b + In.a * s_per_Mpc * ((In.ψ * Q0 + δQ) / th.n_e) - δ_e * aQ_ne
        # STZ eq. 51, solved for the time derivative
        aont = In.a * s_per_Mpc / n_t0
        dTM_dη = -2 * In.ℋ_ * δTM - (2 / 3) * th.T_M * In.θ_b -
                 aont * Q0 * δTM + 2 * th.T_M * In.φ̇ -
                 (2 / 3) * aont * (1.5 * th.T_M * Q0 - Λ0 / kB) *
                 (In.ψ - (th.n_e * δ_e + n_n * In.δ_b) / n_t0) +
                 (2 / 3) * aont * (δΛ / kB - 1.5 * th.T_M * δQ)

        du[1] = dδe_dη / In.ℋ_
        du[2] = dTM_dη / In.ℋ_
        nothing
    end

    x0 = log(1 / (1 + z_start))
    x1 = log(1 / (1 + z_end))
    In0 = inputs(x0)
    # Equilibrium (Saha) initial condition: at z_start the plasma is on the
    # Saha attractor, so δ_e = δ_b + δln x_e^Saha with the logarithmic
    # derivatives taken from the explicit-argument Saha closure. This kills
    # the transient a bare δ_e = δ_b IC would inject through the weakly
    # collisional helium era.
    th0 = thermal(z_start)
    fHe0 = f_He(c)
    εi = 1e-4
    x0S = _saha_xe_of(th0.nH, th0.T_R, fHe0)
    dlnxe_dlnn = (log(_saha_xe_of(th0.nH * (1 + εi), th0.T_R, fHe0)) -
                  log(_saha_xe_of(th0.nH * (1 - εi), th0.T_R, fHe0))) / (2εi)
    dlnxe_dlnT = (log(_saha_xe_of(th0.nH, th0.T_R * (1 + εi), fHe0)) -
                  log(_saha_xe_of(th0.nH, th0.T_R * (1 - εi), fHe0))) / (2εi)
    δT_over_T0 = 0.25 * In0.δ_γ
    δe0 = In0.δ_b * (1 + dlnxe_dlnn) + δT_over_T0 * dlnxe_dlnT
    y0 = [δe0, δT_over_T0 * T_matter(rec, z_start)]
    sol = solve(ODEProblem(rhs!, y0, (x0, x1)), Rodas5P(autodiff=false);
        reltol=1e-8, abstol=1e-10)
    PerturbedRecombination(k, x -> sol(clamp(x, x0, x1))[1],
        x -> sol(clamp(x, x0, x1))[2])
end

"STZ eq. 52: the exact super-horizon target δ_e/δ_b = 1 − ∂_η ln x_e/(3ℋ),
from the background recombination history alone."
function super_horizon_δe_ratio(c::Cosmology, rec::RecombinationSolution, z;
    dz=1.0)
    a = 1 / (1 + z)
    ℋmpc = a * H_Mpc(c, a)                       # conformal ℋ in 1/Mpc
    # dln x_e/dη = (dln x_e/dz)(dz/dη), with dz/dη = −(1+z)ℋ/a·a = −(1+z)ℋ
    dlnxe_dz = (log(x_e(rec, z + dz)) - log(x_e(rec, z - dz))) / (2dz)
    dlnxe_dη = -dlnxe_dz * (1 + z) * ℋmpc
    1 - dlnxe_dη / (3 * ℋmpc)
end
