"""
Second-order Einstein--Boltzmann angular machinery.

This file contains only pieces that are exact and independently testable before
the full second-order solver is exposed: Fourier-triangle geometry, rotations of
first-order scalar multipoles, the complete (ell,m) hierarchy layout, angular
couplings, the linear-in-second-order radiation operator, and the general
line-of-sight projection kernels.  The quadratic sources are deliberately not
represented by zero defaults: the production driver is added only after every
Einstein, Liouville, collision, matter and recombination source has a reference
gate.  This prevents a partial hierarchy from being mistaken for a physical
second-order calculation.

Conventions and equations are those of Pitrou (arXiv:0809.3245) and Pettinari's
thesis (arXiv:1405.2280), chapters 5--6 and appendices A--B.  SONG is used only
as an independent numerical reference at matched flat, adiabatic settings.
"""

using SpecialFunctions: loggamma, sphericalbesselj

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
