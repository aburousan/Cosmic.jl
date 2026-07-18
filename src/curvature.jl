"""
Exact constant-curvature geometry used by the perturbation and transfer layers.

We use the CLASS/Hu--Seljak--White--Zaldarriaga convention

    K = -Omega_k H0^2,

so K>0 is a closed three-sphere and K<0 an open three-hyperboloid.  The
hyperspherical Bessel functions below solve the radial Laplacian eigenproblem
itself.  In particular, this file intentionally has no large-nu replacement by
ordinary spherical Bessel functions: weak curvature is recovered only as the
continuous K -> 0 limit of the curved equations.
"""

@inline spatial_curvature_K(c::Cosmology) =
    -Ω_k(c) * Constants.H0_in_invMpc(c.h)^2

"""
    _curvature_streaming(K, k, ell, m=0)

The free-streaming coefficient `sqrt(q^2 - K*ell^2)/k` of the optimal hierarchy
(arXiv:2005.12119 eq. 20), whose recursion coefficients are the m-independent
`0-kappa_j^0`.  The eigenvalue relation is m-dependent, however:

    q^2 = k^2 + (1 + |m|) K,

so the harmonic wavenumber `k` handed in by the scalar (m=0), vector (m=1) and
tensor (m=2) drivers refers to a different `q` in each case.  Passing `m` is
what keeps that conversion exact; hard-coding the scalar relation costs 2K for
tensors, which is a double-digit error on the coefficient at the smallest q.
"""
@inline function _curvature_streaming(K, k, ell, m::Integer=0)
    r = 1 + K * ((1 + abs(m)) - ell^2) / k^2
    # Closed eigenmodes have a finite angular tower.  Roundoff at its endpoint
    # can make r a tiny negative number; a genuinely negative value means that
    # this multipole does not exist for the selected eigenmode.
    sqrt(max(r, 0.0))
end

@inline function _curvature_cot_closure(K, k, eta)
    if iszero(K)
        return inv(k * eta)
    end
    rootK = sqrt(abs(K))
    chi = rootK * eta
    K > 0 ? rootK / (k * tan(chi)) : rootK / (k * tanh(chi))
end

@inline _sinK(sgnK, x) = sgnK == 1 ? sin(x) : sgnK == -1 ? sinh(x) : x
@inline _cosK(sgnK, x) = sgnK == 1 ? cos(x) : sgnK == -1 ? cosh(x) : 1.0
@inline _cotK(sgnK, x) = _cosK(sgnK, x) / _sinK(sgnK, x)
@inline _mode_root(sgnK, beta, ell) = sqrt(max(beta^2 - sgnK * ell^2, 0.0))

"""
    _hyperspherical_bessels(sgnK, beta, lmax, chi)

Return `(Phi, dPhi, d2Phi)` for all `0:lmax`, where derivatives are with
respect to the dimensionless radial coordinate `chi=sqrt(abs(K))*r` and

    Phi'' + 2 cot_K(chi) Phi'
         + [beta^2 - sgnK - ell(ell+1)/sin_K(chi)^2] Phi = 0.

`beta=q/sqrt(abs(K))`.  Below the turning point we use Miller's stable backward
recurrence, normalized by the analytic monopole.  Above it we use the forward
recurrence.  In a closed model beta is an integer and the tower terminates at
ell=beta-1 exactly.
"""
function _hyperspherical_bessels(sgnK::Int, beta::Real, lmax::Int, chi::Real)
    sgnK in (-1, 0, 1) || throw(ArgumentError("sgnK must be -1, 0, or 1"))
    beta > 0 || throw(ArgumentError("hyperspherical beta must be positive"))
    lmax >= 0 || throw(ArgumentError("lmax must be nonnegative"))
    x = float(chi)
    b = float(beta)
    Phi = zeros(lmax + 1)
    # Closed eigenmodes terminate at ell=beta-1.  Zero initialisation makes
    # that exact termination explicit for the derivative arrays as well.
    dPhi = zeros(lmax + 1)
    d2Phi = zeros(lmax + 1)

    if sgnK == 0
        for ell in 0:lmax
            z = b * x
            Phi[ell+1] = sphericalbesselj(ell, z)
            if z == 0
                dPhi[ell+1] = ell == 1 ? b / 3 : 0.0
                d2Phi[ell+1] = ell == 0 ? -b^2 / 3 : ell == 2 ? 2b^2 / 15 : 0.0
            else
                dPhi[ell+1] = b * (ell * sphericalbesselj(ell, z) / z -
                                          sphericalbesselj(ell + 1, z))
                d2Phi[ell+1] = -2dPhi[ell+1] / x -
                    (b^2 - ell * (ell + 1) / x^2) * Phi[ell+1]
            end
        end
        return Phi, dPhi, d2Phi
    end

    if x == 0
        Phi[1] = 1.0
        lmax >= 1 && (dPhi[2] = _mode_root(sgnK, b, 1) / 3)
        d2Phi[1] = -(b^2 - sgnK) / 3
        lmax >= 2 && (d2Phi[3] = 2 * _mode_root(sgnK, b, 1) *
            _mode_root(sgnK, b, 2) / 15)
        return Phi, dPhi, d2Phi
    end

    # Frobenius solution at the regular origin.  Recurrences determine a tiny
    # Phi_l by cancellation when beta*chi << 1; the spin-2 projectors then divide
    # by sin_K(chi)^2 and magnify that roundoff catastrophically.  Keeping the
    # first curvature-dependent correction gives relative error O((beta*chi)^4)
    # and joins smoothly onto the recurrence well before the turning point.
    if abs(b * x) < 0.1
        coeff = 1.0
        for ell in 0:lmax
            if ell > 0
                coeff *= _mode_root(sgnK, b, ell) / (2ell + 1)
            end
            a2 = -(b^2 - sgnK * (ell^2 + 3ell + 3) / 3) / (2 * (2ell + 3))
            xl = x^ell
            Phi[ell+1] = coeff * xl * (1 + a2 * x^2)
            dPhi[ell+1] = coeff *
                (ell == 0 ? 2a2 * x : ell * x^(ell-1) + a2 * (ell + 2) * x^(ell+1))
            lead2 = ell < 2 ? 0.0 : ell * (ell - 1) * x^(ell-2)
            d2Phi[ell+1] = coeff * (lead2 + a2 * (ell + 2) * (ell + 1) * xl)
        end
        return Phi, dPhi, d2Phi
    end

    sinx = _sinK(sgnK, x)
    cotx = _cotK(sgnK, x)
    # sinc(y)=sin(pi*y)/(pi*y), hence this ratio is stable at the origin.
    phi0 = (x / sinx) * sinc(b * x / pi)
    max_nonzero = sgnK == 1 ? min(lmax, max(round(Int, b) - 1, 0)) : lmax
    if sgnK == 1
        abs(b - round(b)) <= 2e-10 * max(b, 1) || throw(ArgumentError(
            "closed hyperspherical modes require integer beta; got beta=$b"))
    end

    # The l-recurrence changes stability at beta*sin_K(chi) ~ sqrt(l(l+1)).
    forward = b * abs(sinx) > sqrt((max_nonzero + 1) * (max_nonzero + 2.0))
    if forward
        Phi[1] = phi0
        if max_nonzero >= 1
            Phi[2] = phi0 * (cotx - b / tan(b * x)) / _mode_root(sgnK, b, 1)
            for ell in 2:max_nonzero
                Phi[ell+1] = ((2ell - 1) * cotx * Phi[ell] -
                    _mode_root(sgnK, b, ell - 1) * Phi[ell-1]) /
                    _mode_root(sgnK, b, ell)
            end
        end
    else
        # Miller recurrence.  In the closed case Phi_beta=0 is the exact upper
        # boundary.  In the open case a buffer of 64 orders is ample below the
        # turning point and is checked against a larger buffer in the tests.
        top = sgnK == 1 ? max(round(Int, b) - 1, max_nonzero) : max_nonzero + 64
        work = zeros(top + 2) # indices ell+1, with an extra Phi_(top+1)
        work[top+1] = 1.0
        for ell in top:-1:1
            work[ell] = ((2ell + 1) * cotx * work[ell+1] -
                _mode_root(sgnK, b, ell + 1) * work[ell+2]) /
                _mode_root(sgnK, b, ell)
            # Ratios, not the arbitrary Miller normalization, carry physics.
            # Keep the recursion away from overflow without changing them.
            if abs(work[ell]) > 1e150
                work[ell:end] .*= 1e-150
            end
        end
        scale = phi0 / work[1]
        for ell in 0:max_nonzero
            Phi[ell+1] = scale * work[ell+1]
        end
    end

    # Use the lowering identity for derivatives.  Unlike the equivalent
    # raising identity it does not require an unreturned Phi_(lmax+1), and it
    # remains valid at the finite endpoint of a closed eigenmode.
    dPhi[1] = abs(b * x) < 1e-3 ? -(b^2 - sgnK) * x / 3 :
        cos(b * x) / sinx - cotx * Phi[1]
    for ell in 1:max_nonzero
        dPhi[ell+1] = _mode_root(sgnK, b, ell) * Phi[ell] -
            (ell + 1) * cotx * Phi[ell+1]
    end
    for ell in 0:max_nonzero
        d2Phi[ell+1] = -2cotx * dPhi[ell+1] +
            (ell * (ell + 1) / sinx^2 - b^2 + sgnK) * Phi[ell+1]
    end
    Phi, dPhi, d2Phi
end
