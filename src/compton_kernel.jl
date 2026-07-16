"""
Exact isotropic Compton redistribution kernel.

The implementation follows Sarkar, Chluba & Lee, MNRAS 490 (2019) 3705,
eqs. (3)--(14) and (21). Photon and electron energies are in units of
`m_e*c^2`. The single-electron kernel is the exact angular integral of the
tree-level Klein-Nishina matrix element; `thermal_compton_kernel` folds it with
the relativistic Maxwell-Juttner distribution.

No Fokker-Planck expansion is made here. The functions are generic in their
floating-point type so the table generator can evaluate the cancellation-prone
CMB regime with `Double64`, while the package itself does not acquire a
DoubleFloats runtime dependency.
"""

using QuadGK: quadgk

@inline function _ck_S(x)
    ax = abs(x)
    if ax < oftype(x, 1e-8)
        return one(x) - x / 6 + 3x^2 / 40 - 5x^3 / 112 + 35x^4 / 1152
    elseif x >= 0
        y = sqrt(x)
        return asinh(y) / y
    else
        y = sqrt(-x)
        return asin(min(y, one(y))) / y
    end
end

@inline function _ck_F_over_wlambda(w, lambda, kappa)
    x = kappa^2 * lambda / w^2
    if abs(x) < oftype(x, 1e-7)
        a = kappa^2 / w^2
        # [S(x)-sqrt(1+x)]/(w*lambda), with the lambda->0 limit
        # taken before evaluation to avoid an artificial pole.
        return (-2a / 3 + lambda * a^2 / 5 - 3lambda^2 * a^3 / 28 +
                5lambda^3 * a^4 / 72) / w
    end
    (_ck_S(x) - sqrt(max(one(x) + x, zero(x)))) / (w * lambda)
end

function _ck_G(w0, w, kappa, omega0, omega, lambda_plus, lambda_minus)
    xplus = kappa^2 * lambda_plus / w^2
    xminus = kappa^2 * lambda_minus / w0^2
    core = 2 + (w - w0)^2 * (1 + omega * omega0) / (omega^2 * omega0^2)
    core += 2 * (_ck_S(xplus) / w - _ck_S(xminus) / w0)
    core += (1 + omega * omega0) *
            (_ck_F_over_wlambda(w, lambda_plus, kappa) -
             _ck_F_over_wlambda(w0, lambda_minus, kappa))
    kappa * core
end

"""Kinematic output-frequency interval for one electron momentum `p0`."""
function compton_frequency_bounds(omega0, p0)
    gamma0 = sqrt(1 + p0^2)
    wmin = (gamma0 - p0) * omega0 / (gamma0 + p0 + 2omega0)
    wc = (gamma0 + p0) * omega0 / (gamma0 - p0 + 2omega0)
    wt = gamma0 + omega0 - 1
    threshold = (1 + p0 - gamma0) / 2
    wmax = omega0 > threshold ? wt : wc
    wmin, wmax, wc
end

"""Minimum incident electron momentum able to scatter `omega0` to `omega`."""
function minimum_compton_momentum(omega0, omega)
    omega0 > 0 && omega > 0 || throw(ArgumentError("photon energies must be positive"))
    if omega <= omega0
        wcs = omega0 / (1 + 2omega0)
        omega > wcs && return zero(promote(omega0, omega)[1])
        p = (omega0 - omega) / 2 * sqrt((1 + omega * omega0) / (omega * omega0)) -
            (omega0 + omega) / 2
        return max(p, zero(p))
    end

    if omega0 < 1 / 2
        edge = omega0 / (1 - 2omega0)
        if omega <= edge
            return sqrt(max((omega - omega0 + 1)^2 - 1, zero(omega)))
        end
    elseif omega0 == 1 / 2
        return sqrt(max((omega - omega0 + 1)^2 - 1, zero(omega)))
    else
        return sqrt(max((omega - omega0 + 1)^2 - 1, zero(omega)))
    end

    (omega - omega0) / 2 * sqrt((1 + omega * omega0) / (omega * omega0)) +
        (omega0 + omega) / 2
end

"""
    single_electron_compton_kernel(omega0, omega, p0)

Exact redistribution density `P(omega0 -> omega, p0)` per unit `omega` and per
unit Thomson optical depth for an isotropic monoenergetic electron population.
"""
function single_electron_compton_kernel(omega0, omega, p0)
    omega0 > 0 && omega > 0 && p0 > 0 || return zero(promote(omega0, omega, p0)[1])
    gamma0 = sqrt(1 + p0^2)
    p2 = p0^2 + 2gamma0 * (omega0 - omega) + (omega0 - omega)^2
    p2 >= 0 || return zero(p2)
    p = sqrt(max(p2, zero(p2)))
    gamma = gamma0 + omega0 - omega
    gamma >= 1 || return zero(gamma)

    wmin, wmax, wc = compton_frequency_bounds(omega0, p0)
    tol = 64eps(float(one(promote(omega0, omega, p0)[1]))) * max(omega0, omega, one(omega))
    (omega >= wmin - tol && omega <= wmax + tol) || return zero(omega)

    lambda_plus = p0^2 + 2gamma0 * omega0 + omega0^2
    lambda_minus = p0^2 - 2gamma0 * omega + omega^2
    gp = max(gamma + p, eps(float(gamma)))
    g0p0 = gamma0 + p0
    wbar = sqrt(omega * omega0 * gp / g0p0)
    wbar0 = sqrt(omega * omega0 * g0p0 / gp)
    wI, wII = min(wc, omega0), max(wc, omega0)

    G = if omega < wI
        kappa = (p0 - p + omega0 + omega) / 2
        _ck_G(wbar0, wbar, kappa, omega0, omega, lambda_plus, lambda_minus)
    elseif omega < wII
        if p0 <= omega0
            _ck_G(omega, omega0, p0, omega0, omega, lambda_plus, lambda_minus)
        else
            kappa = (p - p0 + omega + omega0) / 2
            _ck_G(wbar, wbar0, kappa, omega0, omega, lambda_plus, lambda_minus)
        end
    else
        _ck_G(omega0, omega, p, omega0, omega, lambda_plus, lambda_minus)
    end

    ans = 3G / (8gamma0 * p0 * omega0^2)
    # Roundoff at a kinematic boundary can leave a tiny negative remnant after
    # cancellations. A materially negative kernel is a hard failure.
    scale = 3 / (8gamma0 * p0 * omega0^2)
    ans < -sqrt(eps(float(one(ans)))) * scale &&
        throw(DomainError(ans, "negative exact Compton kernel"))
    max(ans, zero(ans))
end

function _mj_norm(theta; rtol=1e-12)
    theta > 0 || throw(ArgumentError("electron temperature must be positive"))
    # t=(gamma-1)/theta removes exp(-1/theta) and resolves the nonrelativistic
    # Maxwell-Juttner peak without an asymptotic Bessel approximation.
    f(t) = begin
        gamma = 1 + theta * t
        p = sqrt(max(gamma^2 - 1, zero(gamma)))
        theta * gamma * p * exp(-t)
    end
    quadgk(f, zero(theta), oftype(theta, Inf); rtol)[1]
end

"""
    thermal_compton_kernel(omega0, omega, theta; rtol=1e-10)

Exact Klein-Nishina kernel folded with the normalized relativistic
Maxwell-Juttner electron distribution at `theta=kT_e/(m_e*c^2)`.
"""
function thermal_compton_kernel(omega0, omega, theta; rtol=1e-10)
    norm = _mj_norm(theta; rtol=min(rtol / 10, oftype(theta, 1e-12)))
    _thermal_compton_kernel_normed(omega0, omega, theta, norm; rtol)
end

function _thermal_compton_kernel_normed(omega0, omega, theta, norm; rtol=1e-10)
    pmin = minimum_compton_momentum(omega0, omega)
    tmin = (sqrt(1 + pmin^2) - 1) / theta
    f(q) = begin
        t = tmin + q
        gamma0 = 1 + theta * t
        p0 = sqrt(max(gamma0^2 - 1, zero(gamma0)))
        weight = theta * gamma0 * p0 * exp(-t)
        weight * single_electron_compton_kernel(omega0, omega, p0)
    end
    quadgk(f, zero(theta), oftype(theta, Inf); rtol)[1] / norm
end

"""Thermal detailed-balance residual; should vanish up to quadrature error."""
function compton_detailed_balance_residual(omega0, omega, theta; rtol=1e-10)
    fwd = thermal_compton_kernel(omega0, omega, theta; rtol)
    rev = thermal_compton_kernel(omega, omega0, theta; rtol)
    rhs = (omega0 / omega)^2 * exp((omega - omega0) / theta) * fwd
    (rev - rhs) / max(abs(rev), abs(rhs), eps(float(one(rev))))
end

# --- Cached exact CMB kernel -------------------------------------------------

const _COMPTON_KERNEL_TABLE = Ref{Any}(nothing)

_ck_readvec(io, ::Type{T}, n) where {T} = read!(io, Vector{T}(undef, n))

function _load_compton_kernel_table()
    _COMPTON_KERNEL_TABLE[] !== nothing && return _COMPTON_KERNEL_TABLE[]
    file = joinpath(@__DIR__, "..", "data", "thermalization", "compton_kernel.bin")
    isfile(file) || error("exact Compton cache missing: $file; regenerate it with the exact-kernel tabulator")
    tbl = open(file, "r") do io
        magic = String(read(io, 12))
        magic == "COSMICKN0001" || error("invalid exact Compton cache header: $magic")
        ntheta, nx, ny = Int.(_ck_readvec(io, Int32, 3))
        ltheta = _ck_readvec(io, Float64, ntheta)
        lx = _ck_readvec(io, Float64, nx)
        y = _ck_readvec(io, Float64, ny)
        wy = _ck_readvec(io, Float64, ny)
        logr = reshape(_ck_readvec(io, Float64, ntheta * nx * ny), ny, nx, ntheta)
        (ltheta=ltheta, lx=lx, y=y, wy=wy, logr=logr)
    end
    _COMPTON_KERNEL_TABLE[] = tbl
end

@inline function _ck_bracket(v, q)
    a = clamp(searchsortedlast(v, q), 1, length(v) - 1)
    f = clamp((q - v[a]) / (v[a+1] - v[a]), 0.0, 1.0)
    a, f
end

@inline function _ck_lograte(tbl, iy, lx, ltheta)
    ix, fx = _ck_bracket(tbl.lx, lx)
    it, ft = _ck_bracket(tbl.ltheta, ltheta)
    (1-ft) * ((1-fx) * tbl.logr[iy, ix, it] + fx * tbl.logr[iy, ix+1, it]) +
    ft * ((1-fx) * tbl.logr[iy, ix, it+1] + fx * tbl.logr[iy, ix+1, it+1])
end

@inline function _ck_deposit(x, xo)
    xo <= x[1] && return (1, 1, 1.0, 0.0)
    xo >= x[end] && return (length(x), length(x), 1.0, 0.0)
    j = clamp(searchsortedlast(x, xo), 1, length(x) - 1)
    f = (xo - x[j]) / (x[j+1] - x[j])
    j, j + 1, 1 - f, f
end

function _ck_volumes(x)
    nx = length(x)
    edge = sqrt.(x[1:end-1] .* x[2:end])
    vol = zeros(eltype(x), nx)
    vol[1] = (edge[1] - x[1]) * x[1]^2
    vol[end] = (x[end] - edge[end]) * x[end]^2
    @inbounds for i in 2:nx-1
        vol[i] = (edge[i] - edge[i-1]) * x[i]^2
    end
    vol
end

mutable struct ExactComptonOperator
    table
    x::Vector{Float64}
    lnx::Vector{Float64}
    npl::Vector{Float64}
    bose::Vector{Float64}
    vol::Vector{Float64}
    eweight::Vector{Float64}
    hdelta::Vector{Float64}
    Cdelta::Vector{Float64}
    Ctemp::Vector{Float64}
    Jh::Matrix{Float64}
    Jdelta::Matrix{Float64}
end

function ExactComptonOperator(x::Vector{Float64})
    npl = @. 1 / expm1(x)
    bose = @. npl * (1 + npl)
    z = zeros(length(x))
    vol = _ck_volumes(x)
    ExactComptonOperator(_load_compton_kernel_table(), x, log.(x), npl, bose,
        vol, vol .* x, copy(z), copy(z), copy(z), zeros(length(x), length(x)),
        zeros(length(x), length(x)))
end

function _ck_apply!(out, h, op::ExactComptonOperator, theta)
    tbl, x, lnx, npl = op.table, op.x, op.lnx, op.npl
    lt = log(theta)
    root = sqrt(theta)
    fill!(out, 0.0)
    @inbounds for i in eachindex(x)
        xi = x[i]
        lxi = lnx[i]
        hi = h[i]
        for iy in eachindex(tbl.y)
            s = root * tbl.y[iy]
            xo = xi * exp(s)
            no = xo < 700 ? 1 / expm1(xo) : 0.0
            rate = tbl.wy[iy] * exp(_ck_lograte(tbl, iy, lxi, lt))
            j1, j2, f1, f2 = _ck_deposit(x, xo)

            # Use the detailed-balance weak form. The destination is one
            # linearly interpolated virtual state, h_o=f1*h[j1]+f2*h[j2], not
            # two full-cell jumps. With v=phi(x_o)-phi(x_i), every directed
            # transition contributes -a*v*v' to the discrete operator. Hence
            # sum(vol*C)=0 and sum(vol*h*C)<=0 exactly, while h=x resolves the
            # true sub-cell shift x_o-x_i instead of acquiring grid diffusion.
            # No Fokker-Planck approximation is introduced; a is the cached
            # exact-kernel weight and the 1/2 averages the forward/reverse
            # representations that coincide in the continuum.
            a = 0.5 * op.vol[i] * npl[i] * rate * (1 + no)
            d = a * (f1 * h[j1] + f2 * h[j2] - hi)
            out[i] += d / op.vol[i]
            out[j1] -= f1 * d / op.vol[j1]
            out[j2] -= f2 * d / op.vol[j2]
        end
    end
    out
end

"""Assemble the exact weak-form map from nodal `h` to collision `C`."""
function _ck_matrix!(A, op::ExactComptonOperator, theta)
    tbl, x, lnx, npl = op.table, op.x, op.lnx, op.npl
    lt = log(theta)
    root = sqrt(theta)
    fill!(A, 0.0)
    @inbounds for i in eachindex(x)
        xi = x[i]
        lxi = lnx[i]
        for iy in eachindex(tbl.y)
            s = root * tbl.y[iy]
            xo = xi * exp(s)
            no = xo < 700 ? 1 / expm1(xo) : 0.0
            rate = tbl.wy[iy] * exp(_ck_lograte(tbl, iy, lxi, lt))
            j1, j2, f1, f2 = _ck_deposit(x, xo)
            a = 0.5 * op.vol[i] * npl[i] * rate * (1 + no)
            ids = (i, j1, j2)
            vals = (-1.0, f1, f2)
            # M*C = -a*v*v'*h, with diagonal mass M=diag(vol).
            for r in 1:3, q in 1:3
                A[ids[r], ids[q]] -= a * vals[r] * vals[q] / op.vol[ids[r]]
            end
        end
    end
    A
end

"""Assemble the exact linear map from occupation perturbation `delta_n` to C."""
function _ck_delta_matrix!(op::ExactComptonOperator, theta)
    _ck_matrix!(op.Jh, op, theta)
    @inbounds for j in axes(op.Jh, 2), i in axes(op.Jh, 1)
        op.Jdelta[i, j] = op.Jh[i, j] / op.bose[j]
    end
    op.Jdelta
end

"""
    exact_compton_collision!(out, delta_n, op, theta, fractional_heat_per_tau)

Linearized exact Compton collision term per Thomson optical depth. The returned
electron-temperature excess `phi=(T_e-T_gamma)/T_gamma` is chosen so the photon
energy gain equals `fractional_heat_per_tau`. Photon number is conserved by
construction and a Planck/Bose-Einstein spectrum is stationary.
"""
function exact_compton_collision!(out, delta_n, op::ExactComptonOperator, theta,
    fractional_heat_per_tau=0.0; emission_source=nothing, emission_sink=nothing)
    @. op.hdelta = delta_n / op.bose
    _ck_apply!(op.Cdelta, op.hdelta, op, theta)
    _ck_apply!(op.Ctemp, op.x, op, theta)

    E_delta = dot(op.vol .* op.x, op.Cdelta)
    E_temp = dot(op.vol .* op.x, op.Ctemp)
    I3 = dot(op.vol .* op.x, op.npl)
    target = fractional_heat_per_tau * I3
    # If photon production is active, its phi-dependent source and distortion
    # sink share the same electron energy ledger. Solving the combined balance
    # prevents double-counting injected heat between Compton redistribution and
    # DC/BR emission. Arrays are per-optical-depth coefficients.
    E_source = emission_source === nothing ? 0.0 :
        dot(op.vol .* op.x, emission_source)
    E_sink = emission_sink === nothing ? 0.0 :
        dot(op.vol .* op.x, emission_sink .* delta_n)
    denom = E_temp - E_source
    abs(denom) > eps(Float64) || error("degenerate exact electron-temperature mode")
    phi = (E_delta - E_sink - target) / denom
    @. out = op.Cdelta - phi * op.Ctemp
    phi
end

"""
    _exact_collision_emission_jacobian!(J, penergy, op, theta, source, sink)

Exact derivative of `C_KN + phi*source - sink*delta_n` with respect to the
occupation perturbation. `phi` obeys the same shared electron-energy ledger as
`exact_compton_collision!`.
"""
function _exact_collision_emission_jacobian!(J, penergy, op::ExactComptonOperator,
    theta, emission_source, emission_sink)
    L = _ck_delta_matrix!(op, theta)
    _ck_apply!(op.Ctemp, op.x, op, theta)
    eweight = op.eweight
    E_temp = dot(eweight, op.Ctemp)
    E_source = dot(eweight, emission_source)
    denom = E_temp - E_source
    abs(denom) > eps(Float64) || error("degenerate exact electron-temperature Jacobian")

    @inbounds for j in eachindex(penergy)
        E_column = 0.0
        for i in eachindex(eweight)
            E_column += eweight[i] * L[i, j]
        end
        penergy[j] = (E_column - eweight[j] * emission_sink[j]) / denom
    end
    @inbounds for j in axes(J, 2), i in axes(J, 1)
        J[i, j] = L[i, j] +
            (emission_source[i] - op.Ctemp[i]) * penergy[j]
    end
    @inbounds for i in axes(J, 1)
        J[i, i] -= emission_sink[i]
    end
    J
end
