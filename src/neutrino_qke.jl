"""
Momentum-dependent neutrino flavour-density matrices.

The transport equation used by FortEPiaNO, NEVO, and the precision
neutrino-decoupling calculations is

    dϱ(x,y)/dx = -i A(x) [H_eff(x,y), ϱ(x,y)] + C[ϱ](x,y),

where `x = mₑ a`, `y = p a`, `z = Tγ a`, and `ϱ` is Hermitian in flavour
space.  `H_eff` contains the vacuum, charged-lepton, and neutrino mean-field
terms; `C` is the full weak collision operator.  This file deliberately does
not replace `C` by a hand-tuned relaxation rate.  It reads the public
FortEPiaNO serialization of a solved QKE and exposes it through a checked,
backend-independent Cosmic API.

References:

* Mangano et al., Nucl. Phys. B729 (2005), arXiv:hep-ph/0506164.
* Froustey, Pitrou & Volpe, JCAP 12 (2020) 015, arXiv:2008.01074.
* Bennett et al., JCAP 04 (2021) 073, arXiv:2012.02726.

FortEPiaNO evolves zero primordial lepton asymmetry.  A history loaded through
this interface must therefore not be presented as a solution for unequal
neutrino and antineutrino density matrices.
"""

using DelimitedFiles: readdlm
using LinearAlgebra: Hermitian, SymTridiagonal, diag, eigen, eigvals
using QuadGK: quadgk

"""
    FlavorDensityHistory

A momentum-dependent, Hermitian neutrino density-matrix history.

Array dimensions are `rho[time, momentum, flavour, flavour]`.  `x = mₑa`,
`y = pa`, `z = Tγa`, and `w` is FortEPiaNO's effective comoving neutrino
temperature.  `mass_occupation`, when present, contains the instantaneous
matter/mass-basis diagonal occupations with dimensions
`[time, momentum, state]`.  `momentum_weights`, when present, are the
generalized Gauss--Laguerre weights for directly integrating `y³ ϱ` on the
stored grid.
"""
struct FlavorDensityHistory{T<:AbstractFloat}
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    w::Vector{T}
    rho::Array{Complex{T},4}
    mass_occupation::Union{Nothing,Array{T,3}}
    neff::Union{Nothing,Vector{T}}
    momentum_weights::Union{Nothing,Vector{T}}
    source::String
end

Base.length(h::FlavorDensityHistory) = length(h.x)
flavor_count(h::FlavorDensityHistory) = size(h.rho, 3)

const _STANDARD_QKE_DIR = joinpath(
    @__DIR__, "..", "data", "neutrino_qke", "standard_fortepiano")

"""
    standard_qke_history(; validate=true)

Load Cosmic's bundled Standard-Model three-flavour FortEPiaNO reference
history.  This is a solved, momentum-dependent zero-asymmetry QKE history,
not a fitted relaxation model.  Pass a different output directory to
[`load_flavor_density_history`](@ref) for other flavour counts, mixing
parameters, grids, or collision configurations.
"""
standard_qke_history(; validate=true) =
    load_flavor_density_history(_STANDARD_QKE_DIR; validate)

function _qke_table(file)
    isfile(file) || throw(ArgumentError("missing QKE table: $file"))
    raw = readdlm(file, comments=true, comment_char='#')
    ndims(raw) == 2 || throw(ArgumentError("QKE table is not rectangular: $file"))
    Matrix{Float64}(raw)
end

function _qke_aligned_values(file, x, ny)
    d = _qke_table(file)
    expected = (length(x), ny + 1)
    size(d) == expected ||
        throw(DimensionMismatch("$file has size $(size(d)); expected $expected"))
    all(isapprox.(d[:, 1], x; rtol=2e-9, atol=2e-12)) ||
        throw(ArgumentError("time grid in $file does not match z.dat"))
    d[:, 2:end]
end

function _qke_flavor_count(path)
    n = 0
    while isfile(joinpath(path, "nuDens_diag$(n + 1).dat"))
        n += 1
    end
    n > 0 || throw(ArgumentError("no nuDens_diag*.dat files found in $path"))
    n
end

function _qke_gauss_laguerre_weights(path, y)
    log_file = joinpath(path, "messages.log")
    isfile(log_file) || return nothing
    pattern = r"Gauss-Laguerre,\s*n=\s*(\d+)\s+and selecting the first\s+(\d+)\s+roots"
    for line in eachline(log_file)
        m = match(pattern, line)
        m === nothing && continue
        order = parse(Int, m.captures[1])
        selected = parse(Int, m.captures[2])
        selected == length(y) || continue
        order >= selected ||
            throw(ArgumentError("invalid Gauss--Laguerre order in $log_file"))

        # FortEPiaNO configures get_GLq_vectors(..., alpha=3), so its raw
        # occupation tables are integrated as ∫dy y³ f(y).  Golub--Welsch for
        # the generalized Laguerre weight y^α exp(-y) gives the same nodes and
        # weights as FortEPiaNO's gaulag routine (Numerical Recipes).
        α = 3.0
        d = [2k - 1 + α for k in 1:order]
        e = [sqrt(k * (k + α)) for k in 1:order-1]
        decomposition = eigen(SymTridiagonal(d, e))
        nodes = decomposition.values[1:selected]
        all(isapprox.(nodes, y; rtol=2e-8, atol=2e-10)) ||
            throw(ArgumentError(
                "Gauss--Laguerre metadata in $log_file does not reproduce y_grid.dat"))
        standard_weights = 6 .* decomposition.vectors[1, :].^2 # Γ(α+1)=Γ(4)=6
        return standard_weights[1:selected] .* exp.(nodes)
    end
    nothing
end

"""
    load_flavor_density_history(path; require_coherence=true, validate=true)

Read a FortEPiaNO output directory without assuming a flavour count, momentum
grid, or time grid.  The lower triangle is reconstructed by Hermitian
conjugation from FortEPiaNO's serialized upper triangle.

When `require_coherence=true`, every off-diagonal real and imaginary table must
be present.  Set it to `false` only for a calculation that was intentionally
run without coherences.  With `validate=true`, every matrix is checked for
Hermiticity, positive semidefiniteness, and the fermionic occupation bound
`ϱ ≤ I`.
"""
function load_flavor_density_history(path::AbstractString;
    require_coherence=true, validate=true, atol=2e-7)

    root = abspath(path)
    ztab = _qke_table(joinpath(root, "z.dat"))
    size(ztab, 2) >= 2 ||
        throw(DimensionMismatch("z.dat must contain at least x and z"))
    x = ztab[:, 1]
    z = ztab[:, 2]
    w = size(ztab, 2) >= 3 ? ztab[:, 3] : ones(length(x))
    all(diff(x) .> 0) || throw(ArgumentError("QKE time grid x must be strictly increasing"))
    all(x .> 0) && all(z .> 0) && all(w .> 0) ||
        throw(DomainError((x, z, w), "QKE temperatures and scale variable must be positive"))

    ytab = _qke_table(joinpath(root, "y_grid.dat"))
    size(ytab, 2) == 1 ||
        throw(DimensionMismatch("y_grid.dat must contain one momentum column"))
    y = vec(ytab)
    all(diff(y) .> 0) && all(y .> 0) ||
        throw(ArgumentError("QKE momentum grid y must be positive and strictly increasing"))

    nf = _qke_flavor_count(root)
    nx, ny = length(x), length(y)
    rho = zeros(ComplexF64, nx, ny, nf, nf)
    for i in 1:nf
        rho[:, :, i, i] .= _qke_aligned_values(
            joinpath(root, "nuDens_diag$i.dat"), x, ny)
    end
    for i in 1:nf-1, j in i+1:nf
        rf = joinpath(root, "nuDens_nd_$(i)$(j)_re.dat")
        imf = joinpath(root, "nuDens_nd_$(i)$(j)_im.dat")
        if isfile(rf) && isfile(imf)
            re = _qke_aligned_values(rf, x, ny)
            im = _qke_aligned_values(imf, x, ny)
            rho[:, :, i, j] .= complex.(re, im)
            rho[:, :, j, i] .= complex.(re, .-im)
        elseif require_coherence
            throw(ArgumentError("missing coherence tables for flavours $i,$j in $root"))
        end
    end

    mass_files = [joinpath(root, "nuDens_mass$i.dat") for i in 1:nf]
    mass = if all(isfile, mass_files)
        a = zeros(Float64, nx, ny, nf)
        for i in 1:nf
            a[:, :, i] .= _qke_aligned_values(mass_files[i], x, ny)
        end
        a
    elseif any(isfile, mass_files)
        throw(ArgumentError("incomplete set of nuDens_mass*.dat files in $root"))
    else
        nothing
    end

    neff_file = joinpath(root, "Neff.dat")
    neff = if isfile(neff_file)
        n = _qke_table(neff_file)
        size(n, 1) == nx && size(n, 2) >= 2 ||
            throw(DimensionMismatch("Neff.dat is not aligned with z.dat"))
        all(isapprox.(n[:, 1], x; rtol=2e-9, atol=2e-12)) ||
            throw(ArgumentError("time grid in Neff.dat does not match z.dat"))
        n[:, end]
    else
        nothing
    end

    momentum_weights = _qke_gauss_laguerre_weights(root, y)
    h = FlavorDensityHistory(x, y, z, w, rho, mass, neff,
        momentum_weights, root)
    validate && validate_flavor_density_history(h; atol)
    h
end

"""
    validate_flavor_density_history(history; atol=2e-7)

Check Hermiticity, positivity, and Pauli occupancy bounds for every stored
density matrix.  Returns the extremal eigenvalues and maximum Hermiticity
residual.  A violation beyond `atol` throws `DomainError`.
"""
function validate_flavor_density_history(h::FlavorDensityHistory; atol=2e-7)
    nx, ny, nf, nf2 = size(h.rho)
    nf == nf2 || throw(DimensionMismatch("density matrices must be square"))
    (nx, ny) == (length(h.x), length(h.y)) ||
        throw(DimensionMismatch("density-matrix axes do not match x and y"))

    min_eig = Inf
    max_eig = -Inf
    hermitian_residual = 0.0
    for ix in 1:nx, iy in 1:ny
        r = Matrix(@view h.rho[ix, iy, :, :])
        hermitian_residual = max(hermitian_residual, maximum(abs, r - r'))
        λ = eigvals(Hermitian((r + r') / 2))
        min_eig = min(min_eig, minimum(λ))
        max_eig = max(max_eig, maximum(λ))
    end
    hermitian_residual <= atol ||
        throw(DomainError(hermitian_residual, "QKE density matrix is not Hermitian"))
    min_eig >= -atol ||
        throw(DomainError(min_eig, "QKE density matrix is not positive semidefinite"))
    max_eig <= 1 + atol ||
        throw(DomainError(max_eig, "QKE density matrix violates the fermionic occupancy bound"))

    if h.mass_occupation !== nothing
        minimum(h.mass_occupation) >= -atol ||
            throw(DomainError(minimum(h.mass_occupation),
                "mass-basis occupation is negative"))
        maximum(h.mass_occupation) <= 1 + atol ||
            throw(DomainError(maximum(h.mass_occupation),
                "mass-basis occupation exceeds one"))
    end
    (; min_eigenvalue=min_eig, max_eigenvalue=max_eig,
        hermitian_residual)
end

function _qke_bracket(v, q, name)
    first(v) <= q <= last(v) ||
        throw(DomainError(q, "$name lies outside [$(first(v)), $(last(v))]"))
    i = searchsortedlast(v, q)
    i == length(v) && return (i, i, 0.0)
    i == 0 && return (1, 1, 0.0)
    t = (q - v[i]) / (v[i + 1] - v[i])
    (i, i + 1, t)
end

"""
    flavor_density_matrix(history, x, y)

Bilinearly interpolate the QKE density matrix at `(x,y)`.  Extrapolation is
rejected because neither coherent evolution nor collision freeze-out is
generically linear outside the solved domain.
"""
function flavor_density_matrix(h::FlavorDensityHistory, x, y)
    i0, i1, tx = _qke_bracket(h.x, x, "x")
    j0, j1, ty = _qke_bracket(h.y, y, "y")
    r00 = Matrix(@view h.rho[i0, j0, :, :])
    r01 = Matrix(@view h.rho[i0, j1, :, :])
    r10 = Matrix(@view h.rho[i1, j0, :, :])
    r11 = Matrix(@view h.rho[i1, j1, :, :])
    r0 = (1 - ty) * r00 + ty * r01
    r1 = (1 - ty) * r10 + ty * r11
    r = (1 - tx) * r0 + tx * r1
    (r + r') / 2
end

"Physical photon temperature in MeV at QKE time `x`."
function qke_photon_temperature(h::FlavorDensityHistory, x)
    i0, i1, t = _qke_bracket(h.x, x, "x")
    z = (1 - t) * h.z[i0] + t * h.z[i1]
    _ME_MEV * z / x
end

function qke_scale_variable(h::FlavorDensityHistory, Tγ)
    Tγ > 0 || throw(DomainError(Tγ, "photon temperature must be positive"))
    temperatures = _ME_MEV .* h.z ./ h.x
    all(diff(temperatures) .< 0) ||
        throw(ArgumentError("QKE photon-temperature history must decrease monotonically"))
    if Tγ >= temperatures[1]
        return h.x[1] * temperatures[1] / Tγ
    elseif Tγ <= temperatures[end]
        return h.x[end] * temperatures[end] / Tγ
    end
    ascending_T = reverse(temperatures)
    ascending_x = reverse(h.x)
    exp(_qke_linear(log.(ascending_T), log.(ascending_x), log(Tγ)))
end

"Scale factor normalized to `a=1` at the supplied present CMB temperature."
function qke_scale_factor(h::FlavorDensityHistory, Tγ; Tcmb=2.7255)
    T0 = Tcmb * Constants.K_to_eV * 1e-6
    qke_scale_variable(h, Tγ) / qke_scale_variable(h, T0)
end

"Effective neutrino temperature `mₑ w/x` in MeV at photon temperature `Tγ`."
function qke_neutrino_temperature(h::FlavorDensityHistory, Tγ)
    x = qke_scale_variable(h, Tγ)
    w = if x <= h.x[1]
        h.w[1]
    elseif x >= h.x[end]
        h.w[end]
    else
        _qke_linear(h.x, h.w, x)
    end
    _ME_MEV * w / x
end

"""
    qke_occupation(history, state, energy, Tγ; basis=:flavour)

Return the solved occupation at physical momentum/energy `energy` (MeV) and
photon temperature `Tγ` (MeV).  For relativistic neutrinos,
`y = p a = (energy/mₑ)x`.  Outside the solved temperature interval the first
or last spectrum is free-streamed with `y` fixed; momentum tails use the
positive exponential continuation checked by `_qke_profile`.
"""
function qke_occupation(h::FlavorDensityHistory, state::Integer, energy, Tγ;
    basis=:flavour)

    1 <= state <= flavor_count(h) ||
        throw(BoundsError(1:flavor_count(h), state))
    energy >= 0 || throw(DomainError(energy, "neutrino energy must be non-negative"))
    x = qke_scale_variable(h, Tγ)
    ix0, ix1, tx = if x <= h.x[1]
        (1, 1, 0.0)
    elseif x >= h.x[end]
        (length(h.x), length(h.x), 0.0)
    else
        _qke_bracket(h.x, x, "x")
    end
    values = if basis === :flavour
        (1 - tx) .* collect(real.(@view h.rho[ix0, :, state, state])) .+
            tx .* collect(real.(@view h.rho[ix1, :, state, state]))
    elseif basis === :mass
        h.mass_occupation === nothing &&
            throw(ArgumentError("mass-basis occupations were not saved"))
        (1 - tx) .* collect(@view h.mass_occupation[ix0, :, state]) .+
            tx .* collect(@view h.mass_occupation[ix1, :, state])
    else
        throw(ArgumentError("basis must be :mass or :flavour"))
    end
    y = energy * x / _ME_MEV
    first(_qke_profile(h.y, values, y))
end

function _qke_neff_integral(h::FlavorDensityHistory, x)
    integral = if h.momentum_weights === nothing
        # Interpolate the diagonal occupations onto the node grid at time x, then
        # integrate y³ Trϱ by reconstructing the spectrum from its slowly varying
        # distortion δ = Trϱ/(n_f f_FD) − 1. Integrating flavor_density_matrix
        # directly interpolates the raw occupation between the sparse Gauss--
        # Laguerre nodes, which overestimates ∫y³f at the percent level (the same
        # reason _qke_profile exists). The distortion profile is regular at y=0
        # and holds δ fixed in the exponentially suppressed tail, so a single
        # quadrature over [0, y_tail] needs no separate origin or tail term.
        ix0, ix1, tx = _qke_bracket(h.x, x, "x")
        nf = size(h.rho, 3)
        tr_nodes = zeros(length(h.y))
        for i in 1:nf
            tr_nodes .+= (1 - tx) .* real.(@view h.rho[ix0, :, i, i]) .+
                         tx .* real.(@view h.rho[ix1, :, i, i])
        end
        trace_distortion = tr_nodes ./ (nf .* _qke_fd0.(h.y)) .- 1
        g(y) = y^3 * nf * first(_qke_distortion_profile(h.y, trace_distortion, y))
        quadgk(g, 0.0, last(h.y) + 20; rtol=2e-8)[1]
    else
        sum(eachindex(h.y)) do iy
            h.momentum_weights[iy] *
                real(sum(diag(flavor_density_matrix(h, x, h.y[iy]))))
        end
    end
    z = _qke_linear(h.x, h.z, x)
    (11 / 4)^(4 / 3) / z^4 * (120 / (7π^4)) * integral
end

function _qke_linear(grid, values, q)
    i0, i1, t = _qke_bracket(grid, q, "interpolation coordinate")
    (1 - t) * values[i0] + t * values[i1]
end

"""
    qke_neff(history, x=last(history.x); method=:stored)

Return `N_eff` at QKE time `x`.  `method=:stored` reads the independently
computed backend history.  `method=:integral` reintegrates `Tr ϱ` over
comoving momentum using

`N_eff = (11/4)^(4/3) z^(-4) 120/(7π^4) ∫dy y^3 Tr ϱ`.
"""
function qke_neff(h::FlavorDensityHistory, x=last(h.x); method=:stored)
    if method === :stored
        h.neff === nothing && return _qke_neff_integral(h, x)
        return _qke_linear(h.x, h.neff, x)
    elseif method === :integral
        return _qke_neff_integral(h, x)
    end
    throw(ArgumentError("method must be :stored or :integral"))
end

@inline _qke_fd0(y) = y > 300 ? 0.0 : 1 / (exp(y) + 1)

function _qke_distortion_profile(ygrid, distortions, q)
    q < 0 && throw(DomainError(q, "comoving momentum must be non-negative"))
    n = length(ygrid)
    n >= 2 || throw(ArgumentError("QKE momentum interpolation needs at least two nodes"))
    if q <= ygrid[1]
        slope_δ = (distortions[2] - distortions[1]) / (ygrid[2] - ygrid[1])
        δ = distortions[1] + slope_δ * (q - ygrid[1])
    elseif q <= ygrid[end]
        i = min(searchsortedlast(ygrid, q), n - 1)
        slope_δ = (distortions[i + 1] - distortions[i]) /
                  (ygrid[i + 1] - ygrid[i])
        δ = distortions[i] + slope_δ * (q - ygrid[i])
    else
        # The unresolved tail is exponentially suppressed.  Holding the last
        # fractional distortion fixed preserves that exact asymptotic scale
        # without inventing a polynomial continuation from sparse nodes.
        δ = distortions[end]
        slope_δ = 0.0
    end
    fd = _qke_fd0(q)
    raw = fd * (1 + δ)
    value = clamp(raw, 0.0, 1.0)
    raw == value || return value, 0.0
    derivative = -fd * (1 - fd) * (1 + δ) + fd * slope_δ
    value, derivative
end

function _qke_profile(ygrid, values, q)
    # Direct linear interpolation of f on the sparse Gauss--Laguerre nodes
    # overestimates ∫y³f by percent-level amounts.  The collision solution is
    # serialized in the standard form f=f_FD(1+δ); δ is the slowly varying
    # quantity interpolated by NEVO/PRIMAT and by the decoupling literature.
    distortions = values ./ _qke_fd0.(ygrid) .- 1
    _qke_distortion_profile(ygrid, distortions, q)
end

"""
    qke_distribution(history, state; basis=:mass, time=:final,
                     coordinate=:cosmic, N_eff_share=1)

Convert one solved QKE spectrum into Cosmic's general `TabulatedNu`
distribution.  `basis=:mass` is the correct late-time choice for a massive
free-streaming species; `basis=:flavour` exposes a weak-interaction flavour.
`time` may be `:final`, a stored time index, or a value of `x`.

With `coordinate=:cosmic`, the returned callable uses Cosmic's
`q=p/Tν,0` coordinate.  The exact conversion is
`y=q z_final/(11/4)^(1/3) N_eff_share^(1/4)` and must use the same
`N_eff_share` passed to `massive_neutrino`.  `coordinate=:comoving` exposes
the backend's raw `y` coordinate for diagnostics.
"""
function qke_distribution(h::FlavorDensityHistory, state::Integer;
    basis=:mass, time=:final, coordinate=:cosmic, N_eff_share=1.0)

    1 <= state <= flavor_count(h) ||
        throw(BoundsError(1:flavor_count(h), state))
    i0, i1, tx = if time === :final
        (length(h.x), length(h.x), 0.0)
    elseif time isa Integer
        checkbounds(h.x, time)
        (time, time, 0.0)
    else
        _qke_bracket(h.x, time, "x")
    end
    values = if basis === :mass
        h.mass_occupation === nothing &&
            throw(ArgumentError("mass-basis occupations were not saved"))
        (1 - tx) .* collect(@view h.mass_occupation[i0, :, state]) .+
            tx .* collect(@view h.mass_occupation[i1, :, state])
    elseif basis === :flavour
        (1 - tx) .* collect(real.(@view h.rho[i0, :, state, state])) .+
            tx .* collect(real.(@view h.rho[i1, :, state, state]))
    else
        throw(ArgumentError("basis must be :mass or :flavour"))
    end
    distortions = values ./ _qke_fd0.(h.y) .- 1

    momentum_scale = if coordinate === :cosmic
        N_eff_share > 0 ||
            throw(DomainError(N_eff_share, "N_eff_share must be positive"))
        h.z[end] / (11 / 4)^(1 / 3) * N_eff_share^(1 / 4)
    elseif coordinate === :comoving
        1.0
    else
        throw(ArgumentError("coordinate must be :cosmic or :comoving"))
    end

    f0(q) = first(_qke_distortion_profile(h.y, distortions, momentum_scale * q))
    function dlnf0(q)
        y = momentum_scale * q
        v, dv_dy = _qke_distortion_profile(h.y, distortions, y)
        v <= 0 && return -q
        y * dv_dy / v
    end
    TabulatedNu(f0, dlnf0)
end

function _qke_flavor_spectrum(h::FlavorDensityHistory, flavor, Tγ)
    1 <= flavor <= flavor_count(h) ||
        throw(BoundsError(1:flavor_count(h), flavor))
    x = qke_scale_variable(h, Tγ)
    i0, i1, tx = if x <= h.x[1]
        (1, 1, 0.0)
    elseif x >= h.x[end]
        (length(h.x), length(h.x), 0.0)
    else
        _qke_bracket(h.x, x, "x")
    end
    values = (1 - tx) .* collect(real.(@view h.rho[i0, :, flavor, flavor])) .+
             tx .* collect(real.(@view h.rho[i1, :, flavor, flavor]))
    w = (1 - tx) * h.w[i0] + tx * h.w[i1]
    distortions = values ./ _qke_fd0.(h.y) .- 1
    x, w, distortions
end

"""
    weak_rate(Tγ, history::FlavorDensityHistory, τn;
              direction=1, thermal=true, flavor=1)

Evaluate an `n ↔ p` rate with the solved electron-flavour occupation.  The
spectral contribution follows the NEVO/PRIMAT table prescription: replace the
equilibrium neutrino occupation in the Born kernel by `δf`, with the exact
initial-state/Pauli-blocking sign, and multiply it by the same Coulomb and
zero-temperature radiative factors as the base rate.  PRIMAT's table mode
does not add spectral-distortion cross terms to the finite-nucleon-mass or
finite-temperature QED corrections, so neither does this method.
"""
function weak_rate(Tγ, h::FlavorDensityHistory, τn;
    direction=1, thermal=true, flavor=1, ξ_ν=0.0)

    iszero(ξ_ν) || throw(ArgumentError(
        "a zero-asymmetry QKE history cannot be combined with ξ_ν ≠ 0"))
    direction in (-1, 1) ||
        throw(ArgumentError("direction must be +1 (n→p) or -1 (p→n)"))
    xqke, w, distortions = _qke_flavor_spectrum(h, flavor, Tγ)
    Tν = _ME_MEV * w / xqke

    x = _ME_MEV / Tγ
    xν = xqke / w
    s = direction
    pmax = 50 * (1 + max(Tγ, Tν) / _ME_MEV)
    nodes, weights = _gauss_rule(200)
    half = pmax / 2
    acc = 0.0

    χ(en) = begin
        enν = en - s * _q_np
        enν^2 * _fd(enν, xν) * _fd(-en, x)
    end

    function delta_chi(en)
        enν = en - s * _q_np
        y = abs(enν) * xqke
        actual = first(_qke_distortion_profile(h.y, distortions, y))
        equilibrium = _fd(abs(enν), xν)
        signed_delta = enν < 0 ? equilibrium - actual : actual - equilibrium
        enν^2 * signed_delta * _fd(-en, x)
    end

    @inbounds for i in eachindex(nodes)
        p = half * (nodes[i] + 1)
        en = sqrt(p^2 + 1)
        b = p / en
        cp = _weak_correction(b, en, abs(s * _q_np - en), s, 1)
        cm = _weak_correction(b, en, abs(s * _q_np + en), s, -1)
        acc += half * weights[i] * p^2 * (
            (χ(en) + _chi_fm(en, p, x, xν, s, 0.0) + delta_chi(en)) * cp +
            (χ(-en) + _chi_fm(-en, p, x, xν, s, 0.0) + delta_chi(-en)) * cm)
    end
    thermal && (acc += _thermal_lookup(Tγ, s, 0.0))
    acc / (τn * _FN)
end
