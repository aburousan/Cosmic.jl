"""
CMB lensing: the lensing-potential power spectrum and the lensed TT/TE/EE/BB.

Gravitational lensing remaps the CMB by the deflection field ∇ψ, where the
lensing potential is a weighted line-of-sight integral of the Weyl potential,

    ψ(n̂) = -∫ dχ  W(χ) · (φ + ψ)(χ n̂, η₀ - χ)

The default weight is the single-source-distance kernel
f_K(χ*-χ)/(f_K(χ)f_K(χ*)) with χ* at the visibility peak -- the convention
every published C_ℓ^φφ uses, and the correct one at leading order: nearly all
observed anisotropy was imprinted at recombination, so that is where the
deflection acts, even for the few percent of photons that later re-scatter at
reionization (their incident anisotropy was already lensed along the full
path before the re-scattering).

`kernel = :visibility` instead averages the geometric weight over the
visibility function,

    W(χ) = ∫ dχ_s g(χ_s) f_K(χ_s-χ)/(f_K(χ_s)f_K(χ))
         = cot_K(χ)·G₀(χ) - G₁(χ),      G₀ = ∫_{χ_s>χ} g,  G₁ = ∫_{χ_s>χ} g·cot_K(χ_s)

via two cumulative integrals of g. It is kept as a diagnostic, not a default:
it treats re-scattered photons as *new* sources lensed only along their short
path, which drops the anisotropy they inherited -- measured against CLASS it
shifts C_ℓ^φφ by 3-7%, far above the genuine width-of-last-scattering
correction (~0.1%), so it is not the better answer for the lensed two-point
function. Distance-resolved lensing beyond the single-remap picture is beyond
the correlation-function method itself.

The lensed spectra use the full-sky correlation-function method of Challinor &
Lewis (2005), the same algorithm as CLASS's lensing.c: the lensed correlation
functions ξ, ξX, ξ± are assembled at Gauss-Legendre nodes in the separation
angle from the unlensed C_ℓ, C_ℓ^φφ and the deflection correlations σ²(β),
C_gl,2(β), non-perturbatively in σ² and to second order in C_gl,2, then
transformed back with Wigner d-matrices. It is not an expansion in C_ℓ^φφ ·ℓ²
-- the exponential exp(-ℓ(ℓ+1)σ²/4) is kept whole, which is what makes the
damping-tail smoothing and the lensed BB come out right.
"""

using LinearAlgebra: Tridiagonal
using QuadGK: gauss

# --- geometry helpers ---------------------------------------------------------

@inline function _lens_cotK(K, d)
    d <= 0 && return Inf
    if K > 0
        r = sqrt(K)
        return r * cos(r * d) / sin(r * d)
    elseif K < 0
        r = sqrt(-K)
        return r * cosh(r * d) / sinh(r * d)
    end
    1 / d
end

@inline function _sinKd(K, d)
    if K > 0
        r = sqrt(K)
        return sin(r * d) / r
    elseif K < 0
        r = sqrt(-K)
        return sinh(r * d) / r
    end
    d
end

"""
    _lensing_kernel(c, bg, rec; kernel = :lastscatter, nη = 6000, a_min = 1e-5)

The lensing weight W as an interpolant in conformal time. `:lastscatter` is
the standard single-source-distance kernel at χ* = η₀ - η_rec (visibility
peak, Heaviside cut) -- the default and the published-C_ℓ^φφ convention.
`:visibility` averages the source distance over the visibility function; see
the module docstring for why that is a diagnostic, not an improvement.
"""
function _lensing_kernel(c::Cosmology, bg::BackgroundCache, rec::RecombinationSolution;
    kernel::Symbol=:lastscatter, nη=6000, a_min=1e-5)

    K = spatial_curvature_K(c)
    xs = range(log(a_min), 0.0; length=nη)
    ηs = [bg.η(x) for x in xs]
    η0 = bg.η(0.0)

    W = zeros(nη)
    if kernel === :lastscatter
        # z_* from the visibility maximum, as CLASS's τ_rec.
        zg = [10.0^lz for lz in range(2.6, 3.4; length=400)]
        gs = [visibility(rec, z) for z in zg]
        z_star = zg[argmax(gs)]
        η_rec = bg.η(log(1 / (1 + z_star)))
        χ_star = η0 - η_rec
        fχs = _sinKd(K, χ_star)
        for i in 1:nη
            χ = η0 - ηs[i]
            (χ <= 0 || χ >= χ_star) && continue
            W[i] = _sinKd(K, χ_star - χ) / (_sinKd(K, χ) * fχs)
        end
    elseif kernel === :visibility
        g = zeros(nη)
        for (i, x) in enumerate(xs)
            g[i] = bg.κ̇(x) * exp(-optical_depth(rec, redshift(exp(x))))
        end
        # G₀, G₁ cumulate from early times: χ_s > χ  ⟺  η_s < η.
        G0 = zeros(nη)
        G1 = zeros(nη)
        for i in 2:nη
            dη = ηs[i] - ηs[i-1]
            cot1 = _lens_cotK(K, η0 - ηs[i-1])
            cot2 = _lens_cotK(K, η0 - ηs[i])
            G0[i] = G0[i-1] + 0.5 * (g[i-1] + g[i]) * dη
            G1[i] = G1[i-1] + 0.5 * (g[i-1] * cot1 + (isfinite(cot2) ? g[i] * cot2 : 0.0)) * dη
        end
        for i in 1:nη-1
            χ = η0 - ηs[i]
            W[i] = max(_lens_cotK(K, χ) * G0[i] - G1[i], 0.0)
        end
    else
        throw(ArgumentError("kernel must be :visibility or :lastscatter"))
    end

    linear_interpolation(ηs, W; extrapolation_bc=0.0)
end

"""
    Θφ_ℓ_of_k(S, iℓ, B, η0)

Lensing-potential transfer:  Δ^φ_ℓ(k) = ∫ dη S_L(η) j_ℓ(k(η₀-η)), with the
source S_L = -W(χ)·(φ+ψ) built in `source_function`.
"""
function Θφ_ℓ_of_k(S::SourceFunction, iℓ::Int, B::BesselTable, η0)
    isempty(S.S_L) && return 0.0
    acc = 0.0
    η, SL = S.η, S.S_L
    @inbounds for i in 1:(length(η)-1)
        dη = η[i+1] - η[i]
        f1 = SL[i] * B(iℓ, S.k * (η0 - η[i]))
        f2 = SL[i+1] * B(iℓ, S.k * (η0 - η[i+1]))
        acc += 0.5 * (f1 + f2) * dη
    end
    acc
end

# --- Wigner d-functions -------------------------------------------------------

"""
    _wigner_d!(out, m, mp, μ, lmax)

d^ℓ_{m m'}(β) for ℓ = 0..lmax at μ = cos β, by the Kostelec-Rockmore
recurrence on √((2ℓ+1)/2)·d (stable to ℓ of several thousand). Requires
m ≥ |m'| ≥ 0 or m ≥ m' ≥ -m, i.e. m = max -- every pair used in the lensing
assembly is already in that form. Entries below ℓ = m are zero.
"""
function _wigner_d!(out::AbstractVector{Float64}, m::Int, mp::Int, μ::Float64, lmax::Int)
    (m >= abs(mp) && m >= 0) || throw(ArgumentError("need m ≥ |m'|"))
    fill!(out, 0.0)
    lmax >= m || return out

    if m == 0   # Legendre: seed ℓ = 0, 1 explicitly, mm'-term absent
        out[1] = 1.0
        lmax >= 1 || return out
        out[2] = μ
        fl = μ * sqrt(3 / 2)
        flm1 = 1 / sqrt(2)
        @inbounds for l in 1:(lmax-1)
            flp1 = (sqrt((2l + 3) / (2l + 1)) * (2l + 1) * μ * fl -
                    sqrt((2l + 3) / (2l - 1)) * l * flm1) / (l + 1)
            out[l+2] = flp1 * sqrt(2 / (2l + 3))
            flm1 = fl
            fl = flp1
        end
        return out
    end

    # Closed-form seed at ℓ = m:  d^m_{mm'} = √((2m)!/((m+m')!(m-m')!)) c^{m+m'} (-s)^{m-m'}
    ch = sqrt((1 + μ) / 2)               # cos(β/2), β ∈ [0, π]
    sh = sqrt(max(1 - μ, 0.0) / 2)       # sin(β/2)
    lg = 0.5 * (sum(log, 1:2m) - (m + mp > 1 ? sum(log, 1:(m+mp)) : 0.0) -
                (m - mp > 1 ? sum(log, 1:(m-mp)) : 0.0))
    seed = exp(lg) * ch^(m + mp) * (-sh)^(m - mp)
    out[m+1] = seed

    # Scaled recurrence f_ℓ = √((2ℓ+1)/2)·d_ℓ; f_{ℓ0-1} = 0 starts it.
    fl = seed * sqrt((2m + 1) / 2)
    flm1 = 0.0
    mmp = m * mp
    @inbounds for l in m:(lmax-1)
        A = sqrt(((l + 1)^2 - m^2) * ((l + 1)^2 - mp^2)) / (l + 1)
        B = (2l + 1) * (μ - mmp / (l * (l + 1)))
        flp1 = sqrt((2l + 3) / (2l + 1)) * B * fl / A
        if l > m
            C = sqrt((l^2 - m^2) * (l^2 - mp^2)) / l
            flp1 -= sqrt((2l + 3) / (2l - 1)) * C * flm1 / A
        end
        out[l+2] = flp1 * sqrt(2 / (2l + 3))
        flm1 = fl
        fl = flp1
    end
    out
end

# --- lensed spectra by the correlation-function method ------------------------

"""
    lensed_from_cls(TT, EE, TE, φφ; lmax_out, nmu = 0)

Lens dense unlensed spectra (vectors indexed 1 ⟺ ℓ = 0, through
ℓ_unlensed_max) with C_ℓ^φφ, returning (TT, EE, TE, BB) at ℓ = 2..lmax_out.
This is the Challinor-Lewis full-sky correlation-function algorithm; the
unlensed input should extend a few hundred multipoles beyond `lmax_out`.
"""
function lensed_from_cls(TT::Vector{Float64}, EE::Vector{Float64},
    TE::Vector{Float64}, φφ::Vector{Float64}; lmax_out::Int, nmu::Int=0)

    lmaxu = length(TT) - 1
    (length(EE) == length(TE) == length(φφ) == lmaxu + 1) ||
        throw(ArgumentError("unlensed spectra must share one dense ℓ grid"))
    lmax_out <= lmaxu || throw(ArgumentError("lmax_out exceeds the unlensed input"))
    nmu <= 0 && (nmu = lmaxu + 200)

    μs, w8 = gauss(nmu)

    # σ²(β) needs Cgl(1); d^ℓ_11(0) = 1.
    Cgl1 = sum(l -> (2l + 1) * l * (l + 1) * φφ[l+1], 2:lmaxu) / (4π)

    # ℓ-dependent prefactors of the X functions.
    # ℓ = 0, 1 entries are never used (the sums start at ℓ = 2); zero them
    # rather than let the comprehension take sqrt of a negative.
    sqrt1 = [l >= 2 ? sqrt((l + 2) * (l + 1) * l * (l - 1.0)) : 0.0 for l in 0:lmaxu]
    sqrt2 = [l >= 2 ? sqrt((l + 2) * (l - 1.0)) : 0.0 for l in 0:lmaxu]
    sqrt3 = [l >= 2 ? sqrt((l + 3) * (l - 2.0)) : 0.0 for l in 0:lmaxu]
    sqrt4 = [l >= 3 ? sqrt((l + 4) * (l + 3) * (l - 2.0) * (l - 3)) : 0.0 for l in 0:lmaxu]
    sqrt5 = [sqrt(l * (l + 1.0)) for l in 0:lmaxu]

    # Chunked tasks with task-local buffers: safe under task migration, and the
    # twelve d-vectors are allocated once per chunk, not once per node.
    nchunk = max(Threads.nthreads(), 1)
    ranges = [i:nchunk:nmu for i in 1:nchunk]
    tasks = map(ranges) do idxs
        Threads.@spawn begin
            aTT = zeros(lmax_out + 1)
            aTE = zeros(lmax_out + 1)
            aP = zeros(lmax_out + 1)     # ξ₊ → (EE+BB)
            aM = zeros(lmax_out + 1)     # ξ₋ → (EE-BB)
            d = zeros(lmaxu + 1, 12)
            for imu in idxs
                μ = μs[imu]
                w = w8[imu]
                for (j, (m, mp)) in enumerate(_D_PAIRS)
                    _wigner_d!(view(d, :, j), m, mp, μ, lmaxu)
                end
                d00 = view(d, :, 1); d11 = view(d, :, 2); d1m1 = view(d, :, 3)
                d2m2 = view(d, :, 4); d20 = view(d, :, 5); d3m1 = view(d, :, 6)
                d4m2 = view(d, :, 7); d22 = view(d, :, 8); d31 = view(d, :, 9)
                d3m3 = view(d, :, 10); d40 = view(d, :, 11); d4m4 = view(d, :, 12)

                Cgl = 0.0
                Cgl2 = 0.0
                @inbounds for l in 2:lmaxu
                    f = (2l + 1) * l * (l + 1) * φφ[l+1]
                    Cgl += f * d11[l+1]
                    Cgl2 += f * d1m1[l+1]
                end
                Cgl /= 4π
                Cgl2 /= 4π
                σ2 = Cgl1 - Cgl

                # Lensed correlation functions at this angle, with the unlensed
                # part subtracted term by term (added back exactly at the end)
                # so the quadrature integrates only the small lensing change.
                ξ = 0.0
                ξX = 0.0
                ξp = 0.0
                ξm = 0.0
                @inbounds for l in 2:lmaxu
                    ll = Float64(l)
                    fac = ll * (ll + 1) / 4
                    fac1 = (2ll + 1) / (4π)

                    X000 = exp(-fac * σ2)
                    Xp000 = -fac * X000
                    X220 = 0.25 * sqrt1[l+1] * X000
                    X022 = X000 * (1 + σ2 * (1 + 0.5 * σ2))
                    Xp022 = -(fac - 1) * X022
                    X242 = 0.25 * sqrt4[l+1] * X000
                    X121 = -0.5 * sqrt2[l+1] * X000 * (1 + 2σ2 / 3)
                    X132 = -0.5 * sqrt3[l+1] * X000 * (1 + 5σ2 / 3)

                    ξ += fac1 * TT[l+1] * (X000^2 * d00[l+1] +
                         Xp000^2 * d1m1[l+1] * Cgl2 * 8 / (ll * (ll + 1)) +
                         (Xp000^2 * d00[l+1] + X220^2 * d2m2[l+1]) * Cgl2^2 -
                         d00[l+1])

                    ξX += fac1 * TE[l+1] * (X022 * X000 * d20[l+1] +
                          Cgl2 * 2 * Xp000 / sqrt5[l+1] *
                          (X121 * d11[l+1] + X132 * d3m1[l+1]) +
                          0.5 * Cgl2^2 *
                          ((2 * Xp022 * Xp000 + X220^2) * d20[l+1] +
                           X220 * X242 * d4m2[l+1]) -
                          d20[l+1])

                    # Unlensed scalar BB is zero, so ξ± both carry just EE.
                    ξp += fac1 * EE[l+1] * (X022^2 * d22[l+1] +
                          2 * Cgl2 * X132 * X121 * d31[l+1] +
                          Cgl2^2 * (Xp022^2 * d22[l+1] + X242 * X220 * d40[l+1]) -
                          d22[l+1])

                    ξm += fac1 * EE[l+1] * (X022^2 * d2m2[l+1] +
                          Cgl2 * (X121^2 * d1m1[l+1] + X132^2 * d3m3[l+1]) +
                          0.5 * Cgl2^2 * (2 * Xp022^2 * d2m2[l+1] +
                           X220^2 * d00[l+1] + X242^2 * d4m4[l+1]) -
                          d2m2[l+1])
                end

                @inbounds for l in 2:lmax_out
                    aTT[l+1] += w * ξ * d00[l+1]
                    aTE[l+1] += w * ξX * d20[l+1]
                    aP[l+1] += w * ξp * d22[l+1]
                    aM[l+1] += w * ξm * d2m2[l+1]
                end
            end
            (aTT, aTE, aP, aM)
        end
    end
    parts = fetch.(tasks)

    TTl = zeros(lmax_out + 1)
    EEl = zeros(lmax_out + 1)
    TEl = zeros(lmax_out + 1)
    BBl = zeros(lmax_out + 1)
    for l in 2:lmax_out
        dTT = 2π * sum(p -> p[1][l+1], parts)
        dTE = 2π * sum(p -> p[2][l+1], parts)
        cp = π * sum(p -> p[3][l+1], parts)
        cm = π * sum(p -> p[4][l+1], parts)
        TTl[l+1] = TT[l+1] + dTT
        TEl[l+1] = TE[l+1] + dTE
        EEl[l+1] = EE[l+1] + (cp + cm)
        BBl[l+1] = cp - cm
    end
    (TT=TTl, EE=EEl, TE=TEl, BB=BBl)
end

const _D_PAIRS = ((0, 0), (1, 1), (1, -1), (2, -2), (2, 0), (3, -1),
    (4, -2), (2, 2), (3, 1), (3, -3), (4, 0), (4, -4))

"""
    LensedCMBSpectra

Lensed spectra at every integer ℓ from 2 to `lmax`, in μK². `BB` is the
lensing B-mode generated from E; the unlensed scalar BB is zero.
"""
struct LensedCMBSpectra
    ℓ::Vector{Int}
    TT::Vector{Float64}
    EE::Vector{Float64}
    TE::Vector{Float64}
    BB::Vector{Float64}
end

function D_ℓ(s::LensedCMBSpectra, which::Symbol=:TT)
    C = which === :TT ? s.TT : which === :EE ? s.EE :
        which === :TE ? s.TE : which === :BB ? s.BB :
        throw(ArgumentError("which must be :TT, :EE, :TE or :BB"))
    [ℓ * (ℓ + 1) * C[i] / (2π) for (i, ℓ) in enumerate(s.ℓ)]
end

"""
Natural cubic spline on an irregular grid -- the ℓ-sampling of C_ℓ is
irregular by design, and Interpolations.jl's cubic splines need uniform
spacing. Standard tridiagonal construction; second derivatives vanish at the
ends, as in CLASS's spline of the same quantity.
"""
function _natural_spline(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    h = diff(x)
    dl = h[1:end-1]
    du = h[2:end]
    dg = 2 .* (h[1:end-1] .+ h[2:end])
    rhs = [6 * ((y[i+2] - y[i+1]) / h[i+1] - (y[i+1] - y[i]) / h[i]) for i in 1:n-2]
    M = zeros(n)
    n > 2 && (M[2:n-1] = Tridiagonal(dl[2:end], dg, du[1:end-1]) \ rhs)
    function (q)
        j = clamp(searchsortedlast(x, q), 1, n - 1)
        t = q - x[j]
        hj = h[j]
        a = (M[j+1] - M[j]) / (6hj)
        b = M[j] / 2
        c = (y[j+1] - y[j]) / hj - hj * (2M[j] + M[j+1]) / 6
        y[j] + t * (c + t * (b + t * a))
    end
end

"Spline the sparsely sampled C_ℓ onto every integer ℓ (via ℓ(ℓ+1)C_ℓ)."
function _dense_cls(ℓs::Vector{Int}, C::Vector{Float64}, lmaxu::Int)
    y = [ℓ * (ℓ + 1) * C[i] for (i, ℓ) in enumerate(ℓs)]
    itp = _natural_spline(Float64.(ℓs), y)
    out = zeros(lmaxu + 1)
    for l in 2:lmaxu
        out[l+1] = itp(Float64(l)) / (l * (l + 1))
    end
    out
end

"""
    lensed_cmb_spectra(c, rec; lmax = 2000, dl_buffer = 500, kernel = :lastscatter, ...)

Unlensed spectra plus C_ℓ^φφ in one mode sweep (`cmb_spectra(...; lensing =
true)`), then the correlation-function lensing. Returns `(lensed,
unlensed)` -- the `LensedCMBSpectra` and the `CMBSpectra` carrying φφ and Tφ.

`dl_buffer` sets how far past `lmax` the unlensed spectra are computed, because
lensing couples each multipole to its neighbours. TT/TE/EE converge with the
default 500 across the whole range; lensed BB is a pure off-diagonal EE→BB
transfer that draws power from a broad band of ℓ', so its top ~`dl_buffer/2`
multipoles are under-converged -- with buffer 500, BB is trustworthy to about
ℓ = lmax − 300 and falls ~20% low at ℓ = lmax itself. Ask for BB out to some ℓ
by giving `lmax` a few hundred above it, or raise `dl_buffer`. Reaching φφ (and
thus lensed BB) to high ℓ also needs a large `kmax` in the sweep: the lensing
kernel maps ℓ to k ≈ ℓ/χ over a range of χ down to small values, so C_ℓ^φφ at
ℓ = 2000 wants `kmax ≈ 0.7`, not the `2.2·lmax/η₀` that suffices for TT.
"""
function lensed_cmb_spectra(c::Cosmology, rec::RecombinationSolution;
    lmax=2000, dl_buffer=500, kernel::Symbol=:lastscatter, nmu=0, kwargs...)

    lmaxu = lmax + dl_buffer
    spec = cmb_spectra(c, rec; lmax=lmaxu, lensing=true, lensing_kernel=kernel, kwargs...)
    TT = _dense_cls(spec.ℓ, spec.TT, lmaxu)
    EE = _dense_cls(spec.ℓ, spec.EE, lmaxu)
    TE = _dense_cls(spec.ℓ, spec.TE, lmaxu)
    φφ = _dense_cls(spec.ℓ, spec.φφ, lmaxu)
    lens = lensed_from_cls(TT, EE, TE, φφ; lmax_out=lmax, nmu)
    ℓs = collect(2:lmax)
    (lensed=LensedCMBSpectra(ℓs, lens.TT[3:end], lens.EE[3:end],
            lens.TE[3:end], lens.BB[3:end]),
        unlensed=spec)
end
