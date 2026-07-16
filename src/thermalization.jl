"""
Photon thermalization solver: derive μ and y by evolving the CMB spectrum itself.

Where `distortions.jl` convolves the heating history with Chluba's precomputed
thermalization Green's function, this module *computes that response from scratch*
-- it integrates the photon Boltzmann equation forward in redshift and reads the
distortion off the evolved spectrum. No Green's-function file, no branching-ratio
table: the μ/y split emerges from the competition between Compton scattering and
photon production, exactly as in the real universe.

Governing equation for the occupation number n(x, z), x = hν/(k T_ref),
T_ref = T₀(1+z) (see `docs/src/spectral_distortions.md`):

    dn/dt = 𝒞_KN[n] + ℰ[n]
    𝒞_KN[n] = nₑσ_Tc ∫dω' P_KN(ω→ω'; T_e)
                  [n(ω')(1+n(ω)) - n(ω)(1+n(ω'))]               Compton
    ℰ[n] = (Λ_DC + Λ_BR) e^{-x}/x³ [ 1 - n(e^x - 1) ]                double Compton + bremsstrahlung
    θ_e' = σ_T nₑ c · kT_ref/(m_e c²)

with dt = -dz/[(1+z) H]. `P_KN` is the exact Klein–Nishina kernel folded with
the relativistic Maxwell–Jüttner electron distribution. It is generated with
Double64 quadrature from the FORM-checked tree amplitude; no Fokker–Planck or
Kompaneets approximation is used in the production solve.
"""

using OrdinaryDiffEq: ODEProblem, ODEFunction, solve, Rodas5P, ReturnCode
using LinearAlgebra: dot
using Interpolations: linear_interpolation
using QuadGK: quadgk

const _Cst = Constants

"""
    ThermalizationSolver

Log-spaced frequency grid and the cosmology/recombination it will run against.
Build with [`thermalization_grid`](@ref).
"""
struct ThermalizationSolver{C,R}
    cosmo::C
    rec::R
    x::Vector{Float64}          # dimensionless frequency grid
    lnx::Vector{Float64}
end

"""
    thermalization_grid(c, rec; xmin = 1e-5, xmax = 50, nx = 320)

The frequency grid.

`xmin` must sit **below** the critical frequency x_c = √(Λ/κ) at which photon
production balances Compton scattering. x_c ≈ 4e-3 - 9e-3 across the μ era, so the
default sits an order of magnitude under it. Below x_c the emission term physically
pins the spectrum to a blackbody; if the grid *ends* near x_c instead, the boundary
does the pinning by fiat and acts as a huge artificial photon sink that
over-thermalizes μ toward zero.

Pushing xmin far lower does not help: the emission rate grows as 1/x², so xmin =
1e-4 makes the system stiff enough (~1e7 × H) to break the implicit solver on the
initial transient, while the answer is already converged by xmin = 1e-3.
"""
function thermalization_grid(c::Cosmology, rec::RecombinationSolution;
    xmin=5e-4, xmax=50.0, nx=384)
    nx >= 64 || throw(ArgumentError("exact Compton redistribution requires nx ≥ 64"))
    lnx = collect(range(log(xmin), log(xmax); length=nx))
    ThermalizationSolver(c, rec, exp.(lnx), lnx)
end

# --- Rate coefficients ------------------------------------------------------

"""
    _theta_gamma(c, z)

kT_γ/(m_e c²), the dimensionless photon temperature. Sets the Compton energy-
exchange efficiency.
"""
_theta_gamma(c::Cosmology, z) = _Cst.k_B_SI * c.Tcmb * (1 + z) / (_Cst.m_e_SI * _Cst.c_SI^2)

const _α_fs = 1 / 137.035999
const _I4pl = 4π^4 / 15      # ∫x⁴ n_Pl(n_Pl+1) dx = 25.976

# --- Exact free-free Gaunt factor --------------------------------------------
#
# ⟨g_ff⟩(γ², u): thermally averaged Karzas-Latter Gaunt factor, from the table
# computed offline (exact ₂F₁ evaluation, Maxwell average, 30-digit
# working precision, verified against van Hoof et al 2014). Bilinear in the log-log
# grid; outside the tabulated window the value is clamped to the edge, which is
# safe because the window was chosen to cover the solver's whole (z, x) range.

const _GAUNT_TBL = Ref{Union{Nothing,NamedTuple}}(nothing)

function _load_gaunt()
    if _GAUNT_TBL[] === nothing
        M = readdlm(joinpath(@__DIR__, "..", "data", "thermalization", "gaunt_ff.dat"))
        lg2 = Float64.(M[1, 2:end])
        lu = Float64.(M[2:end, 1])
        _GAUNT_TBL[] = (lg2=lg2, lu=lu, g=Float64.(M[2:end, 2:end]),
            dg2=lg2[2] - lg2[1], du=lu[2] - lu[1])
    end
    _GAUNT_TBL[]
end

function gaunt_ff(γ2, u)
    t = _load_gaunt()
    a = clamp((log10(γ2) - t.lg2[1]) / t.dg2, 0.0, length(t.lg2) - 1.000001)
    b = clamp((log10(u) - t.lu[1]) / t.du, 0.0, length(t.lu) - 1.000001)
    i, fi = floor(Int, a) + 1, a - floor(a)
    j, fj = floor(Int, b) + 1, b - floor(b)
    (1 - fi) * (1 - fj) * t.g[j, i] + fi * (1 - fj) * t.g[j, i+1] +
    (1 - fi) * fj * t.g[j+1, i] + fi * fj * t.g[j+1, i+1]
end

"""
    H_dc(x)

Effective double-Compton correction factor (Chluba & Sunyaev 2012, eq. 13):

    H_dc(x) ≈ e^{-2x} [ 1 + (3/2)x + (29/24)x² + (11/16)x³ + (5/12)x⁴ ]

H_dc(0) = 1; it carries the entire frequency dependence of the DC emissivity
beyond the soft-photon limit, and falls as x⁴e^{-2x} at large x.

Kept as the published reference point; the solver itself uses [`g_dc_exact`](@ref),
which replaces this fit (and its 1/(1+14.16θ) companion) with the full
tree-level QED result.
"""
H_dc(x) = exp(-2x) * (1 + 1.5x + (29 / 24) * x^2 + (11 / 16) * x^3 + (5 / 12) * x^4)

# g_dc(x, θ): exact double-Compton emission factor. Computed from the
# polarization-summed 2→3 matrix element (all six diagrams, FORM traces --
# FORM (exact trace algebra) and an offline phase-space integrator), integrated over the Planck seed phase space
# with the pair counted at its softer member (ω₁ ≥ ω₂; the partner-soft region
# is the tree-level infrared divergence that cancels against the virtual
# correction to elastic Compton and so belongs to the Compton kernel, not
# here). Evaluated in double-double arithmetic -- Float64 loses the trace
# cancellations below θ ~ 1e-3 -- and normalized to the analytic soft limit,
# lim_{x→0} g_dc = I₄^pl, which the table's corner reproduces to 4e-5.
# Generation: an offline emissivity generator; certification: positivity + eikonal
# closure of |M|² (dc_check.jl), exact θ² rate scaling, node-doubling
# convergence ≤1.1e-3, and agreement with the CS2012 fit to ≲0.7% for x ≤ 1
# where that fit is good. Bilinear in (log10 θ, log10 x), edge-clamped: below
# x = 1e-4 the function is its plateau, above x = 20 it is ~e^{-2x} ≈ 0 and
# bremsstrahlung owns the spectrum anyway.

const _DCGDC_TBL = Ref{Union{Nothing,NamedTuple}}(nothing)

function _load_dcgdc()
    if _DCGDC_TBL[] === nothing
        M = readdlm(joinpath(@__DIR__, "..", "data", "thermalization", "dc_gdc.dat"))
        xg = Float64.(M[1, 2:end])
        lx = log10.(xg)
        lθ = log10.(Float64.(M[2:end, 1]))
        # Interpolate ln[g·e^{2x}], not g. The table is exact at its nodes, but
        # g ~ e^{-2x} falls ~30x per 0.075-dex cell at x ~ 10, and a linear
        # chord there overshoots mid-cell by factors of a few -- enough spurious
        # K_DC in the FIRAS band (absorption is amplified by e^x) to eat ~20%
        # of mu. Peeling off the exact e^{-2x} leaves a slowly varying,
        # polynomial-like surface that bilinear interpolation carries at
        # sub-percent accuracy; the exponential is reattached exactly at the
        # query point.
        lh = log.(Float64.(M[2:end, 2:end])) .+ 2 .* xg'
        _DCGDC_TBL[] = (lx=lx, lθ=lθ, lh=lh,
            dx=lx[2] - lx[1], dθ=lθ[2] - lθ[1])
    end
    _DCGDC_TBL[]
end

"""
    g_dc_exact(x, θ)

Exact double-Compton emission factor g_dc(x, θ) (see the table notes above);
g_dc_exact(x→0, θ→0) → I₄^pl = 4π⁴/15.
"""
function g_dc_exact(x, θ)
    t = _load_dcgdc()
    a = clamp((log10(x) - t.lx[1]) / t.dx, 0.0, length(t.lx) - 1.000001)
    b = clamp((log10(θ) - t.lθ[1]) / t.dθ, 0.0, length(t.lθ) - 1.000001)
    i, fi = floor(Int, a) + 1, a - floor(a)
    j, fj = floor(Int, b) + 1, b - floor(b)
    lh = (1 - fi) * (1 - fj) * t.lh[j, i] + fi * (1 - fj) * t.lh[j, i+1] +
         (1 - fi) * fj * t.lh[j+1, i] + fi * fj * t.lh[j+1, i+1]
    # The surface index is clamped to the grid, but the exponential uses the
    # true x: past the table edge g keeps falling as e^{-2x} with the edge
    # polynomial, instead of freezing at the edge value.
    exp(lh - 2x)
end

"""
    K_emission(c, rec, z, x)

Photon-production coefficient K(x) = K_DC + K_BR (Chluba & Sunyaev 2012, §2.2.2),
which enters the photon Boltzmann equation as

    dn/dτ |_{DC+BR} = (1/x³) [ 1 - n(e^{x_e} - 1) ] · K(x)

with τ the Compton optical depth (dτ = σ_T nₑ c dt). Dimensionless.

**Double Compton** (their eqs. 10, 11 give the structure; the factor itself is exact):

    K_DC = (4α/3π) θ² · g_dc_exact(x, θ)

where g_dc_exact is the tabulated exact tree-level QED result (see above) whose
soft limit is I₄^pl = ∫x⁴ n_Pl(n_Pl+1)dx = 4π⁴/15 ≈ 25.976; CS2012's
[I₄^pl/(1+14.16θ)]·H_dc(x) is its fitted stand-in, kept only as a reference.

That integral is the *seed-photon phase space*: double Compton emission is
stimulated by the ambient blackbody, so the rate carries a factor of ~26. Reading
K_DC = (4α/3π)θ² and setting g_dc = 1 -- which is what the headline formula
invites -- underestimates photon production 26-fold and throws the thermalization
redshift badly off. The 14.16θ term is the leading relativistic correction and
only matters for z ≳ few×10⁶.

**Bremsstrahlung** (their eq. 14), fully ionized H+He:

    K_BR = [ α λ_e³ θ^{-7/2} e^{-x} / (2π √(6π)) ] · Σᵢ Zᵢ² Nᵢ · g_ff

with λ_e = h/(m_e c) the *unreduced* Compton wavelength. DC dominates for
z ≳ few×10⁵; BR takes over below that and is what keeps the low-frequency
spectrum a blackbody down to recombination.

The ⁷Li/⁷Be trace ions are not folded into ΣZ²N: even taking every such
nucleus at its largest ionic charge gives
Σ_Li,Be Z²N/(Σ_H,He Z²N) < 5×10⁻⁹ for the network abundances. Their
capture *energy* is retained by `nuclear_deposition_rate`; their bounded
free–free opacity is numerically irrelevant at the stated precision.
"""
function K_emission(c::Cosmology, rec::RecombinationSolution, z, x)
    θ = _theta_gamma(c, z)

    K_DC = (4 * _α_fs / (3π)) * θ^2 * g_dc_exact(x, θ)

    λe = 2π * _Cst.ħ_SI / (_Cst.m_e_SI * _Cst.c_SI)     # h/(m_e c), not ħ/(m_e c)
    # Exact thermally averaged free-free Gaunt factor: Karzas & Latter's
    # Coulomb-wave result (₂F₁ closed form), Maxwell-averaged, computed in
    # an offline generator at 30-digit precision and cached in
    # data/thermalization/gaunt_ff.dat -- verified point by point against van
    # Hoof et al 2014's published table. It is a function of γ² = Z²Ry/kT_e and
    # u = hν/kT_e only, so hydrogen and helium read the same table at γ² and
    # 4γ². Setting g_ff = 1 (or a Born-limit log) misprices photon production
    # by factors of a few exactly where bremsstrahlung pins the spectrum to a
    # blackbody, and μ feels it through the critical frequency x_c ∝ √K_BR.
    Tz = c.Tcmb * (1 + z)
    γ2 = 157887.512 / Tz                              # Ry/k_B in K, over T_e ≈ T_γ
    nH = n_H_of_z(c, z)
    ΣZ2Ng = nH * (gaunt_ff(γ2, x) + 4 * f_He(c) * gaunt_ff(4γ2, x))
    K_BR = (_α_fs * λe^3 * θ^(-7 / 2) * exp(-x) / (2π * sqrt(6π))) * ΣZ2Ng

    K_DC + K_BR
end

"Compton optical-depth rate dτ/dt = σ_T nₑ c, in 1/s."
_dtau_dt(c::Cosmology, rec::RecombinationSolution, z) =
    _Cst.σ_T_SI * n_e(rec, z) * _Cst.c_SI

# --- Reference branching ratios (Chluba 2013) --------------------------------

"""
    chluba_J_y(z), chluba_J_mu(z), chluba_J_vis(z)

Chluba (2013), MNRAS 434, 352, eqs. (5)-(6): fits to the *exact* thermalization
Green's function, and the reference this module is validated against.

    J_y(z)  = [ 1 + ((1+z)/6.0e4)^2.58 ]⁻¹
    J_μ(z)  = 1 - exp[ -((1+z)/5.8e4)^1.88 ]
    𝒥(z)    = exp[ -(z/z_μ)^{5/2} ]            distortion visibility, z_μ ≈ 1.98e6

    G_th(ν,z) ≈ 1.401 J_μ(z) 𝒥(z) M(ν) + (J_y/4) Y_SZ(ν) + ((1-𝒥)/4) G(ν)

so the coefficient of M -- which is exactly what [`distortion_from_spectrum`]
returns as `μ` -- is 1.401·J_μ·𝒥, and the coefficient of Y_SZ is J_y/4. The two
cross at z ≈ 5.3e4, which is what makes 5e4 the conventional μ/y boundary.

**Why not validate against CLASS's `FIRAS_branching_ratios.dat`?** Because those
are not the physical branching ratios. They are a PCA decomposition in a basis
built from FIRAS's 43 frequency channels, and Chluba & Jeong (2014) say outright
that "these ratios are *not unique but depend on the experimental settings*". Their
J_μ reaches 1.896 at z = 1e5, which is above the photon-number ceiling of 1.401 and
is a projection coefficient, not an energy branching. Chluba's own physical value
there is 1.312; this solver gets 1.354.
"""
chluba_J_y(z) = 1 / (1 + ((1 + z) / 6.0e4)^2.58)
chluba_J_mu(z) = 1 - exp(-((1 + z) / 5.8e4)^1.88)
chluba_J_vis(z; z_μ=1.98e6) = exp(-(z / z_μ)^2.5)

# --- Photon-production shapes and forward solve ----------------------------

struct _EmissionOp
    src_sh::Vector{Float64}
    snk_sh::Vector{Float64}
end

function _EmissionOp(x::Vector{Float64})
    # Linearised photon production, Chluba & Sunyaev (2012) eq. 9. With
    # Δx_e = x_e - x ≈ -xφ and x_e = x T_z/T_e,
    #
    #   dΔn/dτ|_em = K(x)/x³ · [ xφ/(1 - e^{-x})  -  Δn (eˣ - 1) ]
    #                          └── source ──┘      └─── sink ───┘
    #
    # The source is the term that deposits thermalized energy into a temperature
    # shift; dropping it destroys that energy instead. Both pieces carry the same
    # K(x), so only these two x-shapes are needed:
    src_sh = @. 1 / (x^2 * (1 - exp(-x)))     # K·src_sh · φ   (from xφ/(1-e^{-x})/x³)
    snk_sh = @. (exp(x) - 1) / x^3            # K·snk_sh · Δn  (sink)
    _EmissionOp(src_sh, snk_sh)
end

"""
    solve_thermalization(grid; z_start = 5e6, z_end = 200, sources...)

Evolve the photon distortion from `z_start` (blackbody; thermalization complete)
down to `z_end`, injecting the guaranteed heating, and return the solution.

Stiff (the Compton rate dwarfs H at high z), so an implicit solver is used. The
Jacobian is tridiagonal plus a rank-1 term (the electron temperature is an
integral over the whole spectrum), so it is taken dense.
"""
function solve_thermalization(grid::ThermalizationSolver;
    z_start=2e6, z_end=500.0, reltol=1e-8, abstol=1e-18, dtmax=0.02,
    source_nodes=600, sources=(;))

    c, rec, x = grid.cosmo, grid.rec, grid.x
    nx = length(x)
    source_nodes >= 2 || throw(ArgumentError("source_nodes must be at least 2"))
    op = _EmissionOp(x)
    cop = ExactComptonOperator(x)
    Ccomp = zeros(nx)
    emsrc = zeros(nx)
    emsnk = zeros(nx)
    eweight = cop.eweight
    penergy = zeros(nx)
    # Solve the trace-nuclear populations once and carry that same history
    # through every RHS evaluation. This is both faster and guarantees that
    # thermalization sees the composition produced by this cosmology's BBN.
    source_opts = if hasproperty(sources, :nuclear_history)
        sources
    elseif get(sources, :nuclear, true)
        merge(sources, (nuclear_history=trace_nuclear_history(c, rec),))
    else
        sources
    end

    # The guaranteed heating terms are functions of redshift alone, but each
    # contains adaptive damping-scale integrals. Evaluating them inside every
    # Newton/GMRES residual repeats the same quadratures thousands of times.
    # Cache their exact evaluations on the integration coordinate; 600/900-node
    # closure is a convergence gate in an offline closure-test script.
    lz_source = collect(range(log(z_end + 1), log(z_start + 1); length=source_nodes))
    q_source = [heating_rate(c, rec, exp(lz) - 1; source_opts...) for lz in lz_source]
    q_of_lz = linear_interpolation(lz_source, q_source)

    # Integration variable lz = ln(1+z); dt/dlz = -1/H.
    function rhs!(dΔn, Δn, p, lz)
        z = exp(lz) - 1
        H = 100 * c.h * 1e3 / _Cst.Mpc_SI * E_z(c, z)
        dτdt = _dtau_dt(c, rec, z)

        # Feed the physical energy injection to the exact collision operator. It
        # solves for T_e by imposing energy balance on the full redistribution
        # integral, rather than using the Kompaneets moment expansion.
        θ = _theta_gamma(c, z)
        q_dz = q_of_lz(lz)
        qdot_frac = q_dz * H * (1 + z)
        heat_per_tau = dτdt > 0 ? qdot_frac / dτdt : 0.0

        @inbounds for i in 1:nx
            K = K_emission(c, rec, z, x[i])
            emsrc[i] = K * op.src_sh[i]
            emsnk[i] = K * op.snk_sh[i]
        end
        φC = exact_compton_collision!(Ccomp, Δn, cop, θ, heat_per_tau;
            emission_source=emsrc, emission_sink=emsnk)

        @inbounds for i in 1:nx
            em = dτdt * (φC * emsrc[i] - Δn[i] * emsnk[i])
            dΔn[i] = (dτdt * Ccomp[i] + em) * (-1 / H)
        end
        nothing
    end

    # The equation is affine-linear in Δn at the stated first-order distortion
    # level. Assemble its exact Jacobian from the same detailed-balance weak
    # form instead of finite-differencing the expensive redistribution integral
    # once per frequency bin. This changes neither the collision physics nor the
    # electron-energy constraint; it only exposes the derivative analytically to
    # the implicit solver.
    function jac!(J, Δn, p, lz)
        z = exp(lz) - 1
        H = 100 * c.h * 1e3 / _Cst.Mpc_SI * E_z(c, z)
        dτdt = _dtau_dt(c, rec, z)
        θ = _theta_gamma(c, z)

        @inbounds for i in 1:nx
            K = K_emission(c, rec, z, x[i])
            emsrc[i] = K * op.src_sh[i]
            emsnk[i] = K * op.snk_sh[i]
        end

        _exact_collision_emission_jacobian!(J, penergy, cop, θ, emsrc, emsnk)
        scale = -dτdt / H
        @inbounds for j in 1:nx, i in 1:nx
            J[i, j] *= scale
        end
        nothing
    end

    f = ODEFunction(rhs!; jac=jac!)
    prob = ODEProblem(f, zeros(nx), (log(z_start + 1), log(z_end + 1)))
    # z_start sits at the thermalization wall (z_th ≈ 2e6): anything injected above
    # it is erased into a blackbody anyway, and starting higher only buys crippling
    # stiffness. z_end = 500 is past recombination, where Silk damping has stopped
    # and the distortion is frozen.
    # The exact redistribution operator is nonlocal in frequency and extremely
    # stiff at the thermalization wall. A concrete finite-difference Jacobian is
    # more expensive than matrix-free GMRES, but the latter loses the initial
    # solve at z=4e6. Source quadratures are cached above, so the robust direct
    # implicit solve remains practical.
    solve(prob, Rodas5P(autodiff=false); reltol, abstol, maxiters=10^8, dtmax)
end

# --- Reading μ and y off the spectrum ---------------------------------------

"""
    distortion_from_spectrum(grid, n_final)

Project the evolved distortion Δn = n_final - n_Planck onto the y and μ spectral
shapes to extract the amplitudes, using the photon-number and energy weights that
define them:

    y = ⟨Δn, Y⟩ / ⟨Y, Y⟩,   μ = ⟨Δn, M⟩ / ⟨M, M⟩

with the same G/M/Y basis as the Green's-function method. Returns `(μ, y, Δρ_ρ)`.
"""
function distortion_from_spectrum(grid::ThermalizationSolver, Δn_final)
    x = grid.x
    # The solver evolves the distortion Δn directly. In *intensity*, ΔI ∝ x³ Δn --
    # this is what the G/Y/M shapes below carry (the extra x³). Projecting the
    # occupation-number Δn against intensity shapes would be a units error.
    ΔI = @. x^3 * Δn_final

    # Intensity spectral shapes (same as CLASS's Gdist/Ydist/Mdist).
    G = @. x^4 * exp(x) / (exp(x) - 1)^2                 # temperature-shift shape
    Y = @. G * (x / tanh(x / 2) - 4)                    # y-distortion shape
    M = @. G * (1 / 2.19229 - 1 / x)                    # μ-distortion shape

    # Gram-Schmidt orthonormalise {Y, M, G} and project ΔI, exactly as the
    # Green's-function method projects G_th (Chluba & Jeong 2014, App. A). The
    # inner product is the plain grid dot product CLASS uses.
    dot_(a, b) = sum(a .* b)
    nrm(v) = sqrt(dot_(v, v))

    e_Y = Y ./ nrm(Y)
    M_Y = dot_(M, e_Y); G_Y = dot_(G, e_Y)
    Mperp = M .- M_Y .* e_Y; e_M = Mperp ./ nrm(Mperp); G_M = dot_(G, e_M)
    Gperp = G .- G_Y .* e_Y .- G_M .* e_M

    f_g = dot_(ΔI, Gperp ./ nrm(Gperp)) / nrm(Gperp)
    f_mu = (dot_(ΔI, e_M) - G_M * f_g) / nrm(Mperp)
    f_y = (dot_(ΔI, e_Y) - M_Y * f_mu - G_Y * f_g) / nrm(Y)

    # Fractional energy in the distortion, Δρ/ρ = ∫x³Δn dx / ∫x³n_pl dx. Energy
    # conservation ties this to the total injected ΔQ/ρ_γ -- the cleanest check
    # on the evolution, independent of the μ/y projection.
    Δρ_ρ = _trap_logx(x, ΔI) / _trap_logx(x, @. x^3 / (exp(x) - 1))

    # Here we project the *actual* distortion intensity ΔI onto the shapes, so the
    # amplitudes come out directly: ΔI = y·Y + μ·M + g·G means f_y = y, etc. (This
    # differs from CLASS's generate_PCA_files, which projects the per-unit-energy
    # Green's function G_th and so carries the extra J_y = 4f_y, J_μ = f_μ/1.401
    # energy-to-amplitude conversions.) Checked: for pure y-era injection this
    # gives y = (1/4)Δρ/ρ, the textbook relation.
    (μ=f_mu, y=f_y, g=f_g, Δρ_ρ=Δρ_ρ)
end

"Trapezoid of f over a log-spaced x grid: ∫ f dx = ∫ f·x dlnx."
function _trap_logx(x, f)
    s = 0.0
    @inbounds for i in 1:(length(x)-1)
        dlnx = log(x[i+1]) - log(x[i])
        s += 0.5 * (f[i] * x[i] + f[i+1] * x[i+1]) * dlnx
    end
    s
end

"""
    thermalization_distortions(c, rec; kws...)

Convenience: build the grid, solve the thermalization equation, and return the
distortion amplitudes `(μ, y, Δρ_ρ)` derived from the evolved spectrum.
"""
function thermalization_distortions(c::Cosmology, rec::RecombinationSolution;
    nx=512, z_start=4e6, z_end=500.0, sources=(;), kws...)
    grid = thermalization_grid(c, rec; nx)
    sol = solve_thermalization(grid; z_start, z_end, sources, kws...)
    sol.retcode == ReturnCode.Success ||
        error("exact thermalization solve failed: $(sol.retcode)")
    # Also hand back the evolved spectrum itself: μ and y are projections and
    # therefore convention-bound, but Δn(x) is the observable.
    merge(distortion_from_spectrum(grid, sol.u[end]), (x=grid.x, Δn=sol.u[end]))
end
