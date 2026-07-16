"""
CMB spectral distortions: μ, y, the temperature shift g, and the residual (r-type)
distortion, by the exact Green's-function / PCA method of Chluba & Jeong (2014).

The CMB is a blackbody to ~1e-4. The residual is not noise: it records how much
energy entered the photon bath and *when*, because whether an injection
thermalizes into a blackbody, a μ-distortion, or a y-distortion depends on the
redshift at which it happened. The physics is in `docs/src/spectral_distortions.md`.

# No approximation

Earlier versions of this file used a sharp μ/y visibility split. That is gone. The
branching between distortion types is now taken from the *exact* thermalization
Green's function G_th(x, z) that Chluba computed with CosmoTherm by solving the full
photon Boltzmann equation (Compton + double Compton + bremsstrahlung). The
redshift-dependent branching ratios J_T(z), J_y(z), J_μ(z) and the residual PCA
eigenvectors E_k(z) -- the projection of G_th onto the orthonormalised
(G, M, Y_SZ) spectral basis (Chluba & Jeong 2014, App. A) -- are shipped in
`data/FIRAS_branching_ratios.dat`, and the matching spectral shapes in
`data/FIRAS_distortions_shapes.dat`. These are cosmology-independent.

The heating history is computed here from first principles: Silk damping of
acoustic modes (reproducing CLASS's `noninjection_rate_acoustic_diss` integrand
exactly) plus adiabatic cooling of the baryons.

The distortion amplitudes are then

    X = ∫ dz  J_X(z)  (dQ/dz)/ργ,    X ∈ {g, y, μ, residual_k}

Validated against CLASS 3.3.4 (`sd_branching_approx = exact`) for the Planck 2018
cosmology: μ = 2.09e-8 vs 1.94e-8 (8%), y = 3.38e-9 vs 3.47e-9 (3%). The
temperature shift g is UV-sensitive (it draws on the thermalization-era heating
above z ≈ 1e6, where CLASS switches to a separate WKB treatment) and is reproduced
only approximately; it is a monopole shift, absorbed into the definition of T_CMB,
and not an observable.
"""

using QuadGK: quadgk
using Interpolations: interpolate, Gridded, Linear

# --- Data tables (cosmology-independent, lazily loaded) ---------------------

const _DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
const _BRANCHING = Ref{Any}(nothing)
const _SHAPES = Ref{Any}(nothing)
const _FIRAS = Ref{Any}(nothing)

"Read a CLASS distortions .dat file (skip comments and the leading `Nrows Ncols` line)."
function _read_class_dat(path)
    lines = filter(l -> !startswith(strip(l), "#") && !isempty(strip(l)), readlines(path))
    reduce(vcat, [permutedims(parse.(Float64, split(strip(l)))) for l in lines[2:end]])
end

"""
Branching ratios J_T(z), J_y(z), J_μ(z) and residual eigenvectors E_k(z), as
interpolants over log(1+z). Columns of the data file: z, J_T, J_y, J_μ, E_1..E_6.
"""
function _branching()
    if _BRANCHING[] === nothing
        M = _read_class_dat(joinpath(_DATA_DIR, "FIRAS_branching_ratios.dat"))
        lnz = log.(M[:, 1] .+ 1)
        itp(col) = interpolate((lnz,), M[:, col], Gridded(Linear()))
        nE = size(M, 2) - 4
        _BRANCHING[] = (
            z=M[:, 1], zmin=M[1, 1], zmax=M[end, 1],
            J_T=itp(2), J_y=itp(3), J_mu=itp(4),
            E=[itp(4 + k) for k in 1:nE], nE=nE,
        )
    end
    _BRANCHING[]
end

"""
Spectral shapes G_T(ν), Y_SZ(ν), M_μ(ν) and residual signals S_k(ν), as
interpolants over frequency in GHz. Units: 10⁻¹⁸ W/(m² Hz sr).
"""
function _shapes()
    if _SHAPES[] === nothing
        M = _read_class_dat(joinpath(_DATA_DIR, "FIRAS_distortions_shapes.dat"))
        ν = M[:, 1]
        itp(col) = interpolate((ν,), M[:, col], Gridded(Linear()))
        nS = size(M, 2) - 4
        _SHAPES[] = (
            ν=ν, νmin=ν[1], νmax=ν[end],
            G_T=itp(2), Y_SZ=itp(3), M_mu=itp(4),
            S=[itp(4 + k) for k in 1:nS], nS=nS,
        )
    end
    _SHAPES[]
end

"FIRAS monopole uncertainty σ(ν) in 10⁻²⁶ W/(m² Hz sr), from Fixsen et al. 1996."
function _firas()
    if _FIRAS[] === nothing
        M = _read_class_dat(joinpath(_DATA_DIR, "FIRAS_nu_delta_I.dat"))
        _FIRAS[] = (ν=M[:, 1], σ=M[:, 2])
    end
    _FIRAS[]
end

# --- Heating sources --------------------------------------------------------

"""
    silk_heating_rate(c, rec, z)

Fractional heating rate (dQ/dz)/ρ_γ from Silk damping of acoustic modes,
reproducing CLASS's `noninjection_rate_acoustic_diss` integrand exactly
(Chluba, Khatri & Sunyaev 2012):

    (dQ/dz)/ργ = ∫ dk  4 A²  k  Δ²_ℛ(k)  e^{-2k²/k_D²}  · 𝒟(z)

    A    = 1 / (1 + 4f_ν/15)                        neutrino anisotropic-stress factor
    Δ²_ℛ = A_s (k/k_pivot)^(n_s-1)                  dimensionless primordial spectrum
    𝒟(z) = 1/(H κ̇) · [16/15 + R²/(1+R)] / [6(1+R)]   how fast the damping scale sweeps

k_D grows as the diffusion length, so the heating at redshift z probes P(k) at
k ≈ k_D(z). In the μ era that is k ~ 10²–10⁴ /Mpc -- scales no other observable
reaches, which is why a μ measurement constrains inflation.
"""
function silk_heating_rate(c::Cosmology, rec::RecombinationSolution, z;
    kmin=0.12, kmax=1.0e6)
    a = scale_factor(z)
    R = R_baryon(c, a)
    Hz = H_Mpc(c, a)                       # proper H in 1/Mpc
    κ̇ = -τ_dot(rec, a)                     # a nₑ σ_T in 1/Mpc (CLASS index_th_dkappa)
    𝒟 = 1 / (Hz * κ̇) * (16 / 15 + R^2 / (1 + R)) / (6 * (1 + R))   # Mpc²
    kD = damping_scale(c, rec, a)

    ρν = _Ω_or_zero(c, MasslessNeutrinos) / a^4
    ργ = Ω_γ(c) / a^4
    f_ν = ρν / (ρν + ργ)
    A = 1 / (1 + 4 * f_ν / 15)

    # Integrate in u = k/k_D, not k. The envelope e^{-2u²} sits at u ~ 1 whatever
    # k_D is, and k_D sweeps five decades with redshift (0.3/Mpc near recombination,
    # 4000/Mpc in the μ era). A fixed k-grid that resolves the envelope in the μ
    # era misses it entirely near recombination -- the adaptive quadrature then
    # samples only k ≫ k_D, sees zero, and reports zero heating, silently killing
    # the whole y-era signal. In u the integrand is always well sampled. The upper
    # cap u = 10 is where e^{-2u²} ≈ 1e-87; nothing physical lives beyond it.
    integrand(u) = begin
        k = u * kD
        4 * A^2 * k * primordial_power(c, k) * exp(-2 * u^2) * 𝒟 * kD
    end
    lo = kmin / kD
    hi = min(kmax / kD, 10.0)
    hi <= lo && return 0.0
    quadgk(integrand, lo, hi; rtol=1e-8)[1]
end

"""
    damping_scale(c, rec, a)

Photon diffusion (Silk) wavenumber k_D, in CLASS's k_D = 2π/r_d normalisation
(envelope e^{-2k²/k_D²}):

    k_D⁻² = ∫₀^a da'/(a'³ H nₑ σ_T) · [16/15 + R²/(1+R)]/[6(1+R)]

Checked against CLASS's `r_d` (index_th_r_d): agrees to <0.5% across 10³ < z < 10⁶.
"""
function damping_scale(c::Cosmology, rec::RecombinationSolution, a; rtol=1e-6)
    integrand(a′) = begin
        z′ = redshift(a′)
        R = R_baryon(c, a′)
        ne = n_e(rec, z′) * Constants.Mpc_SI^3          # 1/Mpc³
        (16 / 15 + R^2 / (1 + R)) / (6 * (1 + R)) /
        (a′^3 * H_Mpc(c, a′) * ne * Constants.σ_T_Mpc2)
    end
    inv_kD2 = quadgk(integrand, 1e-8, a; rtol)[1]
    1 / sqrt(inv_kD2)
end

"""
    adiabatic_cooling_rate(c, rec, z)

Fractional heating (dQ/dz)/ργ from baryon cooling -- a *negative* source. The
baryons cool adiabatically as (1+z)² and, while Compton coupling holds, drain heat
from the photons: the one guaranteed source of the opposite sign to Silk damping.

    (dQ/dz)/ργ = -(3/2) N_tot k_B T_γ / ργ · 1/(1+z),   N_tot = nₑ + n_H + n_He
"""
function adiabatic_cooling_rate(c::Cosmology, rec::RecombinationSolution, z)
    T_γ = c.Tcmb * (1 + z)
    ργ_SI = Constants.a_rad_SI * T_γ^4
    nH = n_H_of_z(c, z)
    N_tot = nH * (1 + f_He(c) + x_e(rec, z))
    -(3 / 2) * N_tot * Constants.k_B_SI * T_γ / ργ_SI / (1 + z)
end

"""
    heating_rate(c, rec, z; silk = true, cooling = true, nuclear = true)

Total fractional heating rate (dQ/dz)/ρ_γ, summing the guaranteed sources.
"""
function heating_rate(c::Cosmology, rec::RecombinationSolution, z;
    silk=true, cooling=true, nuclear=true, nuclear_history=nothing)
    q = 0.0
    silk && (q += silk_heating_rate(c, rec, z))
    cooling && (q += adiabatic_cooling_rate(c, rec, z))
    nuclear && (q += nuclear_heating_rate(c, rec, z; history=nuclear_history))
    q
end

# --- Distortion amplitudes --------------------------------------------------

"""
    distortion_amplitudes(c, rec; silk = true, cooling = true, nuclear = true)

Compute all spectral-distortion amplitudes at once by convolving the heating
history with the exact branching ratios:

    X = ∫ dz  J_X(z)  (dQ/dz)/ργ

Returns a NamedTuple `(g, y, μ, residual)` where `residual` is the vector of PCA
residual amplitudes μ_k. The integral runs over the redshift range of the
branching-ratio table (z ≈ 10³–5×10⁶), which brackets the μ and y eras.
"""
function distortion_amplitudes(c::Cosmology, rec::RecombinationSolution;
    silk=true, cooling=true, nuclear=true, rtol=1e-5)
    br = _branching()
    lo, hi = log(br.zmin + 1), log(br.zmax + 1)
    nh = nuclear ? trace_nuclear_history(c, rec) : nothing

    # ∫ J_X(z) rate(z) dz, integrated in ln(1+z) with the (1+z) Jacobian.
    conv(J) = quadgk(lo, hi; rtol) do lz
        z = exp(lz) - 1
        J(lz) * heating_rate(c, rec, z; silk, cooling, nuclear,
            nuclear_history=nh) * (z + 1)
    end[1]

    (
        g=conv(br.J_T),
        y=conv(br.J_y),
        μ=conv(br.J_mu),
        residual=[conv(Ek) for Ek in br.E],
    )
end

"μ chemical-potential distortion. Planck-cosmology ΛCDM ≈ 2e-8; FIRAS |μ| < 9e-5."
μ_distortion(c::Cosmology, rec::RecombinationSolution; kws...) =
    distortion_amplitudes(c, rec; kws...).μ

"y (Compton) distortion from the pre-recombination sources. FIRAS |y| < 1.5e-5."
y_distortion(c::Cosmology, rec::RecombinationSolution; kws...) =
    distortion_amplitudes(c, rec; kws...).y

"Temperature-shift (g-type) amplitude. A monopole shift, absorbed into T_CMB."
g_distortion(c::Cosmology, rec::RecombinationSolution; kws...) =
    distortion_amplitudes(c, rec; kws...).g

# --- Distortion spectrum ----------------------------------------------------

"""
    distortion_spectrum(c, rec, ν_GHz; kws...)

The spectral distortion ΔI(ν) in 10⁻²⁶ W/(m² Hz sr), summing every component:

    ΔI(ν) = g·G_T(ν) + y·Y_SZ(ν) + μ·M_μ(ν) + Σ_k μ_k S_k(ν)

The residual sum Σ μ_k S_k is the r-type distortion -- the part of the signal that
is neither a pure μ, y, nor temperature shift, and which distinguishes energy-
release *histories* rather than just totals. `ν_GHz` may be a scalar or array.
"""
function distortion_spectrum(c::Cosmology, rec::RecombinationSolution, ν_GHz; kws...)
    amp = distortion_amplitudes(c, rec; kws...)
    sh = _shapes()
    # Shapes are in 1e-18 W/m²/Hz/sr; convert to the 1e-26 units of ΔI (×1e8).
    _one(ν) = 1e8 * (amp.g * sh.G_T(ν) + amp.y * sh.Y_SZ(ν) + amp.μ * sh.M_mu(ν) +
                     sum(amp.residual[k] * sh.S[k](ν) for k in 1:sh.nS))
    ν_GHz isa Number ? _one(ν_GHz) : _one.(ν_GHz)
end

# --- Observational constraints ----------------------------------------------

"""
    firas_chi2(c, rec; kws...)

χ² of the predicted distortion against the COBE/FIRAS null result (Fixsen et al.
1996): Σ_i [ΔI(ν_i)/σ_i]². A ΛCDM distortion is ~4 orders of magnitude below the
FIRAS noise, so this is ≪ 1 -- FIRAS does not constrain standard cosmology, which
is exactly why a next-generation spectrometer (PIXIE) is the goal.
"""
function firas_chi2(c::Cosmology, rec::RecombinationSolution; kws...)
    f = _firas()
    ΔI = distortion_spectrum(c, rec, f.ν; kws...)
    sum((ΔI ./ f.σ) .^ 2)
end

# Observational limits, for convenience.
const FIRAS_μ_limit = 9e-5      # 95% CL, Fixsen et al. 1996
const FIRAS_y_limit = 1.5e-5    # 95% CL
const PIXIE_μ_target = 1e-8     # forecast sensitivity
