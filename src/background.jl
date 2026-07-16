"""
Background (homogeneous, isotropic) cosmology.

A `Cosmology` is a bag of [`Species`](@ref) plus the handful of numbers that
are not species-specific (h, T_cmb, helium fraction, primordial spectrum). The
expansion rate is a sum over the bag, and every distance and time measure is
built from that one function.
"""

using QuadGK: quadgk
using OrdinaryDiffEq: ODEProblem, solve, Tsit5, Vern9
using Interpolations: linear_interpolation, Line
using Roots: find_zero, Bisection

"""
    Cosmology

Container for a background model. Build one with [`cosmology`](@ref) rather
than calling this constructor directly.

Fields beyond the species list:
- `h`      : H0 = 100h km/s/Mpc
- `Tcmb`   : CMB temperature today, K
- `Yp`     : primordial helium mass fraction (needed by recombination)
- `A_s`,`n_s`,`k_pivot` : primordial curvature power spectrum,
  P_R(k) = A_s (k/k_pivot)^(n_s - 1)
"""
struct Cosmology{T<:Real,B,P,S<:Tuple}
    h::T
    Tcmb::T
    Yp::T
    A_s::T
    n_s::T
    k_pivot::T
    bbn_result::B
    primordial::P
    species::S
end

# Backward-compatible positional constructors used by older examples and user
# code. An explicitly supplied Yp does not determine the trace nuclide
# abundances, so the absence of a BBN solution is `nothing`; `primordial` is the
# optional override of the power-law primordial spectrum (see `primordial_power`)
# and defaults to `nothing` (use A_s/n_s/k_pivot).
Cosmology(h, Tcmb, Yp, A_s, n_s, k_pivot, species::Tuple) =
    Cosmology(h, Tcmb, Yp, A_s, n_s, k_pivot, nothing, nothing, species)
Cosmology(h, Tcmb, Yp, A_s, n_s, k_pivot, bbn_result, species::Tuple) =
    Cosmology(h, Tcmb, Yp, A_s, n_s, k_pivot, bbn_result, nothing, species)

# --- Species lookup ---------------------------------------------------------

"""
    get_species(c, Type)

Return the first species of the given type, or `nothing`. Used by the
thermodynamics and perturbation layers, which need e.g. the baryons
specifically rather than just the total density.
"""
function get_species(c::Cosmology, ::Type{S}) where {S<:Species}
    for s in c.species
        s isa S && return s
    end
    return nothing
end

"All species of a given type (there can be several massive neutrinos)."
get_all_species(c::Cosmology, ::Type{S}) where {S<:Species} =
    filter(s -> s isa S, collect(c.species))

_Ω_or_zero(c, S) = (s = get_species(c, S); s === nothing ? 0.0 : s.Ω)

Ω_b(c::Cosmology) = _Ω_or_zero(c, Baryons)
Ω_c(c::Cosmology) = _Ω_or_zero(c, ColdDarkMatter)
Ω_γ(c::Cosmology) = _Ω_or_zero(c, Photons)
Ω_k(c::Cosmology) = _Ω_or_zero(c, Curvature)
Ω_m(c::Cosmology) = Ω_b(c) + Ω_c(c)

"Dark-energy density today, whatever the flavour."
Ω_de(c::Cosmology) = sum(Ω0(s) for s in c.species if s isa AbstractDarkEnergy; init=0.0)

"Total neutrino density today, massless plus massive."
function Ω_ν(c::Cosmology)
    tot = 0.0
    for s in c.species
        (s isa MasslessNeutrinos || s isa MassiveNeutrinos) && (tot += Ω0(s))
    end
    tot
end

"Radiation density today (photons + anything still relativistic)."
Ω_r(c::Cosmology) = Ω_γ(c) + _Ω_or_zero(c, MasslessNeutrinos)

# --- Expansion rate ---------------------------------------------------------

"""
    E(c, a)

Dimensionless expansion rate H(a)/H0, from summing every species:

    E(a)² = Σ_s ρ_s(a)/ρ_c0
"""
function E(c::Cosmology, a)
    tot = zero(a)
    for s in c.species
        tot += ρ_over_ρc0(s, a)
    end
    sqrt(tot)
end

"H(a)/H0 as a function of redshift instead of scale factor."
E_z(c::Cosmology, z) = E(c, scale_factor(z))

"Hubble rate at redshift z, in km/s/Mpc."
H(c::Cosmology, z) = 100 * c.h * E_z(c, z)

"Hubble rate in 1/Mpc (the internal unit), as a function of scale factor."
H_Mpc(c::Cosmology, a) = Constants.H0_in_invMpc(c.h) * E(c, a)

"Conformal Hubble rate ℋ = a H, in 1/Mpc."
ℋ(c::Cosmology, a) = a * H_Mpc(c, a)

scale_factor(z) = 1 / (1 + z)
redshift(a) = 1 / a - 1

"Hubble distance c/H0 in Mpc."
hubble_distance(c::Cosmology) = Constants.hubble_distance_Mpc(c.h)

"Hubble time 1/H0 in Gyr."
hubble_time(c::Cosmology) = Constants.hubble_time_Gyr(c.h)

# --- Distances --------------------------------------------------------------
#
# All distance integrals are done in `a` rather than `z`. The integrand
#   da / (a² E(a))
# is smooth and finite as a → 0 (a²E → √Ω_r there), whereas the z-integrand
# blows up at high z. This matters: the comoving distance to last scattering
# is one of the numbers the CMB layer leans on hardest.

_comoving_integrand(c, a) = 1 / (a^2 * E(c, a))

"""
    comoving_distance(c, z₁, z₂ = nothing)

Line-of-sight comoving distance in Mpc. With one redshift, the distance from
us to `z₁`; with two, the distance between them.
"""
function comoving_distance(c::Cosmology, z::Real; rtol=1e-10)
    hubble_distance(c) * quadgk(a -> _comoving_integrand(c, a), scale_factor(z), 1.0; rtol)[1]
end

function comoving_distance(c::Cosmology, z₁::Real, z₂::Real; rtol=1e-10)
    hubble_distance(c) *
    quadgk(a -> _comoving_integrand(c, a), scale_factor(z₂), scale_factor(z₁); rtol)[1]
end

"""
    transverse_comoving_distance(c, z)

Comoving distance transverse to the line of sight, D_M. This is where spatial
curvature finally shows up: for Ω_k ≠ 0 a comoving separation subtends a
different angle than it would in flat space, so

    D_M = D_H/√Ω_k  · sinh(√Ω_k · D_C/D_H)     (open,   Ω_k > 0)
    D_M = D_C                                   (flat,   Ω_k = 0)
    D_M = D_H/√|Ω_k| · sin (√|Ω_k| · D_C/D_H)   (closed, Ω_k < 0)

Angular-diameter and luminosity distances are built on D_M, not D_C, which is
the piece the old code was missing.
"""
function transverse_comoving_distance(c::Cosmology, z::Real; kws...)
    Ωk = Ω_k(c)
    D_C = comoving_distance(c, z; kws...)
    D_H = hubble_distance(c)
    if Ωk > 0
        sq = sqrt(Ωk)
        return D_H / sq * sinh(sq * D_C / D_H)
    elseif Ωk < 0
        sq = sqrt(-Ωk)
        return D_H / sq * sin(sq * D_C / D_H)
    else
        return D_C
    end
end

"""
    angular_diameter_distance(c, z)

D_A = D_M/(1+z), in Mpc. Physical size = D_A × angle subtended.
"""
angular_diameter_distance(c::Cosmology, z::Real; kws...) =
    transverse_comoving_distance(c, z; kws...) / (1 + z)

"""
    angular_diameter_distance(c, z₁, z₂)

Angular-diameter distance between two redshifts, the quantity gravitational
lensing needs. Note this is *not* D_A(z₂) - D_A(z₁); the curvature-correct
form (Hogg 1999 eq. 19) is used.
"""
function angular_diameter_distance(c::Cosmology, z₁::Real, z₂::Real; kws...)
    Ωk = Ω_k(c)
    D_H = hubble_distance(c)
    D_M1 = transverse_comoving_distance(c, z₁; kws...)
    D_M2 = transverse_comoving_distance(c, z₂; kws...)
    if Ωk ≈ 0
        return (D_M2 - D_M1) / (1 + z₂)
    end
    # Curved case: propagate through the hyperbolic/trigonometric addition law.
    term1 = D_M2 * sqrt(1 + Ωk * D_M1^2 / D_H^2)
    term2 = D_M1 * sqrt(1 + Ωk * D_M2^2 / D_H^2)
    return (term1 - term2) / (1 + z₂)
end

"""
    luminosity_distance(c, z)

D_L = (1+z) D_M, in Mpc. Flux = L / (4π D_L²).
"""
luminosity_distance(c::Cosmology, z::Real; kws...) =
    (1 + z) * transverse_comoving_distance(c, z; kws...)

"""
    distance_modulus(c, z)

μ = 5 log₁₀(D_L / 10 pc), the quantity supernova cosmology is done in.
"""
distance_modulus(c::Cosmology, z::Real; kws...) =
    5 * log10(luminosity_distance(c, z; kws...) * 1e5)  # D_L in Mpc -> units of 10pc

"""
    comoving_volume(c, z)

All-sky comoving volume out to `z`, in Mpc³, with the curvature correction
(Hogg 1999 eq. 29).
"""
function comoving_volume(c::Cosmology, z::Real; kws...)
    Ωk = Ω_k(c)
    D_H = hubble_distance(c)
    D_M = transverse_comoving_distance(c, z; kws...)
    if Ωk ≈ 0
        return 4π / 3 * D_M^3
    end
    x = sqrt(abs(Ωk)) * D_M / D_H
    pref = 4π * D_H^3 / (2 * Ωk)
    if Ωk > 0
        return pref * (D_M / D_H * sqrt(1 + Ωk * D_M^2 / D_H^2) - asinh(x) / sqrt(Ωk))
    else
        return pref * (D_M / D_H * sqrt(1 + Ωk * D_M^2 / D_H^2) - asin(x) / sqrt(-Ωk))
    end
end

# --- Times ------------------------------------------------------------------

"""
    age(c, z; unit = :Gyr)

Age of the universe at redshift `z`: t = ∫₀^a da/(a H(a)).
"""
function age(c::Cosmology, z::Real=0.0; rtol=1e-10)
    a_end = scale_factor(z)
    integral = quadgk(a -> 1 / (a * E(c, a)), 0.0, a_end; rtol)[1]
    hubble_time(c) * integral   # Gyr
end

"""
    lookback_time(c, z)

Time elapsed between redshift `z` and today, in Gyr.
"""
function lookback_time(c::Cosmology, z::Real; rtol=1e-10)
    integral = quadgk(a -> 1 / (a * E(c, a)), scale_factor(z), 1.0; rtol)[1]
    hubble_time(c) * integral
end

"""
    conformal_time(c, a)

Conformal time η = ∫₀^a da'/(a'² H(a')), in Mpc. This is the time variable the
Boltzmann hierarchy is integrated in, and it doubles as the comoving particle
horizon.
"""
function conformal_time(c::Cosmology, a::Real; rtol=1e-10)
    hubble_distance(c) * quadgk(x -> _comoving_integrand(c, x), 0.0, a; rtol)[1]
end

"Comoving distance to the horizon today, η₀, in Mpc."
conformal_time_today(c::Cosmology; kws...) = conformal_time(c, 1.0; kws...)

# --- Characteristic epochs --------------------------------------------------

"""
    ρ_matter(c, a), ρ_radiation(c, a)

Non-relativistic and relativistic energy densities, in units of ρ_c0.

Massive neutrinos are counted in `ρ_radiation` at their *actual* density ρ_ν(a),
not as a⁻³ matter. This is not a fudge: at equality (z ≈ 3400) the neutrino
temperature is ~0.6 eV, so a 0.06 eV neutrino is still ultra-relativistic and
genuinely behaves as radiation. Their mass only matters much later, and by then
radiation is irrelevant. Sorting species into "matter" and "radiation" by *type*
instead of by what they are actually doing at scale factor `a` is what makes the
textbook Ω_m/Ω_r shortcut wrong.
"""
ρ_matter(c::Cosmology, a) =
    sum(ρ_over_ρc0(s, a) for s in c.species if s isa Baryons || s isa ColdDarkMatter; init=0.0)

ρ_radiation(c::Cosmology, a) =
    sum(ρ_over_ρc0(s, a) for s in c.species
        if s isa Photons || s isa MasslessNeutrinos || s isa MassiveNeutrinos; init=0.0)

"""
    z_equality(c)

Redshift of matter-radiation equality, from root-solving ρ_m(a) = ρ_r(a).
"""
function z_equality(c::Cosmology)
    f(loga) = (a = exp(loga); log(ρ_matter(c, a) / ρ_radiation(c, a)))
    loga_eq = find_zero(f, (log(1e-10), log(1.0)), Bisection())
    redshift(exp(loga_eq))
end

"""
    scale_factor_of_time(c, t_Gyr)

Invert t(a) to get a(t). Replaces the old ODE-based `scale_fact`, which
integrated da/dt from a = 1e-8 and inherited whatever error that seeding
introduced; root-solving the age integral is both faster and exact.
"""
function scale_factor_of_time(c::Cosmology, t_Gyr::Real)
    t_Gyr <= 0 && return 0.0
    f(loga) = age(c, redshift(exp(loga))) - t_Gyr
    # Bracket generously: a ∈ [1e-12, 1e3] covers radiation era to far future.
    loga = find_zero(f, (log(1e-12), log(1e3)), Bisection())
    exp(loga)
end

# --- Constructor ------------------------------------------------------------

# Solve the quintessence background: the Klein–Gordon equation for a canonical
# field φ (reduced Planck units) in potential A·Ṽ(φ) (Ṽ in units of ρ_c0),
# with the amplitude A shot by bisection so the field's density today closes
# the budget at Ω_target. In e-folds N with π ≡ φ̇/H0 the system is explicit —
#   dφ/dN = π/E,   dπ/dN = −3π − 3AṼ'(φ)/E,   E² = ρ_static(a) + π²/6 + AṼ(φ)
# — so no self-consistency iteration is needed. IC: frozen field (π = 0, the
# early-time attractor under Hubble friction).
function _solve_scf_history(V, dV, φ0, ρ_static, Ω_target; x_ini=log(1e-9), npts=2000)
    ρscf_today(lnA) = begin
        A = exp(lnA)
        function rhs!(du, u, p, x)
            φ, π = u
            E2 = ρ_static(exp(x)) + π^2 / 6 + A * V(φ)
            E = sqrt(max(E2, 1e-300))
            du[1] = π / E
            du[2] = -3π - 3A * dV(φ) / E
        end
        sol = solve(ODEProblem(rhs!, [float(φ0), 0.0], (x_ini, 0.0)),
            Tsit5(); reltol=1e-10, abstol=1e-12)
        φ1, π1 = sol(0.0)
        (π1^2 / 6 + A * V(φ1), sol, A)
    end
    # bracket the amplitude: frozen field ⇒ ρ_today ≈ A·V(φ0), so seed there
    lnA = log(Ω_target / max(V(φ0), 1e-30))
    lo, hi = lnA - 15, lnA + 15
    (ρlo, _, _) = ρscf_today(lo)
    (ρhi, _, _) = ρscf_today(hi)
    (ρlo < Ω_target < ρhi) || error(
        "quintessence shooting failed to bracket Ω_DE = $(round(Ω_target,digits=4)): " *
        "ρ(A_lo)=$(round(ρlo,sigdigits=3)), ρ(A_hi)=$(round(ρhi,sigdigits=3)). " *
        "Check V_scf/φ0_scf (V in units of ρ_c0, φ in reduced Planck masses).")
    local sol, A
    for _ in 1:60
        mid = 0.5 * (lo + hi)
        ρm, sol, A = ρscf_today(mid)
        (abs(ρm / Ω_target - 1) < 1e-10) && break
        ρm < Ω_target ? (lo = mid) : (hi = mid)
    end
    xs = collect(range(x_ini, 0.0; length=npts))
    lnρ = similar(xs)
    wv = similar(xs)
    for (i, x) in enumerate(xs)
        φ, π = sol(x)
        kin = π^2 / 6
        pot = A * V(φ)
        lnρ[i] = log(kin + pot)
        wv[i] = (kin - pot) / (kin + pot)
    end
    itp(y) = linear_interpolation(xs, y; extrapolation_bc=Line())
    QuintessenceDE(exp(lnρ[end]), itp(lnρ), itp(wv))
end

# Solve the decaying-DM background: coupled ln-densities of the decaying matter
# and its dark-radiation product, self-consistent in H (which contains them) and
# in the dark-energy closure (fixed point — converges in a few passes because DE
# barely feeds back on the early expansion). Γ in 1/Mpc; Ω_ini is the density
# parameter the dcdm *would* have today with Γ = 0.
function _solve_dcdm_history(Γ, Ω_ini, ρ_static, de_shape, h; x_ini=log(1e-9), npass=4)
    H0 = Constants.H0_in_invMpc(h)
    Ω_static0 = ρ_static(1.0)
    Ω_DE = 1.0 - Ω_static0 - Ω_ini            # first guess: as if no decay
    local sol, Ωd0, Ωr0
    for _ in 1:npass
        H(x, ρd, ρr) = H0 * sqrt(ρ_static(exp(x)) + ρd + ρr + Ω_DE * de_shape(exp(x)))
        function rhs!(du, u, p, x)
            ρd, ρr = exp(u[1]), exp(u[2])
            g = Γ / H(x, ρd, ρr)               # aΓ/ℋ = Γ/H
            du[1] = -3 - g
            du[2] = -4 + g * ρd / ρr
        end
        a_i = exp(x_ini)
        ρd_i = Ω_ini / a_i^3
        # early attractor for the product: dρr/dx + 4ρr = gρd with ρd ∝ a⁻³ in RD
        # (g ∝ a²) gives ρr = gρd/3
        g_i = Γ / (H0 * sqrt(ρ_static(a_i) + ρd_i))
        ρr_i = max(g_i * ρd_i / 3, 1e-300)
        sol = solve(ODEProblem(rhs!, [log(ρd_i), log(ρr_i)], (x_ini, 0.0)),
            Tsit5(); reltol=1e-10, abstol=1e-12)
        Ωd0 = exp(sol(0.0)[1])
        Ωr0 = exp(sol(0.0)[2])
        Ω_DE = 1.0 - Ω_static0 - Ωd0 - Ωr0
    end
    s = sol
    xi = x_ini
    lnρd = x -> s(clamp(x, xi, 0.0))[1]
    lnρr = x -> s(clamp(x, xi, 0.0))[2]
    (; lnρd, lnρr, Ω_dcdm0=Ωd0, Ω_dr0=Ωr0)
end

"""
    cosmology(; h, Ω_b, Ω_c, Ω_k, Tcmb, N_eff, m_ν, w0, wa, ...)

Build a `Cosmology`. Defaults are Planck 2018 TT,TE,EE+lowE+lensing+BAO.

Dark energy closes the budget: Ω_de = 1 - Ω_k - Σ(everything else), evaluated
*today* so that massive neutrinos (which are matter now but were radiation
once) are counted correctly.

Neutrinos: `m_ν` is a vector of masses in eV. Species with mass are given the
full phase-space treatment; the remainder of `N_eff` is carried by massless
species. The default `m_ν = [0.06, 0.0, 0.0]` is the minimal normal hierarchy.
Pass `m_ν = Float64[]` for the (unphysical, but common in older codes) case of
three massless neutrinos.

Dark energy: `w0`/`wa` select CPL. For anything more exotic, pass `w_func`,
an arbitrary function of `a`.

# Examples
```julia
c = cosmology()                              # Planck 2018 ΛCDM
c = cosmology(Ω_k = 0.01)                    # open
c = cosmology(w0 = -0.9, wa = -0.3)          # CPL dark energy
c = cosmology(m_ν = [0.1, 0.1, 0.1])         # degenerate 0.3 eV neutrinos
c = cosmology(w_func = a -> -1 + 0.1*(1-a)^2)  # arbitrary w(a)
```
"""
function cosmology(;
    h=0.6766,
    Ω_b=0.04897,
    Ω_c=0.2607,
    Ω_k=0.0,
    Tcmb=2.7255,
    N_eff=3.044,
    m_ν=[0.06, 0.0, 0.0],
    ξ_ν=nothing,
    m_wdm=nothing,
    Ω_wdm=0.0,
    Γ_dcdm=0.0,
    Ω_dcdm_ini=0.0,
    V_scf=nothing,
    dV_scf=nothing,
    φ0_scf=1.0,
    Yp=:bbn,
    w0=-1.0,
    wa=0.0,
    w_func=nothing,
    A_s=2.105e-9,
    n_s=0.9665,
    k_pivot=0.05,
    primordial=nothing,
)
    # Yp = :bbn solves the nuclear network for the helium fraction instead of taking
    # it as an input. It is not really a free parameter -- given ω_b and N_eff, Big
    # Bang nucleosynthesis fixes it -- and treating it as one lets a cosmology quietly
    # carry a helium abundance inconsistent with its own baryon density.
    # N_eff = :derived computes it rather than quoting 3.044 -- see `neff()`.
    N_eff === :derived && (N_eff = N_eff_sm())

    Ωγ = Constants.Ω_γ(Tcmb, h)
    photons = Photons(Ωγ)

    # Split N_eff between massive and massless species. Each of the N_ν slots
    # carries N_eff/N_ν worth of density; slots with m > 0 get the phase-space
    # treatment, the rest stay massless.
    masses = collect(float.(m_ν))
    N_slots = max(length(masses), 1)
    share = N_eff / 3.0   # per-slot share, using the standard 3-neutrino norm

    # Optional per-slot degeneracy parameter ξ = μ/T (lepton asymmetry). A scalar
    # applies to every slot; a vector is taken slot by slot; nothing means ξ = 0.
    ξs = ξ_ν === nothing ? zeros(length(masses)) :
         ξ_ν isa Number ? fill(float(ξ_ν), length(masses)) : collect(float.(ξ_ν))

    ν_species = Species[]
    N_massless_eff = 0.0
    for (i, m) in enumerate(masses)
        ξ = i <= length(ξs) ? ξs[i] : 0.0
        if m > 0
            push!(ν_species, massive_neutrino(m, Tcmb, h; N_eff_share=share, ξ))
        else
            # a degenerate massless slot carries extra radiation density
            N_massless_eff += share * degenerate_Neff_factor(ξ)
        end
    end
    if isempty(masses)
        N_massless_eff = N_eff
    end
    if N_massless_eff > 0
        push!(ν_species, massless_neutrinos(N_massless_eff, Ωγ))
    end

    baryons = Baryons(Ω_b)
    cdm = ColdDarkMatter(Ω_c)

    # Optional thermal warm-dark-matter relic. It rides the massive-neutrino
    # phase-space machinery (full momentum-resolved hierarchy), so it lands in
    # the ν_species list; its Ω today is Ω_wdm by construction.
    if Ω_wdm > 0
        m_wdm === nothing && throw(ArgumentError(
            "Ω_wdm > 0 needs m_wdm (the relic mass in keV)"))
        push!(ν_species, warm_dark_matter(m_wdm, Ω_wdm, Tcmb, h))
    end

    # Close the budget with dark energy, using densities evaluated today.
    Ω_other = Ωγ + Ω_b + Ω_c + Ω_k + sum(Ω0(s) for s in ν_species; init=0.0)

    # Optional decaying dark matter (dcdm → dark radiation). Its ρ histories are
    # solved here (no closed form exists) and its realised today-densities join
    # the budget before dark energy closes it.
    dcdm_hist = nothing
    Γmpc = 0.0
    if Ω_dcdm_ini > 0
        Γ_dcdm > 0 || throw(ArgumentError(
            "Ω_dcdm_ini > 0 needs Γ_dcdm > 0 (km/s/Mpc); for stable DM use Ω_c"))
        Γmpc = Γ_dcdm * Constants.km_SI / Constants.c_SI      # km/s/Mpc → 1/Mpc
        static = (Photons(Ωγ), Baryons(Ω_b), ColdDarkMatter(Ω_c), ν_species...,
            (Ω_k != 0 ? (Curvature(Ω_k),) : ())...)
        ρ_static = a -> sum(ρ_over_ρc0(s, a) for s in static)
        de_shape = if w_func !== nothing
            a -> ρ_over_ρc0(GeneralDarkEnergy(1.0, w_func), a)
        elseif w0 == -1.0 && wa == 0.0
            a -> 1.0
        else
            a -> ρ_over_ρc0(W0WaDarkEnergy(1.0, w0, wa), a)
        end
        dcdm_hist = _solve_dcdm_history(Γmpc, float(Ω_dcdm_ini), ρ_static, de_shape, h)
        Ω_other += dcdm_hist.Ω_dcdm0 + dcdm_hist.Ω_dr0
    end

    Ω_DE = 1.0 - Ω_other

    de = if V_scf !== nothing
        # quintessence replaces the closure species: solve the Klein–Gordon
        # background with the amplitude shot so ρ_scf(a=1) = Ω_DE exactly.
        dV_scf === nothing && throw(ArgumentError(
            "V_scf needs dV_scf (the potential derivative dṼ/dφ)"))
        (w0 == -1.0 && wa == 0.0 && w_func === nothing) || throw(ArgumentError(
            "V_scf is exclusive with w0/wa/w_func — the field determines w(a)"))
        dcdm_hist === nothing || throw(ArgumentError(
            "V_scf together with dcdm is not supported yet"))
        static = (Photons(Ωγ), Baryons(Ω_b), ColdDarkMatter(Ω_c), ν_species...,
            (Ω_k != 0 ? (Curvature(Ω_k),) : ())...)
        _solve_scf_history(V_scf, dV_scf, φ0_scf,
            a -> sum(ρ_over_ρc0(s, a) for s in static), Ω_DE)
    elseif w_func !== nothing
        GeneralDarkEnergy(Ω_DE, w_func)
    elseif w0 == -1.0 && wa == 0.0
        CosmologicalConstant(Ω_DE)
    else
        W0WaDarkEnergy(Ω_DE, w0, wa)
    end

    all_species = Species[photons, baryons, cdm]
    append!(all_species, ν_species)
    if dcdm_hist !== nothing
        push!(all_species, DecayingCDM(Γmpc, dcdm_hist.lnρd))
        push!(all_species, DecayRadiation(dcdm_hist.lnρr))
    end
    Ω_k != 0 && push!(all_species, Curvature(Ω_k))
    push!(all_species, de)

    bbn_result = nothing
    if Yp === :bbn
        if w_func === nothing && V_scf === nothing
            # Cacheable exact composition: matter, curvature and the complete
            # CPL density history all enter the nuclear Friedmann rate.
            # dcdm counts as matter during BBN: Γt(BBN) ≪ 1 for any Γ expressed
            # in km/s/Mpc (1 km/s/Mpc ≈ 1e-19 s⁻¹, BBN at t ~ 10²-10³ s).
            bbn_result = _cached_standard_bbn(Ω_b * h^2,
                (Ω_b + Ω_c + Ω_dcdm_ini) * h^2,
                N_eff, 878.4; ω_k=Ω_k * h^2, ω_de=Ω_DE * h^2, w0, wa)
        else
            # An arbitrary w(a) cannot be keyed by a finite parameter tuple.
            # Build the composition first and let bbn(c) evaluate its exact ρ(a).
            hv = promote(float(h), float(Tcmb), 0.245, float(A_s),
                float(n_s), float(k_pivot))
            seed = Cosmology(hv..., nothing, Tuple(all_species))
            bbn_result = bbn(seed)
        end
        Yp = bbn_result.Y_p
    end

    Cosmology(
        promote(float(h), float(Tcmb), float(Yp), float(A_s), float(n_s), float(k_pivot))...,
        bbn_result,
        primordial,
        Tuple(all_species),
    )
end

function Base.show(io::IO, c::Cosmology)
    println(io, "Cosmology(h = $(c.h), T_cmb = $(c.Tcmb) K)")
    println(io, "  Ω_b  = $(round(Ω_b(c),  sigdigits=5))")
    println(io, "  Ω_c  = $(round(Ω_c(c),  sigdigits=5))")
    println(io, "  Ω_γ  = $(round(Ω_γ(c),  sigdigits=5))")
    println(io, "  Ω_ν  = $(round(Ω_ν(c),  sigdigits=5))")
    println(io, "  Ω_k  = $(round(Ω_k(c),  sigdigits=5))")
    println(io, "  Ω_de = $(round(Ω_de(c), sigdigits=5))")
    print(io, "  species: ", join(string.(nameof.(typeof.(c.species))), ", "))
end
