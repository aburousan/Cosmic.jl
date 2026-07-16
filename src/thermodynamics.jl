"""
Thermal history: recombination, the optical depth, and the visibility function.

The layer computes the free-electron fraction x_e(z) at the level of RECFAST 1.5
(Seager, Sasselov & Scott 1999/2000; Wong, Moss & Scott 2008): an effective
three-level atom for hydrogen and for neutral helium, a separately-evolved
matter temperature, and a set of correction terms calibrated against the full
multi-level-atom codes.

Above z ≈ 5000 the reaction rates beat the expansion by many orders of
magnitude, Saha equilibrium is exact, and we use it. Below that the rates lose
and the ODEs have to be integrated.

Three things a textbook Peebles treatment leaves out, each worth percent-level
accuracy on x_e:

  * **The Peebles C factor.** A recombination straight to the ground state emits
    a photon that immediately re-ionizes a neighbour, so it achieves nothing.
    Net recombination proceeds only through the 2s->1s two-photon decay and
    through Lyman-α photons redshifting out of the line.
  * **Helium**, roughly 24% of the baryon mass. He II -> He I near z ≈ 2000 and
    He III -> He II near z ≈ 6000 both leave their mark on x_e.
  * **The matter temperature**, which tracks the radiation only while Compton
    scattering is efficient, then decouples and falls as (1+z)² rather than
    (1+z). The recombination coefficients depend on it, so this feeds back.

On top of those, the corrections that separate RECFAST 1.5 from the 1999
original -- see [`recombination`](@ref) for what each one does and why: a
retuned hydrogen fudge factor and a double-Gaussian correction to the Lyman-α
escape rate, both matched to HyRec (Ali-Haïmoud & Hirata 2010) and CosmoRec;
and for helium, a Sobolev escape probability, H I continuum opacity (Switzer &
Hirata 2008), and the 2³P₁ -> 1¹S₀ intercombination channel.

Not included: a genuine multi-level atom. HyRec and CosmoRec solve ~100 shells
with resolved angular-momentum substates, and the correction terms above are
*fitted to them*, so what remains is ~0.1% on x_e. If you need that last 0.1%,
use HyRec.
"""

using OrdinaryDiffEq: ODEProblem, solve, Rodas5P
using Interpolations: linear_interpolation, Line
using Roots: find_zero, Bisection
using QuadGK: quadgk
using SpecialFunctions: erfc

const _C = Constants

# --- Rate coefficients ------------------------------------------------------

"""
    α_H(T)

Case-B recombination coefficient for hydrogen, m³/s. Pequignot, Petitjean &
Boisson (1991) fit, as used by RECFAST.

"Case B" means recombinations to the ground state are excluded: those produce
an ionizing photon that is immediately reabsorbed in the optically thick
medium, so they achieve nothing.
"""
function α_H(T, fudge=1.14)
    t = T / 1e4
    fudge * 1e-19 * (4.309 * t^(-0.6166)) / (1 + 0.6703 * t^0.5300)
end

"""
    α_He(T)

Case-B recombination coefficient for neutral helium, m³/s. Verner & Ferland
(1996) fit, as used by RECFAST.
"""
function α_He(T)
    a = 10.0^(-16.744)
    b = 0.711
    T0 = 10.0^0.477121     # 3 K
    T1 = 10.0^5.114
    sq0 = sqrt(T / T0)
    sq1 = sqrt(T / T1)
    a / (sq0 * (1 + sq0)^(1 - b) * (1 + sq1)^(1 + b))
end

"Thermal de Broglie factor (2π m_e k T / h²)^{3/2}, in m⁻³."
function _debroglie(T)
    h_P = 2π * _C.ħ_SI
    (2π * _C.m_e_SI * _C.k_B_SI * T / h_P^2)^1.5
end

"Convert an energy in eV to the equivalent temperature in K."
_eV_to_K(E) = E / _C.K_to_eV

# --- Baryon bookkeeping -----------------------------------------------------

"""
    n_H_of_z(c, z)

Number density of hydrogen nuclei (ionized or not), m⁻³.
"""
function n_H_of_z(c::Cosmology, z)
    ρ_b0 = Ω_b(c) * _C.ρ_crit_SI(c.h)          # kg/m³ today
    (1 - c.Yp) * ρ_b0 / _C.m_H_SI * (1 + z)^3
end

"""
    f_He(c)

Helium-to-hydrogen number ratio, n_He/n_H = Yp / (3.9715 (1 - Yp)).

Yp is a *mass* fraction, so converting it to a number ratio divides by the actual
He/H mass ratio -- 3.9715 (CLASS's `_not4_`), not 4. The α particle's binding
energy is 0.7% of its mass, and that 0.7% lands directly in the helium abundance
recombination sees.
"""
f_He(c::Cosmology) = c.Yp / (3.9715 * (1 - c.Yp))

"Hubble rate in 1/s (the recombination ODEs are naturally written in SI)."
H_SI(c::Cosmology, z) = 100 * c.h * 1e3 / _C.Mpc_SI * E_z(c, z)

# --- Saha equilibrium -------------------------------------------------------

"""
    saha_x_e(c, z)

Free-electron fraction x_e = n_e/n_H from full Saha equilibrium, including
singly and doubly ionized helium.

Three ionization balances are imposed at once,

    n_e n_HII  / n_HI   = S_H
    n_e n_HeII / n_HeI  = S_He1
    n_e n_HeIII/ n_HeII = S_He2

each of the form S = g (2π m_e kT/h²)^{3/2} e^{-E/kT}, with statistical weights
g = 1, 4, 1 respectively. They share the same n_e, so the system is coupled and
solved by root-finding on x_e rather than species by species.

Valid while the rates beat the expansion, i.e. above z ≈ 5000. Below that, real
recombination lags equilibrium badly (Saha would have recombination finish far
too early) and you want [`recombination`](@ref).
"""
function saha_state(c::Cosmology, z)
    T = c.Tcmb * (1 + z)
    nH = n_H_of_z(c, z)
    fHe = f_He(c)
    A = _debroglie(T)

    S_H = A * exp(-_C.T_of_wavenumber(_C.L_H_ion) / T)
    S_He1 = 4 * A * exp(-_C.T_of_wavenumber(_C.L_He1_ion) / T)
    S_He2 = 1 * A * exp(-_C.T_of_wavenumber(_C.L_He2_ion) / T)

    # Each balance gives an ionization state once n_e is known; they all share
    # the same n_e, so close the loop by root-finding on x_e.
    function states(xe)
        ne = xe * nH
        x_HII = S_H / (S_H + ne)                       # n_HII / n_H

        # The helium abundances are in the ratio
        #     HeI : HeII : HeIII  =  1 : S_He1/n_e : (S_He1/n_e)(S_He2/n_e)
        # but written that way the divisions blow up to Inf when n_e is tiny,
        # and Inf/Inf gives NaN. Multiplying the whole chain through by n_e²
        # is algebraically identical and has no small divisors anywhere.
        w_HeI = ne * ne
        w_HeII = S_He1 * ne
        w_HeIII = S_He1 * S_He2
        denom = w_HeI + w_HeII + w_HeIII
        x_HeII = w_HeII / denom
        x_HeIII = w_HeIII / denom
        (x_HII, x_HeII, x_HeIII)
    end

    function residual(log_xe)
        xe = exp(log_xe)
        x_HII, x_HeII, x_HeIII = states(xe)
        log(x_HII + fHe * (x_HeII + 2 * x_HeIII)) - log_xe
    end

    # The bracket has to reach very low: Saha x_e falls off exponentially once
    # recombination is done (~1e-16 by z = 500, and far below that later), so a
    # naive floor like 1e-13 fails to bracket the root at low redshift.
    xe_max = 1 + 2 * fHe
    xe = exp(find_zero(residual, (log(1e-250), log(xe_max)), Bisection()))
    x_HII, x_HeII, x_HeIII = states(xe)
    (x_H=x_HII, x_HeII=x_HeII, x_HeIII=x_HeIII, x_e=xe)
end

"Free-electron fraction from full Saha equilibrium (H, He I, He II)."
saha_x_e(c::Cosmology, z) = saha_state(c, z).x_e

"Hydrogen ionization fraction from full Saha equilibrium."
saha_x_H(c::Cosmology, z) = saha_state(c, z).x_H

# --- Reionization -----------------------------------------------------------

"""
    Reionization(; z_reio = 7.67, Δz = 0.5, ...)

The late-time reionization of the intergalactic medium by the first stars and
quasars, in the tanh parameterisation used by CAMB and CLASS.

Recombination leaves the universe neutral with a residual x_e ≈ 2e-4. It does not
stay that way: by z ≈ 6 the IGM is ionized again, and those free electrons scatter
CMB photons. About 5% of them do, which suppresses the CMB anisotropies by e^{-2τ}
-- roughly 10% in C_ℓ -- so this is not a detail that can be skipped before
computing a power spectrum.

    x_e(z) = (1 + f_He - x_rec) · [1 + tanh(y)]/2 + x_rec + (helium term)
    y      = [(1 + z_reio)^1.5 - (1+z)^1.5] / [1.5 (1+z_reio)^0.5 Δz]

The tanh is taken in (1+z)^1.5 rather than in z because that makes the resulting
optical depth nearly independent of the width Δz, so z_reio and τ are close to
interchangeable. The asymptote is 1 + f_He: hydrogen fully ionized and helium
singly ionized. Helium's *second* ionization needs harder photons and happens
later, around z ≈ 3.5, and is given its own tanh.

The default z_reio = 7.67 reproduces Planck 2018's τ = 0.0544. Use
[`z_reio_from_τ`](@ref) to solve for it given a target τ instead.
"""
Base.@kwdef struct Reionization
    enabled::Bool = true
    z_reio::Float64 = 7.67
    Δz::Float64 = 0.5
    helium2::Bool = true         # He II -> He III, the second helium reionization
    z_helium2::Float64 = 3.5
    Δz_helium2::Float64 = 0.4
end

"No reionization at all. Useful for isolating the recombination history."
const no_reionization = Reionization(enabled=false)

"""
    Injection(; p_ann = 0, Γ_dec = 0, f_dec = 0, f_eff = 1, f_halo = 0, z_halo = 30)

Exotic energy injection into the pre-reionization gas: dark-matter annihilation
and/or decay, in the same parameterisation as CLASS's injection module.

* `p_ann` — annihilation parameter ⟨σv⟩/m_DM in m³ s⁻¹ kg⁻¹. The injected power is
  dE/dV dt = ρ_cdm² c² p_ann, so it scales as (1+z)⁶ and matters earliest.
  (CLASS's `DM_annihilation_efficiency` is this divided by c².)
* `f_halo`, `z_halo` — the erfc((1+z)/(1+z_halo))/(1+z)³ boost once structure forms.
* `Γ_dec`, `f_dec` — decay rate (1/s) and the decaying mass fraction of the DM. The
  power is ρ_cdm c² f_dec Γ e^{−Γt}, scaling as (1+z)³ and mattering latest.
* `f_eff` — deposition efficiency: the fraction of injected energy deposited into
  the gas at the same redshift ("on the spot"). The energy that does deposit is
  split between heating, hydrogen ionization and Lyman-α excitation with the
  Chen–Kamionkowski fractions χ_heat = (1+2x_e)/3, χ_ion = χ_Lyα = (1−x_e)/3.

Injection deposits *into three different equations at once* — T_m directly, x_H
through both ionization and excitation — which is why it lives here and not in a
bolt-on: the CMB constrains p_ann precisely because extra x_e at z ~ 600 widens
the visibility function.
"""
Base.@kwdef struct Injection
    p_ann::Float64 = 0.0
    f_halo::Float64 = 0.0
    z_halo::Float64 = 30.0
    Γ_dec::Float64 = 0.0
    f_dec::Float64 = 0.0
    f_eff::Float64 = 1.0
end

const no_injection = Injection()

_has_injection(inj::Injection) = inj.p_ann > 0 || (inj.Γ_dec > 0 && inj.f_dec > 0)

"""
Injected power in J m⁻³ s⁻¹ at redshift `z`. `t_s` is cosmic time in seconds,
needed only for the decay exponential (pass 0.0 when Γ_dec = 0).
"""
function injection_rate(inj::Injection, c::Cosmology, z, t_s)
    ρ_cdm = Ω_c(c) * _C.ρ_crit_SI(c.h) * (1 + z)^3      # kg/m³
    rate = 0.0
    if inj.p_ann > 0
        boost = inj.f_halo > 0 ?
                inj.f_halo * erfc((1 + z) / (1 + inj.z_halo)) / (1 + z)^3 : 0.0
        rate += ρ_cdm^2 * _C.c_SI^2 * inj.p_ann * (1 + boost)
    end
    if inj.Γ_dec > 0 && inj.f_dec > 0
        rate += ρ_cdm * _C.c_SI^2 * inj.f_dec * inj.Γ_dec * exp(-inj.Γ_dec * t_s)
    end
    rate * inj.f_eff
end

"""
    reionization_x_e(reio, c, z, x_rec)

Free-electron fraction including reionization, smoothly matched onto the residual
`x_rec` left by recombination. At high z the tanh saturates at -1 and this returns
`x_rec` exactly, so it can be applied at every redshift without a switch.
"""
function reionization_x_e(reio::Reionization, c::Cosmology, z, x_rec)
    reio.enabled || return x_rec
    fHe = f_He(c)

    mid = (1 + reio.z_reio)^1.5
    Δ = 1.5 * sqrt(1 + reio.z_reio) * reio.Δz
    y = (mid - (1 + z)^1.5) / Δ
    xe = (1 + fHe - x_rec) * (tanh(y) + 1) / 2 + x_rec

    if reio.helium2 && z < reio.z_helium2 + 5 * reio.Δz_helium2
        y2 = (reio.z_helium2 - z) / reio.Δz_helium2
        xe += fHe * (tanh(y2) + 1) / 2
    end
    xe
end

# --- The RECFAST system -----------------------------------------------------

"""
    RecombinationSolution

Interpolated thermal history. Call the accessors rather than touching fields:
[`x_e`](@ref), [`T_matter`](@ref), [`n_e`](@ref).
"""
struct RecombinationSolution{C,I1,I2,I3,I4,I5}
    cosmo::C
    x_e_interp::I1     # log(1+z) -> x_e, including reionization
    T_m_interp::I2     # log(1+z) -> T_matter [K]
    τ_interp::I3       # log(1+z) -> Thomson optical depth back to z (full)
    τd_interp::I4      # log(1+z) -> baryon drag optical depth
    τrec_interp::I5    # log(1+z) -> optical depth from recombination alone
    z_min::Float64
    z_max::Float64
end

"""
    recombination(c; z_start = 5000, z_end = 0, ...)

Integrate the recombination history and return a [`RecombinationSolution`](@ref).

The state vector is `(x_H, x_HeII, log T_m)`:
- `x_H`    = n_HII / n_H,  hydrogen ionization fraction
- `x_HeII` = n_HeII / n_He, singly-ionized helium fraction
- `T_m`    = matter temperature in K

with x_e = x_H + f_He · x_HeII. Doubly-ionized helium is handled by Saha above
`z_start` (it recombines by z ≈ 5000 and never leaves equilibrium).

# Accuracy

This is RECFAST 1.5 (Wong, Moss & Scott 2008; Scott & Moss 2009), whose
correction terms are calibrated against the full multi-level-atom codes --
HyRec (Ali-Haïmoud & Hirata 2010) and CosmoRec -- and against the Switzer &
Hirata (2008) helium papers. Beyond the textbook Peebles treatment it carries:

  * `fudge_H = 1.125`, retuned in 2012 to match HyRec/CosmoRec. A plain
    three-level hydrogen atom recombines measurably too slowly; this absorbs
    the leading effect of the ~100 bound levels a real atom has.
  * `gaussian_H`, a double-Gaussian correction to the Lyman-α escape rate K.
    This is how the two-photon and radiative-transfer physics of the
    multi-level codes is folded into a two-level equation.
  * Helium via **Sobolev escape probability** rather than a bare λ³/8πH, so the
    optically thin limit is right and not just the thick one.
  * **H I continuum opacity** during helium recombination (Switzer & Hirata
    2008 I; Kholupenko et al. 2007 fit). Neutral hydrogen is already present
    while helium is still recombining, and it absorbs the He I line photons
    that would otherwise have re-ionized helium. This *accelerates* He I
    recombination substantially and is the single largest helium correction.
  * The **2³P₁ → 1¹S₀ intercombination line** (`triplets`), a second, forbidden
    recombination pathway. Its Einstein coefficient is 177 s⁻¹ against the
    singlet's 1.8e9 s⁻¹ -- ten million times weaker -- but the singlet channel
    is so throttled by its own resonant opacity that the forbidden line still
    competes.

Set `fudge_H = 1.14, gaussian_H = false, sobolev_He = false, triplets = false`
to recover the original RECFAST 1.0 of Seager, Sasselov & Scott.

The system is stiff: the reaction rates dwarf H at early times, which is
precisely why Saha holds there. It needs an implicit solver.
"""
function recombination(c::Cosmology;
    z_start=5000.0, z_end=0.0,
    fudge_H=1.125,          # RECFAST 1.5 value, matched to HyRec/CosmoRec
    fudge_He=0.86,          # exponent qb in the H-continuum-opacity fit
    gaussian_H=true,        # double-Gaussian correction to K_H
    sobolev_He=true,        # Sobolev escape + H continuum opacity for He I
    triplets=true,          # 2³P₁ -> 1¹S₀ intercombination channel
    reionization=Reionization(),
    injection=no_injection,
    hydrogen=:recfast,      # :recfast, or :hyrec for the HyRec-2020 SWIFT model
    reltol=1e-10, abstol=1e-12)

    hydrogen in (:recfast, :hyrec) ||
        throw(ArgumentError("hydrogen must be :recfast or :hyrec"))

    fHe = f_He(c)
    Tr(z) = c.Tcmb * (1 + z)

    # SWIFT correction-function sensitivities: offsets of this cosmology from the
    # fit's fiducial, scaled as in HyRec (wrap_hyrec.c / rec_swift_hyrec_dxHIIdlna).
    T0fid3 = (2.7255 / c.Tcmb)^3
    dω_cb = ((Ω_b(c) + Ω_c(c)) * c.h^2 - 0.14175) * T0fid3
    dω_bH = (Ω_b(c) * c.h^2 * (1 - c.Yp) - 0.02242 * (1 - 0.246738546372)) * T0fid3
    dNeff = _N_eff_of(c) - 3.046
    # HyRec's hydrogen model needs its rate tables (0.004 < T_R < 0.4 eV, i.e.
    # z ≲ 1700). Above that hydrogen is pinned to Saha by enormous rates and the
    # RECFAST equation is exact; below the table floor (z ≈ 16) recombination has
    # been frozen for ages. Use the tabulated model exactly where it exists.
    swift_lo, swift_hi = 0.0042, 0.395

    inj_on = _has_injection(injection)
    # Cosmic time, for the decay exponential only. age() runs a quadrature, and a
    # quadrature inside a stiff ODE right-hand side is the recurring mistake of
    # this codebase -- tabulate it once.
    t_of_z = if injection.Γ_dec > 0 && injection.f_dec > 0
        lzs = collect(range(log(1.0), log(1 + z_start); length=200))
        ts = [age(c, expm1(lz)) * _C.Gyr_SI for lz in lzs]
        itp = linear_interpolation(lzs, ts; extrapolation_bc=Line())
        z -> itp(log(1 + z))
    else
        z -> 0.0
    end
    E_H_ion_J = _C.k_B_SI * _C.T_of_wavenumber(_C.L_H_ion)      # 13.6 eV in J
    E_H_lya_J = _C.k_B_SI * _C.T_of_wavenumber(_C.L_H_alpha)    # 10.2 eV in J
    Tw = _C.T_of_wavenumber

    # Level energies as temperatures.
    K_ion_H = Tw(_C.L_H_ion)                       # 157 807 K
    K_Lyα_H = Tw(_C.L_H_alpha)                     # 118 376 K
    K_2s_H = K_ion_H - K_Lyα_H                    #  39 431 K (n=2 binding)
    K_ion_He1 = Tw(_C.L_He1_ion)
    K_2s_He1 = Tw(_C.L_He1_ion - _C.L_He_2s)        # He I 2¹S binding
    K_2s1s_He1 = Tw(_C.L_He_2s)                      # 2¹S -> 1¹S energy
    K_2p2s_He1 = Tw(_C.L_He_2p - _C.L_He_2s)         # singlet 2p-2s splitting
    K_2Pt2St = Tw(_C.L_He_2Pt - _C.L_He_2St)       # triplet 2³P-2³S splitting
    K_2St_ion = Tw(_C.L_He2St_ion)                  # ionization from 2³S
    K_2St1s = Tw(_C.L_He_2St)                     # 2³S -> 1¹S energy

    λ_Lyα = 1 / _C.L_H_alpha                        # m
    CK_H = λ_Lyα^3 / (8π)
    CK_He = (1 / _C.L_He_2p)^3 / (8π)

    m_He = _C.m_He_SI

    function rhs!(du, u, p, z)
        x_H = clamp(u[1], 0.0, 1.0)
        x_HeII = clamp(u[2], 0.0, 1.0)
        T_m = exp(u[3])

        T_r = Tr(z)
        nH = n_H_of_z(c, z)
        n_He_tot = fHe * nH
        Hz = H_SI(c, z)
        x_e = x_H + fHe * x_HeII
        Adb = _debroglie(T_m)

        n_1s = (1 - x_H) * nH             # neutral hydrogen
        n_He1 = (1 - x_HeII) * n_He_tot    # neutral helium

        # ---- Hydrogen ------------------------------------------------------
        aH = α_H(T_m, fudge_H)
        βH = aH * Adb * exp(-K_2s_H / T_m)          # photoionization from n=2
        β_ground_H = aH * Adb * exp(-K_ion_H / T_m)   # detailed balance vs 1s

        # K is the inverse rate at which Lyman-α photons redshift out of the
        # line. The Gaussian correction stands in for the multi-level-atom and
        # two-photon physics that a three-level model cannot represent.
        K_H = CK_H / Hz
        if gaussian_H
            lz = log(1 + z)
            K_H *= 1 - 0.14 * exp(-((lz - 7.28) / 0.18)^2) +
                   0.079 * exp(-((lz - 6.73) / 0.33)^2)
        end

        TR_eV = _HY_kB * T_r
        if hydrogen === :hyrec && swift_lo < TR_eV < swift_hi
            dxlna, w_lya = hyrec_dxHII_dlna(x_e, x_H, nH * 1e-6, Hz,
                _HY_kB * T_m, TR_eV, dω_cb, dω_bH, dNeff)
            dxH = -dxlna / (1 + z)
        else
            C_H = (1 + K_H * _C.Λ_2s1s_H * n_1s) /
                  (1 + K_H * (_C.Λ_2s1s_H + βH) * n_1s)
            dxH = C_H * (x_e * x_H * nH * aH - β_ground_H * (1 - x_H)) / (Hz * (1 + z))
            w_lya = 1 - C_H
        end

        # ---- Exotic energy injection ----------------------------------------
        # Direct ionizations spend 13.6 eV each; Lyman-α excitations ionize only
        # the fraction w_lya that fails to decay back to the ground state. Same
        # structure as CLASS's wrap_recfast.c / wrap_hyrec.c, so the codes are
        # comparable term by term.
        if inj_on
            x = clamp(x_e, 0.0, 1.0)
            dE = injection_rate(injection, c, z, t_of_z(z))     # J/m³/s
            χ_ion = (1 - x) / 3
            χ_lya = (1 - x) / 3
            χ_heat = min((1 + 2x) / 3, 1.0)
            dxH -= (χ_ion * dE / E_H_ion_J + χ_lya * dE * w_lya / E_H_lya_J) /
                   (nH * Hz * (1 + z))
            du_inj = χ_heat * dE / (1.5 * _C.k_B_SI * nH * (1 + fHe + x_e)) /
                     (Hz * (1 + z))
        else
            du_inj = 0.0
        end

        # ---- Helium: singlet channel (2¹P₁ -> 1¹S₀) -------------------------
        # The factor 4 = 2·g_HeII/g_HeI is the Saha degeneracy ratio detailed
        # balance ties β to. Hydrogen's happens to be 1; helium's is not.
        aHe = α_He(T_m)
        βHe = 4 * aHe * Adb * exp(-K_2s_He1 / T_m)
        β_ground_He = βHe * exp(-K_2s1s_He1 / T_m)

        # 2s and 2p are split in helium (unlike hydrogen, where they are
        # degenerate), so the escape term must be reweighted between them.
        He_Boltz = exp(min(K_2p2s_He1 / T_m, 680.0))   # guard the overflow

        if !sobolev_He || x_HeII < 5e-9 || x_HeII > 0.98
            K_He = CK_He / Hz
        else
            # Sobolev escape probability: correct in *both* the optically thick
            # limit (where it reduces to λ³/8πH) and the thin limit.
            τ_He_s = _C.A2P_s * CK_He * 3 * n_He1 / Hz
            p_He_s = (1 - exp(-τ_He_s)) / τ_He_s
            A_esc = _C.A2P_s * p_He_s

            if x_H < 0.9999999
                # Continuum opacity from neutral hydrogen. Once any H I exists
                # it soaks up the He I resonance photons that would otherwise
                # re-ionize helium, opening an escape route the pure-helium
                # calculation does not have. This is the dominant He correction.
                doppler = _C.c_SI * _C.L_He_2p *
                          sqrt(2 * _C.k_B_SI * T_m / (m_He * _C.c_SI^2))
                γ_2Ps = 3 * _C.A2P_s * fHe * (1 - x_HeII) * _C.c_SI^2 /
                        (sqrt(π) * _C.σ_He_2Ps * 8π * doppler * (1 - x_H)) /
                        (_C.c_SI * _C.L_He_2p)^2
                A_Hcon = _C.A2P_s / (1 + 0.36 * γ_2Ps^fudge_He)   # pb = 0.36, KIV 2007
                A_esc += A_Hcon
            end
            K_He = 1 / (A_esc * 3 * n_He1)
        end

        n_He1_eff = n_He1 * He_Boltz
        C_He = (1 + K_He * _C.Λ_2s1s_He1 * n_He1_eff) /
               (1 + K_He * (_C.Λ_2s1s_He1 + βHe) * n_He1_eff)
        dxHe = C_He * (x_e * x_HeII * nH * aHe - β_ground_He * (1 - x_HeII)) /
               (Hz * (1 + z))

        # ---- Helium: triplet channel (2³P₁ -> 1¹S₀) -------------------------
        if triplets && 5e-9 < x_HeII < 0.98
            sq0 = sqrt(T_m / 10.0^0.477121)
            sq1 = sqrt(T_m / 10.0^5.114)
            α_trip = _C.a_trip / (sq0 * (1 + sq0)^(1 - _C.b_trip) *
                                  (1 + sq1)^(1 + _C.b_trip))
            β_trip = α_trip * Adb * exp(-K_2St_ion / T_m) * (4 / 3)

            τ_He_t = _C.A2P_t * n_He1 * 3 / (8π * Hz * _C.L_He_2Pt^3)
            p_He_t = (1 - exp(-τ_He_t)) / τ_He_t
            A_esc_t = _C.A2P_t * p_He_t

            if x_H < 0.99999
                doppler_t = _C.c_SI * _C.L_He_2Pt *
                            sqrt(2 * _C.k_B_SI * T_m / (m_He * _C.c_SI^2))
                γ_2Pt = 3 * _C.A2P_t * fHe * (1 - x_HeII) * _C.c_SI^2 /
                        (sqrt(π) * _C.σ_He_2Pt * 8π * doppler_t * (1 - x_H)) /
                        (_C.c_SI * _C.L_He_2Pt)^2
                A_esc_t += _C.A2P_t / (1 + 0.66 * γ_2Pt^0.9) / 3   # KIV 2007
            end

            Cf_t = A_esc_t * exp(-K_2Pt2St / T_m)
            Cf_t = Cf_t / (β_trip + Cf_t)
            dxHe += (x_e * x_HeII * nH * α_trip -
                     (1 - x_HeII) * 3 * β_trip * exp(-K_2St1s / T_m)) *
                    Cf_t / (Hz * (1 + z))
        end

        # ---- Matter temperature --------------------------------------------
        # Compton scattering off the CMB drags T_m toward T_r; adiabatic
        # expansion cools it as (1+z)². The competition thermally decouples the
        # gas a little after it decouples chemically. Evolved in log to keep it
        # positive -- by z = 0 it has fallen four orders of magnitude and a bare
        # T_m would step straight through zero and poison the rate fits.
        Γ_C = 8 * _C.σ_T_SI * _C.a_rad_SI * T_r^4 /
              (3 * _C.m_e_SI * _C.c_SI) * x_e / (1 + fHe + x_e)
        dTm = Γ_C * (T_m - T_r) / (Hz * (1 + z)) + 2 * T_m / (1 + z) - du_inj

        du[1] = dxH
        du[2] = dxHe
        du[3] = dTm / T_m
        nothing
    end

    st = saha_state(c, z_start)
    u0 = [st.x_H, st.x_HeII, log(Tr(z_start))]

    prob = ODEProblem(rhs!, u0, (z_start, z_end))
    sol = solve(prob, Rodas5P(); reltol, abstol)

    # Resample onto a log(1+z) grid: x_e varies over five decades and a linear
    # z grid wastes almost all its points where nothing happens.
    lz = range(log(1 + z_end), log(1 + z_start); length=4000)
    zs = expm1.(lz)
    xes = similar(zs)
    tms = similar(zs)
    xes_rec = similar(zs)
    for (i, z) in enumerate(zs)
        uu = sol(z)
        x_rec = clamp(uu[1], 0.0, 1.0) + fHe * clamp(uu[2], 0.0, 1.0)
        # Reionization is layered on top of the recombination history rather than
        # evolved: the astrophysics driving it -- when the first stars light up,
        # how their ionizing photons escape -- is not something the atomic-physics
        # ODE above knows about. It is parameterised, not derived.
        xes_rec[i] = x_rec
        xes[i] = reionization_x_e(reionization, c, z, x_rec)
        tms[i] = exp(uu[3])
    end

    # Accumulate the two optical depths once, on the same grid, by trapezoid in
    # log(1+z). Recomputing τ from scratch inside every visibility() call is a
    # nested quadrature and turns any sweep over redshift into an O(N²) crawl --
    # and the perturbation layer evaluates τ and g at every step of every mode.
    #
    #   dτ/dz      = n_e σ_T / [(1+z) H]        (Thomson)
    #   dτ_drag/dz = n_e σ_T / [(1+z) H R]      (baryon drag: R scatterings
    #                                            are needed to move one baryon)
    # and with x = log(1+z), dz = (1+z) dx, so the (1+z) cancels.
    dτ = similar(zs)
    dτd = similar(zs)
    dτr = similar(zs)
    for (i, z) in enumerate(zs)
        a = scale_factor(z)
        pre = n_H_of_z(c, z) * _C.Mpc_SI^3 * _C.σ_T_Mpc2 / H_Mpc(c, a)  # dτ/dlog(1+z) per x_e
        dτ[i] = xes[i] * pre
        dτd[i] = xes[i] * pre / R_baryon(c, a)
        dτr[i] = xes_rec[i] * pre
    end

    τs = similar(zs)
    τds = similar(zs)
    τrs = similar(zs)
    τs[1] = τds[1] = τrs[1] = 0.0
    for i in 2:length(zs)
        h = lz[i] - lz[i-1]
        τs[i] = τs[i-1] + 0.5 * h * (dτ[i] + dτ[i-1])
        τds[i] = τds[i-1] + 0.5 * h * (dτd[i] + dτd[i-1])
        τrs[i] = τrs[i-1] + 0.5 * h * (dτr[i] + dτr[i-1])
    end

    RecombinationSolution(
        c,
        linear_interpolation(lz, xes; extrapolation_bc=Line()),
        linear_interpolation(lz, tms; extrapolation_bc=Line()),
        linear_interpolation(lz, τs; extrapolation_bc=Line()),
        linear_interpolation(lz, τds; extrapolation_bc=Line()),
        linear_interpolation(lz, τrs; extrapolation_bc=Line()),
        float(z_end), float(z_start),
    )
end

"""
    x_e(r, z)

Free-electron fraction x_e = n_e/n_H.

Above the ODE start redshift this falls back to full Saha equilibrium, which is
exact there: the reaction rates beat the expansion by many orders of magnitude,
so the gas has no choice but to sit in equilibrium. Extrapolating the
interpolant instead would be both slower and wrong.
"""
function x_e(r::RecombinationSolution, z)
    z > r.z_max && return saha_x_e(r.cosmo, z)
    r.x_e_interp(log(1 + z))
end

"""
    T_matter(r, z)

Baryon temperature in K. Locked to T_cmb while Compton scattering is efficient,
then falls as (1+z)² once the gas decouples thermally.
"""
function T_matter(r::RecombinationSolution, z)
    z > r.z_max && return r.cosmo.Tcmb * (1 + z)
    r.T_m_interp(log(1 + z))
end

"Free-electron number density in m⁻³."
n_e(r::RecombinationSolution, z) = x_e(r, z) * n_H_of_z(r.cosmo, z)

# --- Optical depth and visibility -------------------------------------------

"""
    τ_dot(r, a)

dτ/dη = -a n_e σ_T, in 1/Mpc. Negative because τ is measured backward from the
observer. This is the quantity the Boltzmann hierarchy actually consumes: it is
the photon scattering rate, and its magnitude relative to k sets whether the
photons are tightly coupled to the baryons or free-streaming.
"""
function τ_dot(r::RecombinationSolution, a)
    z = redshift(a)
    ne_Mpc3 = n_e(r, z) * _C.Mpc_SI^3          # m⁻³ -> Mpc⁻³
    -a * ne_Mpc3 * _C.σ_T_Mpc2
end

"""
    optical_depth(r, z)

Thomson optical depth back to redshift z,

    τ(z) = ∫₀^z  n_e σ_T c  dz'/[(1+z') H(z')]

Dimensionless; τ(z_*) = 1 defines the surface of last scattering. Precomputed
at construction, so this is an interpolation lookup.
"""
optical_depth(r::RecombinationSolution, z) = r.τ_interp(log(1 + z))

"""
    dτ_dz(r, z)

-dτ/dz, the differential scattering probability per unit redshift. Positive.
"""
function dτ_dz(r::RecombinationSolution, z)
    c = r.cosmo
    ne_Mpc3 = n_e(r, z) * _C.Mpc_SI^3
    ne_Mpc3 * _C.σ_T_Mpc2 / ((1 + z) * H_Mpc(c, scale_factor(z)))
end

"""
    visibility(r, z)

The visibility function g(z) = (-dτ/dz)·e^{-τ}, normalised so ∫g dz = 1.

This is the probability density that a CMB photon we see today last scattered
at redshift z. Its peak *is* the last-scattering surface, and its finite width
is why that surface has a thickness rather than being a sharp shell -- which in
turn is what damps the high-ℓ acoustic peaks.
"""
visibility(r::RecombinationSolution, z) = dτ_dz(r, z) * exp(-optical_depth(r, z))

"""
    τ_reio(r)

Optical depth to reionization: the Thomson depth accumulated between us and
z = 50, which is well above the reionization transition but far below
recombination.

This is one of the six parameters of ΛCDM, and it is almost entirely degenerate
with the amplitude A_s in the CMB -- both scale C_ℓ, and only the large-scale
polarization breaks the degeneracy. Planck 2018: τ = 0.0544 ± 0.0073.
"""
τ_reio(r::RecombinationSolution) = optical_depth(r, 50.0)

"""
    z_reio_from_τ(c, τ_target; kws...)

Invert τ(z_reio) to find the reionization redshift giving a target optical depth.
CAMB and CLASS both offer this, because τ is what the CMB measures and z_reio is
what the model is parameterised by.
"""
function z_reio_from_τ(c::Cosmology, τ_target; kws...)
    f(zr) = τ_reio(recombination(c; reionization=Reionization(z_reio=zr), kws...)) - τ_target
    find_zero(f, (2.0, 25.0), Bisection(); xatol=1e-3)
end

"""
    recombination_optical_depth(r, z)

Optical depth back to `z` counting only the electrons left free by recombination,
excluding those liberated again by reionization.
"""
recombination_optical_depth(r::RecombinationSolution, z) = r.τrec_interp(log(1 + z))

"""
    z_star(r)

Redshift of last scattering, defined by τ(z_*) = 1.

Deliberately measured against the *recombination* optical depth, not the total.
Reionization contributes τ ≈ 0.054 of its own, and counting it would push the
τ = 1 crossing down to z ≈ 1085 -- but those low-redshift scatterings are not what
"last scattering" means. The surface that stamps the acoustic peaks onto the CMB
is the recombination one, and this is the convention Planck's 1089.92 uses.
"""
function z_star(r::RecombinationSolution)
    find_zero(z -> recombination_optical_depth(r, z) - 1, (500.0, 2000.0),
        Bisection(); xatol=1e-4)
end

# --- Sound horizon and the drag epoch ---------------------------------------

"""
    R_baryon(c, a)

Baryon-to-photon momentum density ratio, R = 3ρ_b/(4ρ_γ). It is what makes the
photon-baryon fluid heavy: the sound speed is c/√(3(1+R)), so R suppresses the
oscillation frequency and is directly responsible for the alternating peak
heights in the CMB.
"""
R_baryon(c::Cosmology, a) = 3 * Ω_b(c) * a / (4 * Ω_γ(c))

"Sound speed of the photon-baryon fluid, in units of c."
sound_speed(c::Cosmology, a) = 1 / sqrt(3 * (1 + R_baryon(c, a)))

"""
    sound_horizon(c, z)

Comoving distance a sound wave has travelled by redshift z, in Mpc:

    r_s(z) = ∫_z^∞  c_s(z') dz' / H(z')  =  ∫₀^a  c_s da'/(a'² H)

Evaluated at the drag epoch this is *the* standard ruler: it sets the BAO scale
and the CMB acoustic peak spacing.
"""
function sound_horizon(c::Cosmology, z; rtol=1e-10)
    a_end = scale_factor(z)
    integrand(a) = sound_speed(c, a) / (a^2 * E(c, a))
    hubble_distance(c) * quadgk(integrand, 0.0, a_end; rtol)[1]
end

"""
    z_drag(r)

The baryon drag epoch, where the *baryons* stop being dragged by the photons.
Defined by τ_drag(z_d) = 1 with

    τ_drag(z) = ∫₀^z  (dτ/dz') / R  dz'

The 1/R weighting is the whole point, and it is why z_drag ≈ 1060 differs from
z_* ≈ 1090: what decouples a photon is one scattering, but to change a baryon's
momentum takes R scatterings, since the baryons carry R times the momentum. In
our universe R ≈ 0.6 at recombination, so the baryons are released slightly
*after* the photons.
"""
drag_optical_depth(r::RecombinationSolution, z) = r.τd_interp(log(1 + z))

function z_drag(r::RecombinationSolution)
    find_zero(z -> drag_optical_depth(r, z) - 1, (500.0, 2000.0), Bisection(); xatol=1e-4)
end

"Sound horizon at the drag epoch, in Mpc -- the BAO standard ruler."
r_drag(r::RecombinationSolution) = sound_horizon(r.cosmo, z_drag(r))

"Sound horizon at last scattering, in Mpc."
r_star(r::RecombinationSolution) = sound_horizon(r.cosmo, z_star(r))

"""
    θ_star(r)

Angular size of the sound horizon at last scattering, θ_* = r_s(z_*)/D_M(z_*).
The single best-measured number in cosmology (Planck: 100θ_* = 1.04109 ± 0.00030).
"""
function θ_star(r::RecombinationSolution)
    zs = z_star(r)
    sound_horizon(r.cosmo, zs) / transverse_comoving_distance(r.cosmo, zs)
end
