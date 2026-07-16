"""
N_eff from first principles: finite-temperature QED and neutrino decoupling.

N_eff = 3.044 is the last hardcoded number in the code, and it is not a free parameter --
it is a *prediction*, and one that took the field twenty years to nail down (Mangano 2005:
3.046; Bennett et al. 2021: 3.0440). It is bigger than 3 for two reasons, and both are
computable:

  * **Neutrinos do not decouple instantaneously.** When e± annihilate, most of their
    entropy goes to the photons -- but the neutrinos have not *quite* finished decoupling,
    so a little of it leaks into them instead. The neutrinos come out slightly hotter than
    (4/11)^{1/3} T_γ, and since N_eff is *defined* against that ratio, the excess shows up
    as N_eff > 3.

  * **The QED plasma is not an ideal gas.** Photons and electrons interact, so the pressure
    of the e±γ plasma is not the free-field one. That changes how much entropy the
    annihilation has to dump, and therefore the final temperature ratio.

Reference codes do not compute this. PRIMAT reads the neutrino history from precomputed
NEVO tables; PRyMordial ships `nue_ann.txt`/`nue_scatt.txt` from NUDEC_BSM and interpolates.
Here both pieces are computed.

# The QED plasma pressure

This part is closed-form, and is thermal field theory in the strict sense. With x = m_e/T,

    δP = δP_a  [O(e²)]  +  δP_e3  [O(e³)]  +  δP_b  [O(e⁴)]

    δP_a(T)  = (α/π) T⁴ [ -(2/3) I₀₁(x) - (2/π²) I₀₁(x)² ]
    δP_e3(T) = α^{3/2} (4/3) √(2π) T⁴ [ (I₀₁(x) + I₂₋₁(x)) / π² ]^{3/2}
    δP_b(T)  = (α/π³) T⁴ ∫∫ dp₁dp₂  x² p₁p₂/(e₁e₂) · ln|(p₁+p₂)/(p₁-p₂)| / [(e^{e₁}+1)(e^{e₂}+1)]

with e_i = √(p_i² + x²) and the two Fermi integrals

    I₀₁(x)  = ∫₀^∞ dp  p² / [ √(p²+x²) (e^{√(p²+x²)} + 1) ]
    I₂₋₁(x) = ∫₀^∞ dp  √(p²+x²) / (e^{√(p²+x²)} + 1)

δP_a is the one-loop Frenkel-Galitskii-Migdal term and is *negative* -- the interaction
lowers the pressure below the ideal gas. δP_e3 is the ring/plasmon resummation
(Blaizot-Zinn-Justin): it is O(e³), not O(e⁴), because the soft-photon loop needs Debye
resummation, and a naive perturbative expansion in α misses it entirely. It is positive and
about ten times smaller. δP_b is the two-loop exchange term, smaller again by two orders.

The energy density follows from the pressure by the standard thermodynamic identity, with
no extra input:

    δρ = T ∂(δP)/∂T - δP
"""

using QuadGK: quadgk, gauss
using OrdinaryDiffEq: ODEProblem, solve, Tsit5, ContinuousCallback, terminate!
using Interpolations: linear_interpolation, Line
using DelimitedFiles: readdlm

# --- the two Fermi integrals ------------------------------------------------

"I₀₁(x) = ∫₀^∞ dp p²/[√(p²+x²)(e^{√(p²+x²)}+1)]"
function _I01(x)
    quadgk(p -> begin
            e = sqrt(p^2 + x^2)
            e > 300 ? 0.0 : p^2 / (e * (exp(e) + 1))
        end, 0.0, 60 + 5x; rtol=1e-11)[1]
end

"I₂₋₁(x) = ∫₀^∞ dp √(p²+x²)/(e^{√(p²+x²)}+1)"
function _I2m1(x)
    quadgk(p -> begin
            e = sqrt(p^2 + x^2)
            e > 300 ? 0.0 : e / (exp(e) + 1)
        end, 0.0, 60 + 5x; rtol=1e-11)[1]
end

# --- the pressure corrections -----------------------------------------------

"O(e²) one-loop interaction pressure (Frenkel, Galitskii & Migdal). Negative."
function δP_e2(T)
    x = _ME_MEV / T
    I = _I01(x)
    (_ALPHA_FS / π) * T^4 * (-(2 / 3) * I - (2 / π^2) * I^2)
end

"""
O(e³) ring / plasmon pressure (Blaizot & Zinn-Justin).

Order e³, not e⁴: the soft-photon loop is infrared-divergent order by order and has to be
Debye-resummed, which promotes it to half-integer order in α. A straight perturbative
expansion in α never produces this term.
"""
function δP_e3(T)
    x = _ME_MEV / T
    m_D2 = (_I01(x) + _I2m1(x)) / π^2      # Debye mass² in units of T², up to 2α
    _ALPHA_FS^1.5 * (4 / 3) * sqrt(2π) * T^4 * m_D2^1.5
end

"O(e⁴) two-loop exchange pressure. Small (~1% of δP_a); included for completeness."
function δP_e4(T; rtol=1e-7)
    x = _ME_MEV / T
    f(p1, p2) = begin
        (p1 <= 0 || p2 <= 0 || abs(p1 - p2) < 1e-12) && return 0.0
        e1 = sqrt(p1^2 + x^2)
        e2 = sqrt(p2^2 + x^2)
        (e1 > 300 || e2 > 300) && return 0.0
        (_ALPHA_FS / π^3) * x^2 * p1 * p2 / (e1 * e2) *
        log(abs((p1 + p2) / (p1 - p2))) / ((exp(e1) + 1) * (exp(e2) + 1))
    end
    hi = 40 + 5x
    inner(p1) = quadgk(p2 -> f(p1, p2), 0.0, p1 - 1e-9, p1 + 1e-9, hi; rtol)[1]
    T^4 * quadgk(inner, 0.0, hi; rtol)[1]
end

"""
    qed_pressure(T; order = 3)

Total QED interaction pressure of the e±γ plasma, in MeV⁴, for T in MeV.
`order = 2`, `3` or `4` selects how far the expansion in e is carried.
"""
function qed_pressure(T; order=3)
    p = δP_e2(T)
    order >= 3 && (p += δP_e3(T))
    order >= 4 && (p += δP_e4(T))
    p
end

"""
    qed_energy_density(T; order = 3)

δρ = T ∂(δP)/∂T - δP, the thermodynamic identity. Nothing new goes in here: given the
pressure, the energy density is fixed.
"""
function qed_energy_density(T; order=3, h=1e-5)
    dP = (qed_pressure(T * (1 + h); order) - qed_pressure(T * (1 - h); order)) / (2 * h * T)
    T * dP - qed_pressure(T; order)
end

# =============================================================================
# Neutrino decoupling
# =============================================================================
#
# Three coupled temperatures -- T_γ (the e±γ plasma), T_νe, and T_νμ (shared by ν_μ and
# ν_τ, which have identical couplings) -- exchanging energy through the weak interaction
# while the universe expands. Energy conservation for each sector:
#
#   dρ_plasma/dt + 3H(ρ_plasma + P_plasma) = -(δρ_νe + 2 δρ_νμ)
#   dρ_νe/dt     + 4H ρ_νe                 = +δρ_νe
#   dρ_νμ/dt     + 4H ρ_νμ                 = +δρ_νμ
#
# (neutrinos are relativistic, so P = ρ/3 and 3H(ρ+P) = 4Hρ). The plasma carries the QED
# interaction pressure computed above, which is the point: without it the entropy the
# annihilation has to dump is wrong, and so is the final temperature ratio.
#
# If the weak rates were infinite the neutrinos would track the photons exactly and
# decouple instantaneously, giving T_ν/T_γ = (4/11)^{1/3} and N_eff = 3. They are not
# infinite. The δρ terms keep leaking a little e± annihilation energy into the neutrinos
# after they have mostly decoupled, and the leak is what makes N_eff > 3.

# Electroweak couplings, at tree level from α, G_F and m_Z -- not typed in.
const _GF = 1.1663787e-11              # MeV^-2, PDG
const _MZ = 91.1876e3                  # MeV
const _SW2 = 0.5 * (1 - sqrt(1 - 2 * sqrt(2) * π * _ALPHA_FS / (_GF * _MZ^2)))
const _gL_e = 0.5 + _SW2               # ν_e: charged + neutral current
const _gR_e = _SW2
const _gL_mu = -0.5 + _SW2             # ν_μ, ν_τ: neutral current only
const _gR_mu = _SW2
const _GN = 6.70883e-45                # MeV^-2

# Pauli-blocking corrections to the Maxwell-Boltzmann collision integrals
# (Escudero, JCAP 05 (2020) 048). The MB result overestimates the rates because it
# ignores the fact that the final states are already partly occupied.
const _F_ANN_FD = 0.884
const _F_SCAT_FD = 0.829

"Energy density of one neutrino flavour (ν + ν̄) at temperature T, in MeV⁴."
_ρ_nu(T) = 2 * (7 / 8) * (π^2 / 30) * T^4

"Plasma energy density and pressure: photons + e± + the QED interaction terms."
_ρ_plasma(T; order=3) = (π^2 / 15) * T^4 + _ρ_epm(T) + qed_energy_density(T; order)
_P_plasma(T; order=3) = (π^2 / 45) * T^4 + _P_epm(T) + qed_pressure(T; order)

"""
Weak energy-transfer rate between a bath at T1 and one at T2, in the massless-electron
limit (Escudero 2020):

    f(T1,T2) = 32·f_ann·(T1⁹ - T2⁹)·F_ann(T1) + 56·f_scat·F_scat(T1)·T1⁴T2⁴·(T1 - T2)

The two terms are the two ways energy moves: pair annihilation (e⁺e⁻ ↔ νν̄), which
carries the T⁹ scaling, and elastic scattering (νe± → νe±), which is proportional to the
temperature *difference* and so shuts off the moment the two baths equalise.

`F_ann` and `F_scat` are the finite-electron-mass corrections: the ratio of the rate with
a massive electron to the massless one. They tend to 1 for T ≫ m_e and to 0 once the
electrons are gone.
"""
function _f_transfer(T1, T2, F_ann, F_scat)
    32 * _F_ANN_FD * (T1^9 - T2^9) * F_ann +
    56 * _F_SCAT_FD * F_scat * T1^4 * T2^4 * (T1 - T2)
end

# ν_μ <-> ν_e transfer: neutrinos are massless, so no F factors.
_f_nunu(T1, T2) = 32 * _F_ANN_FD * (T1^9 - T2^9) + 56 * _F_SCAT_FD * T1^4 * T2^4 * (T1 - T2)

# --- finite-electron-mass corrections to the collision rates ------------------
#
# F_ann and F_scat are the ratio of the energy-transfer rate with a massive electron to
# the same rate with a massless one. They are what PRyMordial ships as `nue_ann.txt` and
# NUDEC_BSM computed; here they come out of the matrix elements.
#
# FeynCalc, for e⁻(p₁)e⁺(p₂) → ν(p₃)ν̄(p₄) in the effective four-fermion theory:
#
#   |M|² = 128 G_F² [ g_L²(p₁·p₄)(p₂·p₃) + g_R²(p₁·p₃)(p₂·p₄) + g_L g_R m_e² (p₃·p₄) ]
#
# and for ν(k₁)e(p₁) → ν(k₂)e(p₂) the same structure with the mass term *negative*:
#
#   |M|² = 128 G_F² [ g_L²(k₁·p₁)(k₂·p₂) + g_R²(k₁·p₂)(k₂·p₁) - g_L g_R m_e² (k₁·k₂) ]
#
# The sign flip on the m_e² term between the two channels is crossing symmetry, and it is
# not something to reconstruct from memory -- it is the whole finite-mass effect, and
# getting it backwards would push F the wrong way.
#
# Integrating the annihilation matrix element over the massless two-body final state
# (the LIPS tensor ∫dΠ p₃^μ p₄^ν = (s/96π) g^{μν} + (1/48π) P^μP^ν, checked against its
# own trace) collapses it to a function of s alone:
#
#   W(s) = (4 G_F² s / 3π) · [ (g_L² + g_R²)(s - m_e²) + 6 g_L g_R m_e² ]
#
# which is all the ratio needs, since every constant cancels.

"Integrated annihilation matrix element, up to constants that cancel in the ratio."
_W_ann(s, m2, gL, gR) = s * ((gL^2 + gR^2) * (s - m2) + 6 * gL * gR * m2)

"""
Energy-transfer rate for e⁺e⁻ → νν̄ at temperature T with electron mass m, in
Maxwell-Boltzmann statistics (the Fermi-Dirac correction is the separate constant
`_F_ANN_FD`). Everything constant is dropped: only the ratio is wanted.
"""
function _Q_ann(T, m, gL, gR; n=64)
    m2 = m^2
    pmax = 30T + 5m
    p, w = gauss(n, 0.0, pmax)
    c, wc = gauss(n, -1.0, 1.0)
    acc = 0.0
    @inbounds for i in eachindex(p), j in eachindex(p)
        E1 = sqrt(p[i]^2 + m2)
        E2 = sqrt(p[j]^2 + m2)
        E1 + E2 > 400T && continue
        pref = p[i]^2 * p[j]^2 / (E1 * E2) * exp(-(E1 + E2) / T) * (E1 + E2) * w[i] * w[j]
        for k in eachindex(c)
            s = 2m2 + 2 * (E1 * E2 - p[i] * p[j] * c[k])
            acc += pref * wc[k] * _W_ann(s, m2, gL, gR)
        end
    end
    acc
end

"""
    F_ann(T; flavour = :e)

Finite-electron-mass suppression of the annihilation energy-transfer rate,
F = rate(m_e)/rate(0). Tends to 1 for T ≫ m_e and to 0 once the electrons are gone.
"""
function F_ann(T; flavour=:e, n=64)
    gL, gR = flavour === :e ? (_gL_e, _gR_e) : (_gL_mu, _gR_mu)
    q0 = _Q_ann(T, 0.0, gL, gR; n)
    q0 <= 0 && return 0.0
    _Q_ann(T, _ME_MEV, gL, gR; n) / q0
end

# --- elastic scattering ν e± → ν e± ------------------------------------------
#
# Scattering is harder than annihilation, and for a reason worth stating: in annihilation
# the pair *disappears* and deposits its whole energy E₁+E₂ into the neutrinos, so the
# weight is just the total energy. In scattering nothing is created or destroyed -- the
# neutrino leaves with a different energy than it arrived with, and what matters is the
# *difference*. So the weight is ΔE = E_k2 - E_k1, which is frame-dependent and only
# survives the phase-space integral because the two baths are at different temperatures.
# Set T_e = T_ν and it cancels identically, which is exactly the detailed-balance check
# used below.
#
# The matrix element (exact tree-level trace algebra) is
#
#   |M|² ∝ g_L²(k₁·p₁)(k₂·p₂) + g_R²(k₁·p₂)(k₂·p₁) - g_L g_R m_e² (k₁·k₂)
#
# Rather than reduce this analytically, the collision integral is done directly: sample the
# incoming electron and neutrino, boost to the centre of momentum, scatter, boost back, and
# weight by the energy actually transferred. Deterministic quadrature over the incoming
# energies, Gauss-Legendre over the CM angle.

"Minkowski dot product, (+,-,-,-)."
@inline _mdot(a, b) = a[1] * b[1] - a[2] * b[2] - a[3] * b[3] - a[4] * b[4]

"""
Net energy-transfer rate for ν e± → ν e± at electron temperature `Te` and neutrino
temperature `Tn`, with electron mass `m`. Maxwell-Boltzmann statistics; all constants that
cancel in the ratio F = rate(m_e)/rate(0) are dropped.
"""
function _Q_scat(Te, Tn, m, gL, gR; n=40, nang=24)
    m2 = m^2
    pe, we = gauss(n, 0.0, 25Te + 5m)      # electron momentum
    kn, wn = gauss(n, 1e-6, 25Tn)          # neutrino momentum
    ct, wct = gauss(nang, -1.0, 1.0)       # angle between them
    cs, wcs = gauss(nang, -1.0, 1.0)       # CM scattering angle
    acc = 0.0

    @inbounds for i in eachindex(pe), j in eachindex(kn)
        E1 = sqrt(pe[i]^2 + m2)             # electron
        w1 = pe[i]^2 * we[i]
        w2 = kn[j]^2 * wn[j]
        for a in eachindex(ct)
            # incoming four-momenta: electron along z, neutrino at angle ct
            p1 = (E1, 0.0, 0.0, pe[i])
            st = sqrt(max(1 - ct[a]^2, 0.0))
            k1 = (kn[j], kn[j] * st, 0.0, kn[j] * ct[a])
            s = m2 + 2 * _mdot(p1, k1)
            s <= m2 * 1.0000001 && continue

            # CM boost
            P = (p1[1] + k1[1], p1[2] + k1[2], p1[3] + k1[3], p1[4] + k1[4])
            Ecm = sqrt(s)
            β = (P[2] / P[1], P[3] / P[1], P[4] / P[1])
            β2 = β[1]^2 + β[2]^2 + β[3]^2
            β2 >= 1 && continue
            γ = 1 / sqrt(1 - β2)
            # neutrino CM momentum (massless nu, massive e)
            kcm = (s - m2) / (2Ecm)

            # incoming neutrino direction in CM (boost k1)
            bk = β[1] * k1[2] + β[2] * k1[3] + β[3] * k1[4]
            k1cm0 = γ * (k1[1] - bk)
            f = (γ - 1) / β2 * bk - γ * k1[1]
            k1cm = (k1[2] + f * β[1], k1[3] + f * β[2], k1[4] + f * β[3])
            nk = sqrt(k1cm[1]^2 + k1cm[2]^2 + k1cm[3]^2)
            nk <= 0 && continue
            u = (k1cm[1] / nk, k1cm[2] / nk, k1cm[3] / nk)
            # an orthonormal pair to u
            tmp = abs(u[1]) < 0.9 ? (1.0, 0.0, 0.0) : (0.0, 1.0, 0.0)
            v1 = (u[2] * tmp[3] - u[3] * tmp[2], u[3] * tmp[1] - u[1] * tmp[3], u[1] * tmp[2] - u[2] * tmp[1])
            nv = sqrt(v1[1]^2 + v1[2]^2 + v1[3]^2)
            v = (v1[1] / nv, v1[2] / nv, v1[3] / nv)

            for b in eachindex(cs)
                sc = sqrt(max(1 - cs[b]^2, 0.0))
                d = (cs[b] * u[1] + sc * v[1], cs[b] * u[2] + sc * v[2], cs[b] * u[3] + sc * v[3])
                k2cm = (kcm, kcm * d[1], kcm * d[2], kcm * d[3])
                # boost back to lab
                bk2 = β[1] * k2cm[2] + β[2] * k2cm[3] + β[3] * k2cm[4]
                E_k2 = γ * (k2cm[1] + bk2)
                g2 = (γ - 1) / β2 * bk2 + γ * k2cm[1]
                k2 = (E_k2, k2cm[2] + g2 * β[1], k2cm[3] + g2 * β[2], k2cm[4] + g2 * β[3])
                p2 = (p1[1] + k1[1] - k2[1], p1[2] + k1[2] - k2[2],
                    p1[3] + k1[3] - k2[3], p1[4] + k1[4] - k2[4])

                M2 = gL^2 * _mdot(k1, p1) * _mdot(k2, p2) +
                     gR^2 * _mdot(k1, p2) * _mdot(k2, p1) -
                     gL * gR * m2 * _mdot(k1, k2)

                ΔE = k2[1] - k1[1]                  # energy handed to the neutrino
                # forward minus backward, Maxwell-Boltzmann. Identical temperatures give
                # exactly zero, which is the detailed-balance test.
                fwd = exp(-E1 / Te - k1[1] / Tn)
                bwd = exp(-p2[1] / Te - k2[1] / Tn)
                acc += w1 * w2 * wct[a] * wcs[b] * M2 * ΔE * (fwd - bwd) / (E1 * k1[1])
            end
        end
    end
    acc
end

"""
    F_scat(T; flavour = :e)

Finite-electron-mass suppression of the *scattering* energy-transfer rate,
F = rate(m_e)/rate(0), evaluated in the linear-response regime (a small temperature
difference between the two baths, which is where the T⁴T⁴(T₁-T₂) form is valid anyway).
"""
function F_scat(T; flavour=:e, n=40, nang=24, ε=0.02)
    gL, gR = flavour === :e ? (_gL_e, _gR_e) : (_gL_mu, _gR_mu)
    Te, Tn = T * (1 + ε), T           # small offset drives the transfer
    q0 = _Q_scat(Te, Tn, 0.0, gL, gR; n, nang)
    abs(q0) < 1e-300 && return 0.0
    _Q_scat(Te, Tn, _ME_MEV, gL, gR; n, nang) / q0
end

# --- the decoupling solve ----------------------------------------------------

"""
    neff(; order = 3, T_start = 10.0, T_end = 5e-3)

Solve for N_eff from scratch: evolve T_γ, T_νe and T_νμ through e± annihilation with the
QED-corrected plasma and the weak energy-transfer rates, and read off

    N_eff = (ρ_νe + 2 ρ_νμ) / ρ_γ / [ (7/8) (4/11)^{4/3} ]

Nothing about the answer is put in. If the weak rates were switched off entirely, the
neutrinos would decouple instantaneously and this would return exactly 3.
"""
function neff(; order=3, T_start=10.0, T_end=5e-3, nT=24, reltol=1e-10, abstol=1e-14,
    n_ann=48, n_scat=28, nang=16)
    # F_ann and F_scat are expensive collision integrals and depend only on T. Tabulate.
    Ts = exp.(range(log(T_end / 2), log(T_start * 1.5); length=nT))
    Fa_e = zeros(nT); Fs_e = zeros(nT); Fa_m = zeros(nT); Fs_m = zeros(nT)
    Threads.@threads for i in eachindex(Ts)
        Fa_e[i] = F_ann(Ts[i]; flavour=:e, n=n_ann)
        Fa_m[i] = F_ann(Ts[i]; flavour=:mu, n=n_ann)
        Fs_e[i] = F_scat(Ts[i]; flavour=:e, n=n_scat, nang=nang)
        Fs_m[i] = F_scat(Ts[i]; flavour=:mu, n=n_scat, nang=nang)
    end
    lnT = log.(Ts)
    ia_e = linear_interpolation(lnT, Fa_e; extrapolation_bc=Line())
    ia_m = linear_interpolation(lnT, Fa_m; extrapolation_bc=Line())
    is_e = linear_interpolation(lnT, Fs_e; extrapolation_bc=Line())
    is_m = linear_interpolation(lnT, Fs_m; extrapolation_bc=Line())
    clamp01(v) = clamp(v, 0.0, 1.0)

    # Tabulate the plasma thermodynamics. ρ_plasma calls qed_energy_density, which runs
    # two quadgk integrations, and dρ/dT calls it twice more -- roughly twenty adaptive
    # quadratures per right-hand-side evaluation, thousands of times over. Computing
    # background integrals inside an ODE RHS is the single most reliable way to turn a
    # two-second problem into one that never finishes, and it is the third time in this
    # codebase. Tabulate once; the plasma depends on T_γ and nothing else.
    lnTp = collect(range(log(T_end / 3), log(T_start * 2); length=400))
    ρp_t = [_ρ_plasma(exp(x); order) for x in lnTp]
    Pp_t = [_P_plasma(exp(x); order) for x in lnTp]
    ρp_i = linear_interpolation(lnTp, log.(ρp_t); extrapolation_bc=Line())
    Pp_i = linear_interpolation(lnTp, log.(Pp_t); extrapolation_bc=Line())
    ρp(T) = exp(ρp_i(log(T)))
    Pp(T) = exp(Pp_i(log(T)))
    # dρ/dT from the tabulated log-log slope: dρ/dT = (ρ/T)·dlnρ/dlnT
    function dρp_dT(T)
        h = 1e-4
        x = log(T)
        dlnρ = (ρp_i(x + h) - ρp_i(x - h)) / (2h)
        ρp(T) / T * dlnρ
    end

    pref = _GF^2 / π^5

    function rhs!(du, u, _, lna)
        Tg, Te, Tm = u
        (Tg <= 0 || Te <= 0 || Tm <= 0) && (fill!(du, 0.0); return nothing)
        lg = log(Tg)

        fe = _f_transfer(Tg, Te, clamp01(ia_e(lg)), clamp01(is_e(lg)))
        fm = _f_transfer(Tg, Tm, clamp01(ia_m(lg)), clamp01(is_m(lg)))
        fnn = _f_nunu(Tm, Te)

        δρ_e = pref * (4 * (_gL_e^2 + _gR_e^2) * fe + 2 * fnn)
        δρ_m = pref * (4 * (_gL_mu^2 + _gR_mu^2) * fm - fnn)

        ρ_e = _ρ_nu(Te)
        ρ_m = _ρ_nu(Tm)
        ρ_tot = ρp(Tg) + ρ_e + 2ρ_m
        H = sqrt(8π * _GN * ρ_tot / 3)

        du[1] = (-3 * (ρp(Tg) + Pp(Tg)) - (δρ_e + 2δρ_m) / H) / dρp_dT(Tg)
        du[2] = -Te + δρ_e * Te / (4 * ρ_e * H)
        du[3] = -Tm + δρ_m * Tm / (4 * ρ_m * H)
        nothing
    end

    # Deep in the coupled regime everything shares one temperature.
    u0 = [T_start, T_start, T_start]
    # Integrate in ln a. T_γ falls roughly as 1/a, so ~ln(T_start/T_end) of e-folds, plus
    # slack for the photon reheating during annihilation.
    span = (0.0, log(T_start / T_end) + 1.0)
    cb_T = ContinuousCallback((u, t, i) -> u[1] - T_end, terminate!)
    sol = solve(ODEProblem(rhs!, u0, span), Tsit5(); reltol, abstol, callback=cb_T)

    Tg, Te, Tm = sol.u[end]
    ργ = (π^2 / 15) * Tg^4
    ρν = _ρ_nu(Te) + 2 * _ρ_nu(Tm)
    N = ρν / ργ / ((7 / 8) * (4 / 11)^(4 / 3))

    # The whole trajectory, not just the endpoint. This is what BBN actually needs: the
    # *history* of T_ν/T_γ through the annihilation, rather than a ramp fitted to reproduce
    # its endpoints.
    # Carry ln a as well. It is needed because a·T_ν is *not* exactly constant while the
    # neutrinos are still being heated -- which is precisely the regime this solve exists to
    # describe. Anything that infers the scale factor from T_ν alone (as n_b ∝ T_ν³ does)
    # silently assumes the instantaneous limit it is meant to be correcting.
    traj = [(sol.t[i], sol.u[i][1], sol.u[i][2], sol.u[i][3]) for i in eachindex(sol.u)]
    (; N_eff=N, T_gamma=Tg, T_nue=Te, T_numu=Tm,
        Tnu_over_Tgamma=(Te + 2Tm) / 3 / Tg, trajectory=traj)
end

# --- the decoupling history, for BBN ----------------------------------------
#
# BBN needs T_ν(T_γ) and ρ_ν(T_γ) through e± annihilation. Previously it got them from
# instantaneous decoupling plus a hand-built ramp of N_eff from 3.000 to 3.044 -- which
# reproduced the endpoints by construction and the *path between them* by assumption.
#
# The decoupling solve produces that path directly, so use it. The ramp is gone.

const _DECOUP_CACHE = joinpath(_BBN_DIR, "neutrino_decoupling.dat")

"Compute (or load) the tabulated neutrino decoupling history T_νe(T_γ), T_νμ(T_γ)."
function decoupling_history()
    if !isfile(_DECOUP_CACHE)
        @info "Solving neutrino decoupling (once; cached afterwards)"
        r = neff()
        open(_DECOUP_CACHE, "w") do io
            println(io, "# Neutrino decoupling history, from Cosmic.neff().")
            println(io, "# N_eff = ", r.N_eff)
            println(io, "# ln a   T_gamma [MeV]   T_nue [MeV]   T_numu [MeV]")
            for (l, a, b, c) in r.trajectory
                println(io, "$l  $a  $b  $c")
            end
        end
    end
    d = readdlm(_DECOUP_CACHE, comments=true, comment_char='#')
    p = sortperm(Float64.(d[:, 2]))              # ascending in T_γ
    lnTg = log.(Float64.(d[p, 2]))
    lna = Float64.(d[p, 1])
    Te = Float64.(d[p, 3])
    Tm = Float64.(d[p, 4])

    # The terminating callback writes the final state twice, leaving two knots at the
    # same T_γ. Interpolations nudges one apart ("deduplicated" warning) and the last
    # interval acquires a slope of zero -- which Line() extrapolation then continues
    # forever, freezing T_ν *itself* (not the ratio) below the table and inflating ρ_ν
    # by an order of magnitude where the last of the deuterium is still burning. Keep
    # strictly increasing knots only.
    keep = [true; diff(lnTg) .> 1e-10]
    lnTg = lnTg[keep]; lna = lna[keep]; Te = Te[keep]; Tm = Tm[keep]

    # The solve's ln a carries an arbitrary offset. Anchor it with the *photon* temperature,
    # not the neutrino one: once e± annihilation is over, a·T_γ is exactly constant, and
    # today a = 1 with T_γ = T_cmb.
    #
    # Anchoring on a·T_ν = T_ν0 with T_ν0 = (4/11)^{1/3}·T_cmb is the obvious thing to write
    # and it is circular -- (4/11)^{1/3} *is* the instantaneous-decoupling answer, and this
    # whole calculation exists to say the neutrinos are hotter than that (0.7169, not
    # 0.7138). Using it here puts a 0.4% error into the scale factor, 1.3% into n_b, and
    # ~2% into D/H, while leaving every temperature ratio looking perfectly correct.
    T_γ0 = 2.7255 * 8.617333262e-11              # MeV
    i_end = argmin(lnTg)                          # coldest point = latest time
    δ = log(T_γ0) - lnTg[i_end] - lna[i_end]
    lna_abs = lna .+ δ

    # Interpolate the temperatures in log-log. Outside the table both ends are exact
    # power laws (T_ν ∝ T_γ with a frozen ratio), so a Line() extrapolation of ln T_ν
    # against ln T_γ continues them correctly. Linear-in-T_ν extrapolation does not:
    # below the table floor it crosses zero and hands ρ_ν a negative temperature.
    (linear_interpolation(lnTg, log.(Te); extrapolation_bc=Line()),
        linear_interpolation(lnTg, log.(Tm); extrapolation_bc=Line()),
        linear_interpolation(lnTg, lna_abs; extrapolation_bc=Line()))
end

const _DEC_E = Ref{Any}(nothing)
const _DEC_M = Ref{Any}(nothing)
const _DEC_A = Ref{Any}(nothing)

"""
    nu_temperatures(T_γ)

(T_νe, T_νμ) at photon temperature `T_γ` (MeV), from the solved decoupling history.
This *is* the non-instantaneous decoupling: no ramp, no (4/11)^{1/3} assumption.
"""
function nu_temperatures(T_γ)
    _load_decoupling()
    x = log(T_γ)
    (exp(_DEC_E[](x)), exp(_DEC_M[](x)))
end

function _load_decoupling()
    if _DEC_E[] === nothing
        a, b, c = decoupling_history()
        _DEC_E[] = a
        _DEC_M[] = b
        _DEC_A[] = c
    end
    nothing
end

"""
    scale_factor_at(T_γ)

Scale factor (a = 1 today) at photon temperature `T_γ` (MeV), from the decoupling solve.
Normalised by a·T_γ = const *after* e± annihilation, which is the only epoch in which
that statement is exactly true.
"""
function scale_factor_at(T_γ)
    _load_decoupling()
    exp(_DEC_A[](log(T_γ)))
end


# --- cached result ------------------------------------------------------------
#
# The solve takes minutes (the collision integrals dominate) and depends on nothing a user
# can vary within the Standard Model -- it is a fixed prediction of the SM plus gravity. So
# the number is cached, and `neff()` regenerates it on demand.

const _NEFF_CACHE = joinpath(_BBN_DIR, "neff_sm.dat")

"""
    N_eff_sm()

The Standard-Model prediction for N_eff, from [`neff`](@ref). Cached: the underlying solve
is a several-minute collision-integral computation whose answer never changes.

    cosmology(N_eff = :derived)

uses this instead of the literature value.
"""
function N_eff_sm()
    if isfile(_NEFF_CACHE)
        return parse(Float64, strip(split(read(_NEFF_CACHE, String), "\n")[end-1]))
    end
    r = neff()
    open(_NEFF_CACHE, "w") do io
        println(io, "# Standard-Model N_eff, computed by Cosmic.neff().")
        println(io, "# QED plasma pressure (closed form) + neutrino decoupling with")
        println(io, "# FeynCalc-derived collision integrals. Delete to recompute.")
        println(io, r.N_eff)
    end
    r.N_eff
end
