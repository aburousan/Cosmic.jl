"""
Big Bang nucleosynthesis: the primordial helium fraction, solved rather than assumed.

Y_He enters recombination through f_He = n_He/n_H, and every Boltzmann code needs it.
CAMB and CLASS both get it by interpolating a precomputed table (CLASS ships
`external/bbn/sBBN.dat`). That is a lookup, not a calculation, and it locks you to
whatever nuclear physics the table's author assumed. Here the nuclear network is
integrated directly, so Y_p and D/H come out of the same code that produced them.

Three pieces:

  * **The weak rates n <-> p.** These are *derived*: the Born rate is a single
    integral over the electron momentum (Pitrou, Coc, Uzan & Vangioni, Phys. Rep.
    754 (2018), eqs 77-78 and 81),

        Λ(x, x_ν) = ∫₀^∞ dp p² [ χ(E) + χ(-E) ],   E = √(p²+1)
        χ(E)      = E_ν² · g_ν(E_ν, x_ν) · g(-E, x),   E_ν = E - s·q

    with s = +1 for n→p and -1 for p→n, q = (m_n - m_p)/m_e, and g, g_ν the
    Fermi-Dirac occupations at the photon and neutrino temperatures. The two terms
    χ(E) + χ(-E) are not a trick of notation: evaluating χ at negative E swaps the
    Pauli-blocking factor [1 - f_e] for the positron occupation f_e, so the single
    integral covers all three channels at once -- free decay n → p e ν̄, electron
    capture p e → n ν, and positron capture n e⁺ → p ν̄ -- each with the correct
    blocking. Writing them out separately is a well-known way to get a sign wrong.

    The overall constant is fixed by the measured neutron lifetime: as T → 0 the
    integral reduces to the vacuum phase-space factor ∫₁^q dE E√(E²-1)(q-E)², so
    K = 1/(τ_n · that). Every T-dependence is then a prediction. G_F, V_ud and g_A
    never appear -- they are exactly what τ_n measures.

  * **The nuclear network.** Sixty-three reactions among twelve nuclides (n, p, d, t,
    ³He, ⁴He, ⁶He, ⁶Li, ⁷Li, ⁸Li, ⁷Be, ⁸B) -- the twelve-reaction core that fixes Y_p
    and D/H (`network = :small` runs just that) plus the full PRIMAT-compilation set
    that fixes the lithium isotopes. Nuclear cross-sections are *measured*, not
    derivable -- no code on earth computes ⁴He photodisintegration from QCD -- so the
    forward rates are tabulated compilations: PRIMAT throughout by default (one
    consistent chain; `rates = :parthenope` swaps the 12 core reactions for
    PArthENoPE 3.0's post-LUNA fits). They are verified as far as theory allows: the
    three deuterium-burning rates are recomputed in the test suite by direct Gamow
    integration of the published LUNA-era S-factors (0.1% agreement, checked
    independently with the Wolfram engine), and every reverse rate follows from
    detailed balance,

        λ_rev / λ_fwd = α · T9^β · exp(γ/T9),    γ = -Q/(k·10⁹K)

    which *is* derivable: `detailed_balance` computes (α, β, γ) from AME masses and
    ground-state spins, and all 63 stored coefficient sets are checked against it.
    β = (3/2)·Δn: 3/2 for radiative capture, 0 for 2→2 rearrangements, -3/2 and
    beyond for the breakup channels.

    The tabulated 1σ error envelopes propagate: `bbn_mc` redraws every rate
    lognormally within its envelope (and τ_n within its uncertainty) and returns
    distributions for Y_p, D/H, ³He/H and ⁷Li/H.

  * **The background.** ρ = ρ_γ + ρ_e± + ρ_ν, with the electron-positron term the
    full Fermi-Dirac integral, not a g_* step function -- e± annihilation happens
    *during* nucleosynthesis and its width is what sets T_ν/T_γ. That ratio comes
    from entropy conservation in the coupled γ-e± plasma, which reproduces
    (4/11)^{1/3} in the limit without it being put in by hand.

Validated against CLASS's `sBBN_2025_primat.dat`.
"""

using OrdinaryDiffEq: ODEProblem, solve, TRBDF2, ReturnCode
using ADTypes: AutoFiniteDiff
using Interpolations: linear_interpolation, Line
using DelimitedFiles: readdlm
using QuadGK: quadgk, gauss
using SpecialFunctions: gamma
import Random
import Statistics

const _BBN_DIR = joinpath(@__DIR__, "..", "data", "bbn")

# --- Nuclide bookkeeping -----------------------------------------------------
#
# Order matters only in that the ODE state vector follows it.

const NUCLIDES = (:n, :p, :d, :t, :He3, :He4, :Li6, :Li7, :Li8, :Be7, :B8, :He6)
const NUC_A = (1, 1, 2, 3, 3, 4, 6, 7, 8, 7, 8, 6)    # mass number
const NUC_Z = (0, 1, 1, 1, 2, 2, 3, 3, 3, 4, 5, 2)    # charge
const IDX = NamedTuple{NUCLIDES}(ntuple(i -> i, length(NUCLIDES)))

# Atomic masses in u, for the mass fraction at the end. m_H is the *atom*, which is
# the convention CLASS's sBBN table uses; getting this wrong shifts Y_p by 0.6% --
# enough to matter and easy to miss, since both numbers "look like" 0.246.
const M_H_ATOM = 1.00782503
const M_HE4_ATOM = 4.00260325

# --- Reaction network --------------------------------------------------------

"""
One thermonuclear reaction: reactants → products, with a tabulated forward rate and
detailed-balance coefficients for the reverse.

`radiative` marks a+b → c+γ, whose reverse is a one-body photodisintegration and so
carries different units (1/s rather than cm³/mol/s) and the T9^{3/2} phase-space factor.
"""
struct Reaction{I}
    name::String
    reactants::Vector{Int}     # indices into NUCLIDES, repeated for multiplicity
    products::Vector{Int}
    rate::I                    # log-log interpolant: T9 -> N_A<σv>
    α::Float64
    β::Float64
    γ::Float64
    radiative::Bool
end

const _RATE_TABLES = Dict{String,NTuple{3,Vector{Float64}}}()

function _rate_table(name, subdir="")
    get!(_RATE_TABLES, joinpath(subdir, name)) do
        raw = readdlm(joinpath(_BBN_DIR, "rates", subdir, name * ".txt"), comments=true, comment_char='#')
        T9 = Float64.(raw[:, 1])
        R = max.(Float64.(raw[:, 2]), 1e-300)
        # Column 3 is the 1σ multiplicative envelope: rate·σ^±1 brackets the band.
        σ = size(raw, 2) >= 3 ? max.(Float64.(raw[:, 3]), 1.0) : ones(length(T9))
        (log.(T9), log.(R), log.(σ))
    end
end

"""
Forward-rate interpolant for `name`, displaced by `ξ` standard deviations of its
tabulated error envelope: λ(T9) = median(T9)·σ(T9)^ξ. ξ = 0 is the median rate;
drawing ξ ~ N(0,1) per reaction is the PArthENoPE/PRIMAT Monte-Carlo convention.
"""
function _load_rate(name; ξ=0.0, subdir="")
    lnT9, lnR, lnσ = _rate_table(name, subdir)
    # Interpolate in log-log: the rates span ~20 decades and are smooth powers of T9,
    # so linear interpolation of the logs is far more faithful than of the values.
    vals = ξ == 0.0 ? lnR : lnR .+ ξ .* lnσ
    itp = linear_interpolation(lnT9, vals; extrapolation_bc=Line())
    lo, hi = exp(lnT9[1]), exp(lnT9[end])
    T -> exp(itp(log(clamp(T, lo, hi))))
end

function _network(which::Symbol=:full; rates::Symbol=:primat, ξ=nothing)
    db = readdlm(joinpath(_BBN_DIR, "detailed_balance.csv"), ','; header=true)[1]
    coef = Dict(string(db[i, 1]) => (Float64(db[i, 3]), Float64(db[i, 4]), Float64(db[i, 5]))
                for i in axes(db, 1))

    spec = [
        # name             reactants           products          radiative
        ("n_p__d_g", [:n, :p], [:d], true),
        ("d_p__He3_g", [:d, :p], [:He3], true),
        ("d_d__He3_n", [:d, :d], [:He3, :n], false),
        ("d_d__t_p", [:d, :d], [:t, :p], false),
        ("t_p__a_g", [:t, :p], [:He4], true),
        ("t_d__a_n", [:t, :d], [:He4, :n], false),
        ("t_a__Li7_g", [:t, :He4], [:Li7], true),
        ("He3_n__t_p", [:He3, :n], [:t, :p], false),
        ("He3_d__a_p", [:He3, :d], [:He4, :p], false),
        ("He3_a__Be7_g", [:He3, :He4], [:Be7], true),
        ("Be7_n__Li7_p", [:Be7, :n], [:Li7, :p], false),
        ("Li7_p__a_a", [:Li7, :p], [:He4, :He4], false),
    ]

    # The 51 additional reactions of the full network (PRIMAT compilation, via
    # PRyMordial): lithium-6/8, beryllium, boron-8 and helium-6 channels. They
    # move Y_p and D/H by well under 0.1% -- their job is Li6, Li7 and the tail
    # isotopes, and the honest error bars that come from having the full set.
    spec_full = [
        ("Be7_n__a_a", [:Be7, :n], [:He4, :He4], false),
        ("Be7_d__a_a_p", [:Be7, :d], [:He4, :He4, :p], false),
        ("d_a__Li6_g", [:d, :He4], [:Li6], true),
        ("Li6_p__Be7_g", [:Li6, :p], [:Be7], true),
        ("Li6_p__He3_a", [:Li6, :p], [:He3, :He4], false),
        ("B8_n__a_a_p", [:B8, :n], [:He4, :He4, :p], false),
        ("Li6_He3__a_a_p", [:Li6, :He3], [:He4, :He4, :p], false),
        ("Li6_t__a_a_n", [:Li6, :t], [:He4, :He4, :n], false),
        ("Li6_t__Li8_p", [:Li6, :t], [:Li8, :p], false),
        ("Li7_He3__Li6_a", [:Li7, :He3], [:Li6, :He4], false),
        ("Li8_He3__Li7_a", [:Li8, :He3], [:Li7, :He4], false),
        ("Be7_t__Li6_a", [:Be7, :t], [:Li6, :He4], false),
        ("B8_t__Be7_a", [:B8, :t], [:Be7, :He4], false),
        ("B8_n__Li6_He3", [:B8, :n], [:Li6, :He3], false),
        ("B8_n__Be7_d", [:B8, :n], [:Be7, :d], false),
        ("Li6_t__Li7_d", [:Li6, :t], [:Li7, :d], false),
        ("Li6_He3__Be7_d", [:Li6, :He3], [:Be7, :d], false),
        ("Li7_He3__a_a_d", [:Li7, :He3], [:He4, :He4, :d], false),
        ("Li8_He3__a_a_t", [:Li8, :He3], [:He4, :He4, :t], false),
        ("Be7_t__a_a_d", [:Be7, :t], [:He4, :He4, :d], false),
        ("Be7_t__Li7_He3", [:Be7, :t], [:Li7, :He3], false),
        ("B8_d__Be7_He3", [:B8, :d], [:Be7, :He3], false),
        ("B8_t__a_a_He3", [:B8, :t], [:He4, :He4, :He3], false),
        ("Be7_He3__p_p_a_a", [:Be7, :He3], [:p, :p, :He4, :He4], false),
        ("d_d__a_g", [:d, :d], [:He4], true),
        ("He3_He3__a_p_p", [:He3, :He3], [:He4, :p, :p], false),
        ("Be7_p__B8_g", [:Be7, :p], [:B8], true),
        ("Li7_d__a_a_n", [:Li7, :d], [:He4, :He4, :n], false),
        ("d_n__t_g", [:d, :n], [:t], true),
        ("t_t__a_n_n", [:t, :t], [:He4, :n, :n], false),
        ("He3_n__a_g", [:He3, :n], [:He4], true),
        ("He3_t__a_d", [:He3, :t], [:He4, :d], false),
        ("He3_t__a_n_p", [:He3, :t], [:He4, :n, :p], false),
        # Two independent determinations of Li7(t,2n)2α ship in the compilation
        # (as for Be7+He3 below); both enter, exactly as in PRyMordial's network.
        ("Li7_t__a_a_n_n_alt", [:Li7, :t], [:He4, :He4, :n, :n], false),
        ("Li7_He3__a_a_n_p", [:Li7, :He3], [:He4, :He4, :n, :p], false),
        ("Li8_d__Li7_t", [:Li8, :d], [:Li7, :t], false),
        ("Be7_t__a_a_n_p", [:Be7, :t], [:He4, :He4, :n, :p], false),
        ("Be7_He3__a_a_p_p", [:Be7, :He3], [:He4, :He4, :p, :p], false),
        ("Li6_n__t_a", [:Li6, :n], [:t, :He4], false),
        ("He3_t__Li6_g", [:He3, :t], [:Li6], true),
        ("a_n_p__Li6_g", [:He4, :n, :p], [:Li6], true),
        ("Li6_n__Li7_g", [:Li6, :n], [:Li7], true),
        ("Li6_d__Li7_p", [:Li6, :d], [:Li7, :p], false),
        ("Li6_d__Be7_n", [:Li6, :d], [:Be7, :n], false),
        ("Li7_n__Li8_g", [:Li7, :n], [:Li8], true),
        ("Li7_d__Li8_p", [:Li7, :d], [:Li8, :p], false),
        ("Li8_p__a_a_n", [:Li8, :p], [:He4, :He4, :n], false),
        ("a_n_n__He6_g", [:He4, :n, :n], [:He6], true),
        ("p_p_n__d_p", [:p, :p, :n], [:d, :p], false),
        ("Li7_t__a_a_n_n", [:Li7, :t], [:He4, :He4, :n, :n], false),
        ("Li7_p__a_a_g", [:Li7, :p], [:He4, :He4], true),
    ]

    # One compilation for the whole chain. The 51 auxiliary reactions exist only in
    # the PRIMAT compilation, so `rates = :primat` (default) takes the 12 core
    # reactions from PRIMAT as well -- a consistent chain, and the configuration in
    # which this network reproduces PRyMordial to 0.004% (Y_p), 0.03% (D/H) and
    # 1.1% (Li7). `rates = :parthenope` swaps the core 12 for PArthENoPE 3.0's
    # post-LUNA fits (D/H moves up ~2.6%, Li7 down ~13% -- that is the genuine
    # spread between the two compilations, not an error in either).
    rates in (:primat, :parthenope) ||
        throw(ArgumentError("rates must be :primat or :parthenope"))
    ncore = length(spec)
    chosen = which === :small ? spec : vcat(spec, spec_full)
    [begin
         α, β, γ = coef[nm]
         sub = (rates === :primat && j <= ncore) ? "primat_core" : ""
         Reaction(nm, [IDX[s] for s in rs], [IDX[s] for s in ps],
             _load_rate(nm; ξ=ξ === nothing ? 0.0 : ξ[j], subdir=sub), α, β, γ, rad)
     end for (j, (nm, rs, ps, rad)) in enumerate(chosen)]
end

# --- Detailed balance, derived -------------------------------------------------
#
# The reverse-rate coefficients stored in detailed_balance.csv are not free data:
# they are fixed by statistical mechanics. For Σr → Σp (+γ) in kinetic and chemical
# equilibrium,
#
#   λ_rev/λ_fwd = (sym_p/sym_r) · (Πg_r/Πg_p) · (Πm_r/Πm_p)^{3/2} · C₀^Δ ·
#                 T9^{3Δ/2} · e^{-Q/kT9},     Δ = n_reactants - n_products
#
# with C₀ = [(m_u·k_B·10⁹K)/(2πħ²)]^{3/2}/N_A the one-particle phase-space density
# per mole at T9 = 1. Nothing else survives: the photon carries no chemical
# potential, and every power of density cancels against the N_A^{n-1} convention of
# the rate tables. `detailed_balance` derives (α, β, γ) from nuclear masses and
# ground-state spins alone, and the test suite checks all 63 stored coefficient
# sets against it -- which is as much verification as the *reverse* half of a rate
# file can ever get from theory. (The forward halves are measured cross-sections;
# no amount of QFT reproduces an R-matrix fit.)

# Nuclear masses in u: AME2020 atomic masses minus Z electron masses (electron
# binding is ~1e-8 u, far below the 5-digit precision of the stored coefficients).
const NUC_M_U = (1.00866491595, 1.00727646688, 2.01355321, 3.01550071, 3.01493216,
    4.00150609, 6.01347715, 7.01435770, 8.02084050, 7.01473440,
    8.02186442, 6.01778873)
const NUC_G = (2, 2, 3, 2, 2, 1, 3, 4, 5, 4, 5, 1)   # 2J+1, ground states
const _C0_MOL = 9.86787e9        # [(m_u·k_B·10⁹K)/(2πħ²)]^{3/2}/N_A, mol/cm³
const _T9_MEV = 0.086173332621   # k_B·10⁹K in MeV

"""
    detailed_balance(reactants, products) -> (α, β, γ)

Reverse-rate coefficients derived from masses and spins (see the note above).
Takes nuclide indices as stored in a `Reaction`; the photon, if any, is simply
absent from the product list.
"""
function detailed_balance(rs::Vector{Int}, ps::Vector{Int})
    Δ = length(rs) - length(ps)
    G = prod(i -> NUC_G[i], rs) / prod(i -> NUC_G[i], ps)
    M = prod(i -> NUC_M_U[i], rs) / prod(i -> NUC_M_U[i], ps)
    sym = _symmetry(ps) / _symmetry(rs)
    Q_MeV = (sum(i -> NUC_M_U[i], rs) - sum(i -> NUC_M_U[i], ps)) * 931.49410242
    (sym * G * M^1.5 * _C0_MOL^Δ, 1.5Δ, -Q_MeV / _T9_MEV)
end

"Number of ways to pick the (possibly identical) reactants: the 1/n! double-counting factor."
_symmetry(idxs) = prod(factorial, values(_counts(idxs)); init=1)
function _counts(idxs)
    c = Dict{Int,Int}()
    for i in idxs
        c[i] = get(c, i, 0) + 1
    end
    c
end

# --- Weak rates --------------------------------------------------------------

const _Q_NP_MEV = 1.29333236          # m_n - m_p
const _ME_MEV = 0.51099895
const _q_np = _Q_NP_MEV / _ME_MEV     # 2.5310, in units of m_e

"Fermi-Dirac occupation 1/(e^{x·E}+1), safe for either sign of E."
@inline _fd(E, x) = (a = x * E; a > 300 ? 0.0 : a < -300 ? 1.0 : 1 / (exp(a) + 1))

# --- Zero-temperature QED corrections to the weak rates ----------------------
#
# The Born rate alone leaves Y_p about 1.7% low, and it is worth being precise about
# why, because the obvious argument says it should not matter at all: if the
# correction were a constant factor C, then K = 1/(τ_n·C·λ₀) and Γ = K·C·Λ would have
# C cancel exactly, and the whole thing would be invisible.
#
# It does not cancel, because C is *not* constant -- and specifically because the
# Coulomb factor is large exactly where the decay integral lives and small exactly
# where the thermal rates live. In free decay the electron is slow (E_e ≤ q ≈ 2.53 m_e,
# so b = v/c reaches at most 0.92 and is often much less) and F(b) → 2πα/b as b → 0.
# In the thermal plasma at T ~ 1 MeV the electrons are relativistic, b → 1, and
# F → 1.02. So the normalisation integral is enhanced *more* than the rates it
# normalises, the net weak rate goes down, freeze-out happens earlier, and n/p -- and
# therefore Y_p -- goes up.

const _ALPHA_FS = 7.2973525693e-3
const _MP_MEV = 938.27208816
const _R_PROTON = 0.8414e-15          # m, CODATA 2018 charge radius
const _LAMBDA_C = 3.8615926796e-13    # m, ħ/(m_e c)

"Dilogarithm Li₂(x) for 0 ≤ x ≤ 1, by series with reflection near 1."
function _li2(x)
    x < 0 && return -_li2(-x / (1 - x)) - 0.5 * log1p(-x)^2
    if x > 0.5
        # Li₂(x) = π²/6 - ln x ln(1-x) - Li₂(1-x)
        y = 1 - x
        y <= 0 && return π^2 / 6
        return π^2 / 6 - log(x) * log(y) - _li2(y)
    end
    s = 0.0
    t = x
    for k in 1:60
        s += t / k^2
        t *= x
        abs(t) < 1e-18 && break
    end
    s
end

"""
Fermi-Coulomb factor F(b) for an electron emitted into the field of the daughter
proton (Sirlin 1967; Pitrou et al. Phys. Rep. §III.D):

    F(b) = (1 + Γ/2)·4·(2 r_p b/λ_C)^{2Γ} / Γ(3+2Γ)² · e^{πα/b} / (1-b²)^Γ · |Γ(1+Γ+iα/b)|²

with Γ = √(1-α²) - 1. It diverges as b → 0, but p² = b²E² vanishes faster, so the
product is integrable. Reduces to the familiar 2πη/(1 - e^{-2πη}) when α/b is small.
"""
function _fermi_coulomb(b)
    b <= 0 && return 1.0
    Γc = sqrt(1 - _ALPHA_FS^2) - 1
    η = _ALPHA_FS / b
    (1 + Γc / 2) * 4 * (2 * _R_PROTON * b / _LAMBDA_C)^(2Γc) / gamma(3 + 2Γc)^2 *
    exp(π * η) / (1 - b^2)^Γc * abs2(gamma(complex(1 + Γc, η)))
end

"""
Resummed zero-temperature radiative correction R(b, y, E) (Sirlin's outer function
plus the Czarnecki et al. 2004 short-distance constants; Pitrou et al. eqs 101-105).
`y` is the neutrino energy and `E` the electron energy, both in units of m_e.

The two constant factors below cancel in the τ_n normalisation -- they are carried
anyway so that the rate is right in absolute terms, not just in ratio.
"""
function _rad_corr(b, y, E)
    b <= 0 && return 1.0
    Rd = atanh(b) / b
    Q = _Q_NP_MEV
    sirlin = 3 * log(_MP_MEV / _ME_MEV) - 3 / 4 +
             4 * (Rd - 1) * (y / (3E) - 3 / 2 + log(2y)) +
             Rd * (2 * (1 + b^2) + y^2 / (6E^2) - 4 * b * Rd) -
             (4 / b) * _li2(2b / (1 + b))
    outer = 1 + _ALPHA_FS / (2π) * (sirlin - 3 * log(_MP_MEV / (2Q)))

    mA = 1.2e3                                   # GeV-scale hadronic matching
    long = 1.02094 + (_ALPHA_FS / π) * 0.891 - 0.00043
    short = 1.02248 + (log(_MP_MEV / mA) - 0.34) / (134 * 2π) - 0.0001
    outer * long * short
end

"""
Correction factor multiplying the Born integrand, for one term of χ(±E).

`sgn_E = +1` for the χ(E) term, `-1` for χ(-E). The Coulomb factor applies only when
the produced charged lepton is the *electron*, which is the case when
`direction·sgn_E > 0`; a positron is repelled, not attracted, and PRIMAT drops the
factor there rather than flipping its sign.
"""
@inline function _weak_correction(b, E, y, direction, sgn_E)
    F = (direction * sgn_E > 0) ? _fermi_coulomb(b) : 1.0
    F * _rad_corr(b, y, E)
end

# --- Finite-nucleon-mass corrections ----------------------------------------
#
# The Born rate treats the nucleon as an infinitely heavy, static source: the neutron
# does not recoil, and its thermal motion is ignored. Both are corrections of order
# T/m_N and E_e/m_N -- a few times 1e-3 -- and they are the single largest remaining
# effect on Y_p, worth about +0.4%.
#
# Working them out means expanding the collision integral to first order in 1/m_N
# (a Fokker-Planck expansion, Pitrou et al. Phys. Rep. §III.G). What comes out is a
# combination of Fermi-Dirac kernels weighted by three couplings built from g_A and the
# weak-magnetism constant δκ = κ_p - κ_n. The nucleon structure enters here and nowhere
# else in the Born rate, which is why g_A finally appears explicitly: unlike the overall
# strength, this piece is *not* absorbed into τ_n.

const _GA = 1.2756                                  # axial coupling, PDG
const _DELTA_KAPPA = (2.79284734463 - 1) - (-1.91304273)   # κ_p - κ_n, weak magnetism
const _MN_MEV = 939.56542052

# Fermi-Dirac kernels with μ_ν = 0, from the Fokker-Planck expansion. Each is written
# in whichever of the two algebraically identical forms does not overflow: for u ≥ 0 the
# e^{2u}, e^{3u} in the naive expression blow past Float64 long before the function
# itself becomes negligible, so numerator and denominator are divided through by e^{3u}.
@inline function _fdk(k::Symbol, E, x)
    u = x * E
    u > 300 && return 0.0
    if u >= 0
        t = exp(-u)
        d2 = (1 + t)^2
        d3 = (1 + t)^3
        k === :e2p0 && return E^2 * t / (1 + t)
        k === :e3p0 && return E^3 * t / (1 + t)
        k === :e2p1 && return E * (2t^2 + (2 - u) * t) / d2
        k === :e3p1 && return E^2 * (3t^2 + (3 - u) * t) / d2
        k === :e4p1 && return E^3 * (4t^2 + (4 - u) * t) / d2
        k === :e2p2 && return ((u * (u - 4) + 2) * t + (4 - u * (u + 4)) * t^2 + 2t^3) / d3
        k === :e3p2 && return E * ((12 - u * (u + 6)) * t^2 + (u * (u - 6) + 6) * t + 6t^3) / d3
        k === :e4p2 && return E^2 * ((24 - u * (u + 8)) * t^2 + (u - 6) * (u - 2) * t + 12t^3) / d3
    else
        e = exp(u)
        d2 = (e + 1)^2
        d3 = (e + 1)^3
        k === :e2p0 && return E^2 / (e + 1)
        k === :e3p0 && return E^3 / (e + 1)
        k === :e2p1 && return E * (2 + e * (2 - u)) / d2
        k === :e3p1 && return E^2 * (3 + e * (3 - u)) / d2
        k === :e4p1 && return E^3 * (4 + e * (4 - u)) / d2
        k === :e2p2 && return ((u * (u - 4) + 2) * e^2 + (4 - u * (u + 4)) * e + 2) / d3
        k === :e3p2 && return E * ((12 - u * (u + 6)) * e + e^2 * (u * (u - 6) + 6) + 6) / d3
        k === :e4p2 && return E^2 * ((24 - u * (u + 8)) * e + e^2 * (u - 6) * (u - 2) + 12) / d3
    end
    0.0
end

"The three nucleon-structure couplings, for n→p (s=+1) or p→n (s=-1)."
@inline function _fm_couplings(s)
    g = _GA
    dk = _DELTA_KAPPA
    den = 1 + 3g^2
    f1 = ((1 + s * g)^2 + 2 * dk * s * g) / den
    f2 = ((1 - s * g)^2 - 2 * dk * s * g) / den
    f3 = (g^2 - 1) / den
    f1, f2, f3
end

"Finite-mass contribution to χ, Pitrou et al. §III.G (as implemented in PRyMordial)."
function _chi_fm(en, pe, x, znu, s)
    f1, f2, f3 = _fm_couplings(s)
    M = (_MP_MEV + _MN_MEV - s * _Q_NP_MEV) / (2 * _ME_MEV)
    Eν = en - s * _q_np
    F2 = _fd(-en, x)

    (f1 * _fdk(:e2p0, Eν, znu) * F2 * (pe^2 / (M * en))
     +
     f2 * _fdk(:e3p0, Eν, znu) * F2 * (-1 / M)
     +
     (f1 + f2 + f3) / (2 * x * M) * (_fdk(:e4p2, Eν, znu) * F2 + _fdk(:e2p2, Eν, znu) * F2 * pe^2)
     +
     (f1 + f2 + f3) / (2 * M) * (_fdk(:e4p1, Eν, znu) * F2 + _fdk(:e2p1, Eν, znu) * F2 * pe^2)
     -
     (f1 + f2) / (x * M) * (_fdk(:e3p1, Eν, znu) * F2 + _fdk(:e2p1, Eν, znu) * F2 * pe^2 / (-en))
     -
     f3 * 3 / (x * M) * _fdk(:e2p0, Eν, znu) * F2
     +
     f3 / (3M) * _fdk(:e3p1, Eν, znu) * F2 * pe^2 / en
     +
     f3 * 2 / (2 * x * 3 * M) * _fdk(:e3p2, Eν, znu) * F2 * pe^2 / en
     -
     (f1 + f2 + f3) * 3 / (2x) * (1 - (_MN_MEV / _MP_MEV)^s) * (_fdk(:e2p1, Eν, znu) * F2))
end

"""
Neutron-decay phase-space factor, with every correction the thermal rates carry.

    F_n = F_n^rad + F_n^FM

This is what fixes K = 1/(τ_n·F_n), and it *must* carry the same corrections as the
thermal integrand. Normalising a corrected rate with an uncorrected phase-space factor
is precisely the mistake that shifts Y_p by ~1.7%, in the direction that looks fine.
"""
const _FN_RAD = quadgk(E -> begin
        p = sqrt(E^2 - 1)
        b = p / E
        E * p * (_q_np - E)^2 * _weak_correction(b, E, _q_np - E, 1, 1)
    end, 1.0, _q_np; rtol=1e-11)[1]

const _FN_FM = quadgk(pe -> begin
        E = sqrt(pe^2 + 1)
        b = pe / E
        mn_me = _MN_MEV / _ME_MEV
        f1, f2, f3 = _fm_couplings(1)
        d = E - _q_np
        χ = (f1 * d^2 * (pe^2 / (mn_me * E)) - f2 / mn_me * d^3 +
             (f1 + f2 + f3) / (2 * mn_me) * (4d^3 + 2d * pe^2) +
             f3 / mn_me * d^2 * pe^2 / E)
        pe^2 * χ * _weak_correction(b, E, abs(E - _q_np), 1, 1)
    end, 0.0, sqrt(_q_np^2 - 1); rtol=1e-11)[1]

const _FN = _FN_RAD + _FN_FM

"Born-only phase-space factor, for comparison: the literature value is 1.6360."
const _LAMBDA0 = quadgk(E -> E * sqrt(E^2 - 1) * (_q_np - E)^2, 1.0, _q_np; rtol=1e-12)[1]

"""
    weak_rate(T_γ, T_ν, τ_n; direction)

n→p (`direction = +1`) or p→n (`direction = -1`) conversion rate in 1/s.

    Γ = (1/(τ_n λ₀)) ∫₀^∞ dp p² [ χ(E) + χ(-E) ],   E = √(p²+1)
    χ(E) = E_ν² · g_ν(E_ν, x_ν) · g(-E, x),   E_ν = E - s·q

Temperatures in MeV. See the module docstring for why the χ(E) + χ(-E) form is the
right way to write this.
"""
function weak_rate(T_γ, T_ν, τ_n; direction=1, thermal=true)
    x = _ME_MEV / T_γ
    xν = _ME_MEV / T_ν
    s = direction

    χ(E) = begin
        Eν = E - s * _q_np
        Eν^2 * _fd(Eν, xν) * _fd(-E, x)
    end

    # Gauss-Legendre on p ∈ [0, pmax]. The integrand is cut off by the Fermi factors
    # at a few × T/m_e, so pmax = 50 + 50·(T/m_e) is generous at every temperature
    # that matters and the rule converges to machine precision.
    pmax = 50 * (1 + max(T_γ, T_ν) / _ME_MEV)
    p, w = gauss(200, 0.0, pmax)
    acc = 0.0
    @inbounds for i in eachindex(p)
        E = sqrt(p[i]^2 + 1)
        b = p[i] / E
        cp = _weak_correction(b, E, abs(s * _q_np - E), s, 1)
        cm = _weak_correction(b, E, abs(s * _q_np + E), s, -1)
        # Born + finite-nucleon-mass, each carrying the Coulomb and radiative factors.
        acc += w[i] * p[i]^2 * (
            (χ(E) + _chi_fm(E, p[i], x, xν, s)) * cp +
            (χ(-E) + _chi_fm(-E, p[i], x, xν, s)) * cm)
    end
    # Finite-temperature (plasma) QED correction. It vanishes at T = 0, so it does not
    # enter F_n -- it changes the rates without changing the decay they are normalised
    # against, and therefore does not cancel.
    thermal && (acc += _thermal_lookup(T_γ, s))

    acc / (τ_n * _FN)
end

# --- Background during BBN ---------------------------------------------------
#
# ρ and s for e± need the full Fermi-Dirac integrals: e± annihilation is *not*
# instantaneous, it happens right across the nucleosynthesis window, and the width of
# that transition is what sets T_ν/T_γ. A g_*(T) step function would put the
# annihilation in the wrong place.

"Energy density of e± (both charges, g=4) in units of MeV⁴, for T in MeV."
function _ρ_epm(T)
    T < 1e-4 && return 0.0
    m = _ME_MEV / T
    f, _ = quadgk(u -> begin
            E = sqrt(u^2 + m^2)
            u^2 * E / (exp(E) + 1)
        end, 0.0, 60 + 10m; rtol=1e-10)
    (2 / π^2) * T^4 * f
end

"Pressure of e± in MeV⁴."
function _P_epm(T)
    T < 1e-4 && return 0.0
    m = _ME_MEV / T
    f, _ = quadgk(u -> begin
            E = sqrt(u^2 + m^2)
            u^4 / (E * (exp(E) + 1))
        end, 0.0, 60 + 10m; rtol=1e-10)
    (2 / (3π^2)) * T^4 * f
end

_ρ_γ(T) = (π^2 / 15) * T^4
_s_γe(T) = (4 / 3) * _ρ_γ(T) / T + (_ρ_epm(T) + _P_epm(T)) / T

"""
    T_nu_over_T_gamma(T_γ)

Neutrino-to-photon temperature ratio, from entropy conservation in the γ-e± plasma.

a·T_ν is constant (neutrinos are decoupled and free), and s_{γe}·a³ is constant
(photons and e± stay coupled), so T_ν/T_γ = [ŝ(T_γ)/ŝ(∞)]^{1/3} with ŝ = s/T³.
Taking T → 0 gives back (4/11)^{1/3} without it ever being typed in.
"""
function T_nu_over_T_gamma(T_γ)
    ŝ∞ = (2π^2 / 45) * (11 / 2)          # γ + e± while relativistic: g_s = 2 + 7/8·4
    (_s_γe(T_γ) / T_γ^3 / ŝ∞)^(1 / 3)
end

# --- The solve ---------------------------------------------------------------

"""Standard GR: H = √(8πG ρ/3), with ρ in g/cm³. The default `friedmann` hook."""
_gr_friedmann(ρ_gcm3) = sqrt(8π * 6.67430e-8 * ρ_gcm3 / 3)

"""No extra energy density beyond γ, e±, ν. The default `ρ_extra` hook."""
_zero_ρ(_T_γ) = 0.0

"""
    BBNSolution

Primordial abundances. `Y` holds n_i/n_b for each of [`NUCLIDES`](@ref).
"""
struct BBNSolution
    ω_b::Float64
    N_eff::Float64
    Y::Vector{Float64}
    Y_p::Float64        # ⁴He mass fraction, CLASS/sBBN convention (real atomic masses)
    Y_p_nucleon::Float64 # 4·Y(⁴He), the convention BBN papers usually headline
    DH::Float64         # D/H by number
    He3H::Float64
    Li7H::Float64       # (⁷Li + ⁷Be)/H -- ⁷Be electron-captures to ⁷Li much later
end

"""
    bbn(; ω_b = 0.02242, ω_m = 0, N_eff = 3.044, τ_n = 878.4,
        T9_weak = 100, T9_start = 10, T9_end = 0.005)

Integrate the nuclear network and return the primordial abundances.

`τ_n` is the neutron lifetime in seconds (PDG 2024: 878.4 ± 0.5). It is the single
measured number the weak rates need, and Y_p is genuinely sensitive to it: ±0.5 s
moves Y_p by about ±0.0002.

The integration variable is ln T9 rather than time. Time never appears explicitly:
a·T_ν is constant, so H = -dln T_ν/dt gives dt/dln T_γ = -(dln T_ν/dln T_γ)/H, and the
whole thermal history follows from the entropy of the γ-e± plasma.
`ω_m` supplies the present physical matter density to the Friedmann rate;
`bbn(cosmology)` passes the cosmology's baryon-plus-CDM value automatically.
"""
function bbn(; ω_b=0.02242, ω_m=0.0, N_eff=3.044, τ_n=878.4,
    ρ_extra=_zero_ρ, friedmann=_gr_friedmann,
    network=:full,          # :full (63 reactions, 12 nuclides) or :small (12, 8)
    rates=:primat,          # rate compilation: :primat (consistent chain) or :parthenope core
    T9_weak=100.0, T9_start=10.0, T9_end=0.005,
    reltol=1e-10, abstol=1e-25, return_sol=false,
    _net=nothing)           # a pre-sampled network, for bbn_mc

    net = _net === nothing ? _network(network; rates) : _net
    sym = [_symmetry(r.reactants) for r in net]
    symp = [_symmetry(r.products) for r in net]

    K_MEV_PER_T9 = 8.617333262e-2 / 1000 * 1e9 / 1e6   # k_B·10⁹K in MeV = 0.0861733
    T_of_T9(T9) = T9 * 0.086173332621

    # Nucleon number density today, from ω_b. The subtlety: ω_b measures the baryon
    # *mass* density, and today most nucleons sit inside helium, whose mass per nucleon
    # is below m_H by the binding energy. Dividing by the atomic mass unit miscounts the
    # nucleons by 0.6%, which the network then carries as a 0.6% error in η -- and D/H
    # responds to η with a power of about -1.6. Divide by the actual mean mass per
    # nucleon of the post-BBN composition (H + He4, helium mass fraction ≈ 0.245; the
    # residual dependence of the mean on Y_p is a 1e-5 relative effect, far below η's
    # other uncertainties). Atomic masses, so the electrons ω_b also weighs are counted.
    ρ_c0_over_h2 = 1.87834e-29                  # g/cm³
    m_u = 1.66053907e-24                        # g
    Yp_today = 0.245
    m_per_nucleon = ((1 - Yp_today) * M_H_ATOM + Yp_today * M_HE4_ATOM / 4) * m_u
    n_b0 = ω_b * ρ_c0_over_h2 / m_per_nucleon   # cm⁻³, nucleon count
    T_ν0 = 2.7255 * 8.617333262e-11 * (4 / 11)^(1 / 3)   # MeV  (k_B·T_cmb·(4/11)^⅓)

    N_A = 6.02214076e23
    MeV4_to_g_cm3 = 2.3201e17 * 1e-3            # not used; kept explicit below instead

    # ρ_b(T) from n_b ∝ a⁻³, with a taken from the decoupling solve.
    #
    # The tempting shortcut is n_b ∝ T_ν³, on the grounds that a·T_ν is constant. But it is
    # only constant *after* the neutrinos stop being heated -- and the heating is exactly
    # what the decoupling solve is for. Using T_ν as a clock through the annihilation
    # quietly assumes the instantaneous limit one is trying to correct.
    function baryon_density(T_γ)                 # returns ρ_b in g/cm³
        a = scale_factor_at(T_γ)
        n_b0 * a^-3 * m_u
    end

    # Expansion rate. G in cgs; energy densities converted from MeV⁴ to g/cm³.
    # 1 MeV⁴ = (1 MeV)⁴/(ħc)³ /c² ... do it once, explicitly:
    #   ρ[g/cm³] = ρ[MeV⁴] · (1.602176634e-6 erg/MeV) / (197.3269804e-13 MeV·cm / MeV)³ / c²
    ħc_MeVcm = 197.3269804e-13
    erg_per_MeV = 1.602176634e-6
    c_cgs = 2.99792458e10
    MeV4_to_gcm3 = erg_per_MeV / ħc_MeVcm^3 / c_cgs^2
    G_cgs = 6.67430e-8

    # The expansion rate is the *only* way the rest of physics reaches BBN, and it is
    # deliberately left as a hook rather than a hardcoded Friedmann equation.
    #
    # Nucleosynthesis is a race between the weak rates and H: everything the universe
    # is made of, and every modification to gravity, enters Y_p through H(T) and
    # through nothing else. So `ρ_extra` takes any additional energy density (dark
    # radiation, a relic, whatever) in MeV⁴, and `friedmann` maps total ρ to H --
    # which is exactly where a modified gravity model substitutes its own H(ρ).
    # Default is standard GR.
    #
    # Matter is retained explicitly when the caller supplies ω_m. Its contribution is
    # very small during BBN, but carrying it makes the background bookkeeping exact and
    # avoids silently switching cosmologies at the nuclear-network boundary. Vacuum energy
    # or a modified expansion law belongs in `ρ_extra`/`friedmann`.
    # The neutrino sector comes from the *solved* decoupling history, not from a model of
    # it. `nu_temperatures(T_γ)` returns T_νe and T_νμ as computed by `neff()`: the actual
    # non-instantaneous decoupling, with the e± annihilation energy leaking into the
    # neutrinos exactly as fast as the collision integrals say it does.
    #
    # This replaces two approximations at once. There used to be an instantaneous-decoupling
    # T_ν/T_γ from entropy conservation, and on top of it a hand-built ramp taking N_eff
    # from 3.000 to 3.044 across the annihilation. The ramp reproduced the endpoints by
    # construction and the path between them by assumption -- and since freeze-out happens
    # *during* that path, the path is the part that matters. Now there is no ramp: the
    # temperatures are read off the trajectory.
    #
    # Anything beyond the three Standard-Model neutrinos (dark radiation, a light relic) is
    # a genuine extra species, present throughout and not part of this history.
    ΔN_extra = N_eff - 3.044

    function hubble(T_γ)                          # 1/s
        Tνe, Tνμ = nu_temperatures(T_γ)
        # ρ_ν from the actual per-flavour temperatures. ν_e runs hotter than ν_μ (it has a
        # charged-current coupling as well as a neutral one), and that split is a result of
        # the decoupling solve, not an input.
        ρν = (7 / 8) * (π^2 / 15) * (Tνe^4 + 2 * Tνμ^4)
        # Extra relativistic species, in the convention of the sBBN tables: ΔN = 1 is
        # one full species mimicking a decoupled neutrino, so it carries the *average*
        # per-flavour ν energy density -- (Tνe⁴ + 2Tνμ⁴)/3 -- not a third of one
        # flavour. Getting this factor wrong leaves ΔN = 0 untouched and quietly
        # triples out of the Y_p response to dark radiation.
        ρν += ΔN_extra * (7 / 8) * (π^2 / 15) * (Tνe^4 + 2 * Tνμ^4) / 3
        ρ_MeV4 = _ρ_γ(T_γ) + _ρ_epm(T_γ) + ρν + ρ_extra(T_γ)
        ρ_m_gcm3 = ω_m * ρ_c0_over_h2 / scale_factor_at(T_γ)^3
        friedmann(ρ_MeV4 * MeV4_to_gcm3 + ρ_m_gcm3)
    end

    # The clock: dt/dln T_γ = (dln a/dln T_γ)/H, with a(T_γ) from the decoupling solve.
    # It is tempting to write dt = -dln T_ν/H instead, on the grounds that a·T_ν is
    # constant -- but it is only constant *after* the heating stops. a·T_ν grows by 0.14%
    # through the annihilation, and a clock that ignores that undercounts the time
    # between freeze-out and the deuterium bottleneck. Fewer seconds means fewer neutron
    # decays means Y_p high, at the level of a part in a thousand -- the same
    # instantaneous-decoupling assumption the comment above n_b warns about, hiding in
    # the time variable instead of the density.
    function dlna_dlnTγ(T_γ)
        δ = 1e-4
        (log(scale_factor_at(T_γ * exp(δ))) - log(scale_factor_at(T_γ * exp(-δ)))) / (2δ)
    end

    # Everything above depends only on temperature, and the right-hand side is
    # evaluated thousands of times. The e± Fermi-Dirac integrals and the 200-point
    # Gauss rule for each weak rate are far too expensive to run per step -- and
    # quadgk inside an RHS is not differentiable, so a Rosenbrock solver cannot even
    # build its Jacobian through it. Tabulate once, on a grid in ln T9, exactly as
    # `BackgroundCache` does for the perturbations.
    lnT9s = collect(range(log(T9_end) - 0.1, log(T9_weak) + 0.1; length=900))
    _ρb = similar(lnT9s)
    _dtdlnT = similar(lnT9s)
    _Γnp = similar(lnT9s)
    _Γpn = similar(lnT9s)
    Threads.@threads for i in eachindex(lnT9s)
        Tγ = T_of_T9(exp(lnT9s[i]))
        # The n <-> p rates are charged-current: they see ν_e specifically, which runs
        # hotter than the flavour average, and hotter still than the instantaneous
        # entropy-conservation ratio. Feeding them anything colder than the solved T_νe
        # undercounts the n -> p conversions and leaves Y_p high.
        Tνe, _ = nu_temperatures(Tγ)
        _ρb[i] = baryon_density(Tγ)
        _dtdlnT[i] = dlna_dlnTγ(Tγ) / hubble(Tγ)
        _Γnp[i] = weak_rate(Tγ, Tνe, τ_n; direction=1)
        _Γpn[i] = weak_rate(Tγ, Tνe, τ_n; direction=-1)
    end
    ρb_of = linear_interpolation(lnT9s, log.(_ρb); extrapolation_bc=Line())
    dt_of = linear_interpolation(lnT9s, _dtdlnT; extrapolation_bc=Line())
    Γnp_of = linear_interpolation(lnT9s, log.(_Γnp); extrapolation_bc=Line())
    Γpn_of = linear_interpolation(lnT9s, log.(max.(_Γpn, 1e-300)); extrapolation_bc=Line())

    function rhs!(dY, Y, _, lnT9)
        T9 = exp(lnT9)
        ρ_b = exp(ρb_of(lnT9))
        dtdlnT = dt_of(lnT9)                    # negative: T falls as t grows

        fill!(dY, 0.0)

        # weak n <-> p
        Γnp = exp(Γnp_of(lnT9))
        Γpn = exp(Γpn_of(lnT9))
        flow = Γnp * Y[IDX.n] - Γpn * Y[IDX.p]
        dY[IDX.n] -= flow
        dY[IDX.p] += flow

        # thermonuclear
        @inbounds for (j, r) in enumerate(net)
            # Tables store N_A^{n-1}<σv>: cm³/mol/s for two-body, cm⁶/mol²/s for
            # three-body, so each extra reactant beyond the first costs one power
            # of ρ_b. Same bookkeeping on the reverse side, where the photon (not
            # in the products list) is free: a radiative capture's reverse is a
            # one-body photodisintegration with no density factor at all.
            λf = r.rate(T9)
            fwd = λf / sym[j]
            for _ in 2:length(r.reactants)
                fwd *= ρ_b
            end
            for i in r.reactants
                fwd *= Y[i]
            end

            revrate = λf * r.α * T9^r.β * exp(r.γ / T9)
            rev = revrate / symp[j]
            for _ in 2:length(r.products)
                rev *= ρ_b
            end
            for i in r.products
                rev *= Y[i]
            end

            net_flux = fwd - rev
            for i in r.reactants
                dY[i] -= net_flux
            end
            for i in r.products
                dY[i] += net_flux
            end
        end

        dY .*= dtdlnT
        nothing
    end

    # ---- Stage 1: the weak era, n <-> p only, T9_weak -> T9_start ----
    #
    # The nuclear network cannot be started at T9 = 10, even though that is where the
    # rate tables begin. At T9 = 10 the conversion rate is only Γ/H ≈ 1.6: the weak
    # interactions are *already freezing out*, so imposing exact weak equilibrium
    # there forces n/p below its true value and every neutron missing from the count
    # is two missing helium nucleons. It costs ~7% of Y_p and it looks entirely
    # plausible while doing it.
    #
    # Above T9 = 10 nothing but n <-> p happens anyway -- deuterium photodisintegration
    # is ~1e17/s there, so the bottleneck is absolute and no nucleus survives. So run
    # the weak pair alone from T9 = 100, where Γ/H ≈ 900 and equilibrium is genuinely
    # exact, and hand over to the network at T9 = 10. This is the same split PRIMAT
    # makes ("the high-temperature era integrates only n and p directly").
    function weak_rhs!(du, u, _, lnT9)
        Γnp = exp(Γnp_of(lnT9))
        Γpn = exp(Γpn_of(lnT9))
        flow = (Γnp * u[1] - Γpn * u[2]) * dt_of(lnT9)
        du[1] = -flow
        du[2] = +flow
        nothing
    end
    rnp0 = exp(-_Q_NP_MEV / T_of_T9(T9_weak))
    w0 = [rnp0 / (1 + rnp0), 1 / (1 + rnp0)]
    wsol = solve(ODEProblem(weak_rhs!, w0, (log(T9_weak), log(T9_start))),
        TRBDF2(autodiff=AutoFiniteDiff()); reltol, abstol=1e-14)

    Y0 = zeros(length(NUCLIDES))
    Y0[IDX.n], Y0[IDX.p] = wsol.u[end]

    # ---- Stage 2: the full network ----
    #
    # Seed deuterium at its quasi-static equilibrium rather than at zero. With Y_d = 0
    # the network must relax across the n+p <-> d balance on a timescale 1/λ_rev ~ 2e-15 s,
    # which in ln T9 is below machine epsilon: the solver cannot take a first step and
    # gives up, silently returning the initial condition -- a Y_p of exactly zero that
    # still conserves baryon number perfectly. Seeding the balance puts the system on
    # the slow manifold from the outset.
    r1 = net[1]
    λf1 = r1.rate(T9_start)
    λr1 = λf1 * r1.α * T9_start^r1.β * exp(r1.γ / T9_start)
    Y0[IDX.d] = exp(ρb_of(log(T9_start))) * λf1 * Y0[IDX.n] * Y0[IDX.p] / λr1

    prob = ODEProblem(rhs!, Y0, (log(T9_start), log(T9_end)))
    sol = solve(prob, TRBDF2(autodiff=AutoFiniteDiff()); reltol, abstol)
    return_sol && return sol
    sol.retcode == ReturnCode.Success ||
        @warn "BBN integration did not converge" sol.retcode
    Y = max.(sol.u[end], 0.0)

    yHe4 = Y[IDX.He4]
    yH = Y[IDX.p]
    # Mass fraction with real atomic masses -- the convention CLASS's sBBN table uses.
    Yp = yHe4 * M_HE4_ATOM / (yH * M_H_ATOM + yHe4 * M_HE4_ATOM)

    BBNSolution(ω_b, N_eff, Y, Yp, 4 * yHe4,
        Y[IDX.d] / yH, Y[IDX.He3] / yH, (Y[IDX.Li7] + Y[IDX.Be7]) / yH)
end

"""
    Y_He_bbn(ω_b, N_eff; ω_m = 0, ω_k = 0, ω_de = 0,
        w0 = -1, wa = 0, τ_n = 878.4)

Primordial helium mass fraction, for feeding straight into recombination.
This is what `cosmology()` calls when `Y_He` is not given.
"""
const _STANDARD_BBN_CACHE = Dict{NTuple{8,Float64},BBNSolution}()
const _STANDARD_BBN_CACHE_LOCK = ReentrantLock()

_copy_bbn(b::BBNSolution) = BBNSolution(
    b.ω_b, b.N_eff, copy(b.Y), b.Y_p, b.Y_p_nucleon, b.DH, b.He3H, b.Li7H,
)

function _cached_standard_bbn(ω_b, ω_m, N_eff, τ_n;
    ω_k=0.0, ω_de=0.0, w0=-1.0, wa=0.0)
    key = (Float64(ω_b), Float64(ω_m), Float64(N_eff), Float64(τ_n),
        Float64(ω_k), Float64(ω_de), Float64(w0), Float64(wa))
    stored = lock(_STANDARD_BBN_CACHE_LOCK) do
        get!(_STANDARD_BBN_CACHE, key) do
            ħc_MeVcm = 197.3269804e-13
            MeV4_to_gcm3 = 1.602176634e-6 / ħc_MeVcm^3 / (2.99792458e10)^2
            function composition_extra(Tγ)
                a = scale_factor_at(Tγ)
                de = key[6] * a^(-3 * (1 + key[7] + key[8])) *
                    exp(-3 * key[8] * (1 - a))
                1.87834e-29 * (key[5] / a^2 + de) / MeV4_to_gcm3
            end
            bbn(; ω_b=key[1], ω_m=key[2], N_eff=key[3], τ_n=key[4],
                ρ_extra=composition_extra)
        end
    end
    # BBNSolution contains a mutable abundance vector. Never expose the canonical
    # cached vector through a Cosmology: mutating one model must not alter the
    # primordial composition of a model constructed later with the same parameters.
    _copy_bbn(stored)
end

function Y_He_bbn(ω_b, N_eff; ω_m=0.0, ω_k=0.0, ω_de=0.0,
    w0=-1.0, wa=0.0, τ_n=878.4)
    _cached_standard_bbn(ω_b, ω_m, N_eff, τ_n; ω_k, ω_de, w0, wa).Y_p
end

"The full primordial abundance solution carried by `c`, or `nothing` for manual Yp."
primordial_bbn(c::Cosmology) = c.bbn_result

"""
    bbn(c::Cosmology; kwargs...)

Run nucleosynthesis for the cosmology `c`, closing the loop: ω_b and N_eff are read
off the species list rather than passed separately, so the helium that comes out is
the helium that *that* baryon density and radiation content actually make.

Any relativistic species beyond γ and the standard neutrinos -- a dark-radiation
component, a light relic -- shows up in the expansion rate during nucleosynthesis and
therefore in Y_p, which is the whole reason Y_p is a useful probe of N_eff. Passing a
`friedmann` hook is the entry point for modified gravity: BBN feels a modified H(ρ)
through the freeze-out of the weak rates and through nothing else.

Matter, curvature, and the cosmology's exact dark-energy `ρ(a)` are retained too;
their standard-ΛCDM contributions are tiny at BBN, but early-dark-energy models are
therefore handed to the network without a radiation-era substitution.
"""
function bbn(c::Cosmology; τ_n=878.4, kwargs...)
    ω_b = Ω_b(c) * c.h^2
    ω_m = (Ω_b(c) + Ω_c(c)) * c.h^2
    # N_eff from the actual neutrino content, not from a remembered input.
    N_eff = _N_eff_of(c)
    if haskey(kwargs, :ρ_extra)
        # An explicit hook takes precedence, as documented by the low-level API.
        return bbn(; ω_b, ω_m, N_eff, τ_n, kwargs...)
    end

    # Preserve the cosmology's exact curvature and dark-energy history during
    # BBN too. Curvature is an effective Friedmann density (possibly negative),
    # while every dark-energy subtype supplies its own exact ρ(a). Standard Λ is
    # fantastically small here, but an early-dark-energy model need not be.
    ħc_MeVcm = 197.3269804e-13
    erg_per_MeV = 1.602176634e-6
    c_cgs = 2.99792458e10
    MeV4_to_gcm3 = erg_per_MeV / ħc_MeVcm^3 / c_cgs^2
    ρc0_gcm3 = c.h^2 * 1.87834e-29
    function composition_extra(Tγ)
        a = scale_factor_at(Tγ)
        rel = sum((ρ_over_ρc0(s, a) for s in c.species
            if s isa Curvature || s isa AbstractDarkEnergy); init=0.0)
        ρc0_gcm3 * rel / MeV4_to_gcm3
    end
    bbn(; ω_b, ω_m, N_eff, τ_n, ρ_extra=composition_extra, kwargs...)
end

"""
    bbn_mc(; n = 400, στ_n = 0.5, seed = 0x5eed, kwargs...) -> NamedTuple

Monte-Carlo propagation of the nuclear-rate and neutron-lifetime uncertainties
through the network. Each draw displaces every forward rate by ξᵢ standard
deviations of its own tabulated 1σ envelope (λ = median·σ(T9)^ξᵢ, ξᵢ ~ N(0,1) --
the PArthENoPE/PRIMAT convention; the error columns in the rate files finally
earn their keep) and draws τ_n ~ N(τ_n, στ_n). Reverse rates follow each draw
automatically through detailed balance.

Returns means, standard deviations and 16/50/84 percentiles for Y_p, D/H, ³He/H
and ⁷Li/H. Draws run threaded; `kwargs` are passed through to `bbn`.
"""
function bbn_mc(; n=400, ω_b=0.02242, N_eff=3.044, τ_n=878.4, στ_n=0.5,
    network=:full, rates=:primat, seed=0x5eed, kwargs...)
    nrx = length(_network(network; rates))        # also warms the table cache
    out = Matrix{Float64}(undef, n, 4)
    Threads.@threads for i in 1:n
        rng = Random.Xoshiro(seed * 1000003 + i)
        ξ = randn(rng, nrx)
        τ = τ_n + στ_n * randn(rng)
        b = bbn(; ω_b, N_eff, τ_n=τ, network,
            _net=_network(network; rates, ξ), kwargs...)
        out[i, :] .= (b.Y_p, b.DH, b.He3H, b.Li7H)
    end
    q(v, p) = Statistics.quantile(v, p)
    stat(k) = (mean=Statistics.mean(out[:, k]), std=Statistics.std(out[:, k]),
        p16=q(out[:, k], 0.16), p50=q(out[:, k], 0.5), p84=q(out[:, k], 0.84))
    (Y_p=stat(1), DH=stat(2), He3H=stat(3), Li7H=stat(4), draws=out)
end

"""
Effective number of neutrino species implied by a cosmology's ν content.

Measured from the density ratio deep in the radiation era, not from Ω_ν today. Those
are wildly different things once a neutrino has mass: Ω_ν(today) is its *rest-mass*
density (Σm_ν/93.14 eV), which for Σm_ν = 0.06 eV is ~40x its relativistic
contribution. Reading N_eff off it gives 117 instead of 3.044.

At a = 1e-12 every neutrino, massive or not, is ultra-relativistic (T_ν ~ MeV >> m_ν),
which is precisely the epoch BBN cares about -- so the ratio there is the one that
matters, and it falls straight out of the species machinery.
"""
function _N_eff_of(c::Cosmology; a=1e-12)
    ργ = ρ_over_ρc0(get_species(c, Photons), a)
    ρν = 0.0
    for s in c.species
        (s isa MasslessNeutrinos || s isa MassiveNeutrinos) && (ρν += ρ_over_ρc0(s, a))
    end
    ρν / ργ / ((7 / 8) * (4 / 11)^(4 / 3))
end
