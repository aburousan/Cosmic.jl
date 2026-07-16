"""
Post-BBN trace-species energy deposition.

The primordial ³H and ⁷Be abundances are taken from the `BBNSolution` carried by
the cosmology. Tritium beta decay and hydrogenic ⁷Be recombination/electron
capture are evolved in the actual cosmological background. Only electromagnetic
energy is deposited: the tritium antineutrino and both ⁷Be neutrino lines are
explicitly excluded.

The ⁷Be kinetics and energies follow Khatri & Sunyaev, Astronomy Letters 37
(2011), arXiv:1009.3932, eqs. (1)--(9) and section 4.
"""

using OrdinaryDiffEq: ODEProblem, solve, Rodas5P, ReturnCode

const _DAY_S = 86400.0
const _YEAR_S = 365.25 * _DAY_S
const _LAMBDA_T = log(2.0) / (12.32 * _YEAR_S)
const _LAMBDA_BE3 = 0.5 * log(2.0) / (53.2 * _DAY_S)

# Evaluated mean beta-electron kinetic energy. The remaining beta endpoint
# energy is carried predominantly by the antineutrino and is not deposited.
const _E_T_BETA_EV = 5685.0

# ⁷Be electron-capture branches (Khatri-Sunyaev section 4): the neutrino takes
# 861.8 keV (ground branch) or 384.2 keV (10.44% excited branch) and is excluded. The
# electromagnetic/recoil products below thermalize in the plasma.
const _BE_EXCITED_BRANCH = 0.1044
const _E_BE_GAMMA_EV = 477_600.0
const _E_BE_RECOIL_GROUND_EV = 56.5
const _E_BE_RECOIL_EXCITED_EV = 11.2
const _E_BE_DEPOSIT_EV = (1 - _BE_EXCITED_BRANCH) * _E_BE_RECOIL_GROUND_EV +
    _BE_EXCITED_BRANCH * (_E_BE_GAMMA_EV + _E_BE_RECOIL_EXCITED_EV)

"Total hydrogenic ⁷Be⁴⁺ recombination coefficient, in m³/s."
function _alpha_be4(T)
    Z = 4.0
    a, b, c, d = 5.596, -0.6038, 0.3436, 0.4479
    t = (T / 1.0e4) / Z^2
    1.0e-13 * Z * a * t^b / (1 + c * t^d) * 1.0e-6
end

"⁷Be³⁺ photoionization coefficient fixed by detailed balance, in 1/s."
function _beta_be3(T)
    Z = 4.0
    χ = Z^2 * Constants.k_B_SI * Constants.T_of_wavenumber(Constants.L_H_ion)
    hP = 2π * Constants.ħ_SI
    saha = (2π * Constants.m_e_SI * Constants.k_B_SI * T)^1.5 / hP^3
    _alpha_be4(T) * saha * exp(-χ / (Constants.k_B_SI * T))
end

struct TraceNuclearHistory{ST,SB}
    tritium::ST
    beryllium::SB
    tritium_to_H::Float64
    be7_to_H::Float64
    zt_hi::Float64
    zt_lo::Float64
    zbe_hi::Float64
    zbe_lo::Float64
end

const _TRACE_HISTORY_CACHE = IdDict{Any,Any}()
const _TRACE_HISTORY_LOCK = ReentrantLock()

function _trace_temperature_electrons(c, rec, z)
    # The production recombination solution starts near z=5000. ⁷Be evolves an
    # order of magnitude earlier, where tight thermal coupling and full H/He
    # ionization are exact to exponentially small corrections.
    if z > 5000
        T = c.Tcmb * (1 + z)
        ne = n_H_of_z(c, z) * (1 + 2f_He(c))
    else
        T = T_matter(rec, z)
        ne = n_e(rec, z)
    end
    T, ne
end

"""
    trace_nuclear_history(c, rec)

Solve the post-BBN ³H and ⁷Be histories. Returns `nothing` when `c` was built
with a manual helium abundance and therefore carries no trace-nuclide solution.
"""
function trace_nuclear_history(c::Cosmology, rec::RecombinationSolution;
    zt_hi=2.0e7, zt_lo=1.0e3, zbe_hi=1.0e5, zbe_lo=1.0e4)
    b = primordial_bbn(c)
    b === nothing && return nothing

    lock(_TRACE_HISTORY_LOCK) do
        get!(_TRACE_HISTORY_CACHE, rec) do
            Hs(z) = H_SI(c, z)

            function trhs!(du, u, _, z)
                du[1] = _LAMBDA_T * u[1] / (Hs(z) * (1 + z))
            end
            st = solve(ODEProblem(trhs!, [1.0], (zt_hi, zt_lo)), Rodas5P();
                reltol=2e-11, abstol=1e-13)
            st.retcode == ReturnCode.Success || error("tritium history failed: $(st.retcode)")

            function brhs!(du, u, _, z)
                X4, X3, _ = u
                T, ne = _trace_temperature_electrons(c, rec, z)
                ar = _alpha_be4(T)
                bi = _beta_be3(T)
                den = Hs(z) * (1 + z)
                du[1] = (ne * X4 * ar - bi * X3) / den
                du[2] = (-ne * X4 * ar + (bi + _LAMBDA_BE3) * X3) / den
                du[3] = -_LAMBDA_BE3 * X3 / den
            end
            sb = solve(ODEProblem(brhs!, [1.0, 0.0, 0.0], (zbe_hi, zbe_lo)),
                Rodas5P(); reltol=2e-11, abstol=1e-14)
            sb.retcode == ReturnCode.Success || error("⁷Be history failed: $(sb.retcode)")

            yp = b.Y[IDX.p]
            TraceNuclearHistory(st, sb, b.Y[IDX.t] / yp, b.Y[IDX.Be7] / yp,
                zt_hi, zt_lo, zbe_hi, zbe_lo)
        end
    end
end

@inline function _tritium_fraction(h::TraceNuclearHistory, z)
    z >= h.zt_hi && return 1.0
    z <= h.zt_lo && return h.tritium(h.zt_lo)[1]
    clamp(h.tritium(z)[1], 0.0, 1.0)
end

@inline function _be_fractions(h::TraceNuclearHistory, z)
    z >= h.zbe_hi && return (1.0, 0.0, 0.0)
    z <= h.zbe_lo && return Tuple(clamp.(h.beryllium(h.zbe_lo), 0.0, 1.0))
    Tuple(clamp.(h.beryllium(z), 0.0, 1.0))
end

"""
    nuclear_deposition_rate(c, rec, z; history=nothing)

Electromagnetic power from primordial ³H and ⁷Be in J m⁻³ s⁻¹. Neutrino and
antineutrino energy is absent by construction.
"""
function nuclear_deposition_rate(c::Cosmology, rec::RecombinationSolution, z;
    history=nothing)
    h = history === nothing ? trace_nuclear_history(c, rec) : history
    h === nothing && return 0.0
    nH = n_H_of_z(c, z)
    Ft = _tritium_fraction(h, z)
    _, X3, _ = _be_fractions(h, z)
    Pt = nH * h.tritium_to_H * Ft * _LAMBDA_T * _E_T_BETA_EV * Constants.eV_SI
    Pbe = nH * h.be7_to_H * X3 * _LAMBDA_BE3 * _E_BE_DEPOSIT_EV * Constants.eV_SI
    Pt + Pbe
end

"Fractional positive heating magnitude (dQ/dz)/ρ_γ from trace nuclei."
function nuclear_heating_rate(c::Cosmology, rec::RecombinationSolution, z;
    history=nothing)
    P = nuclear_deposition_rate(c, rec, z; history)
    ρgamma = Constants.a_rad_SI * (c.Tcmb * (1 + z))^4
    P / (ρgamma * H_SI(c, z) * (1 + z))
end
