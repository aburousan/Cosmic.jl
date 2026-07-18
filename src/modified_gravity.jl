"""
Linear Horndeski gravity in the EFT-of-DE stable basis.

The α-parametrisation of Bellini & Sawicki (arXiv:1404.3713) describes every
linear-order Horndeski model with α_T = 0 (GW170817) through four functions of
time — the kineticity α_K, braiding α_B, Planck-mass run-rate α_M = d ln M*²/
d ln a — plus the background. Sampling the α's directly is unstable territory:
most choices violate the no-ghost/no-gradient conditions somewhere, and c_s²
computed from them suffers catastrophic cancellation (mochi_class,
arXiv:2407.11968, Fig. 1). The stable basis of Kennedy, Lombriser & Taylor
(arXiv:1804.04582) inverts the logic: the user supplies

    D_kin(x) > 0,  M*²(x) > 0,  c_s²(x) > 0,  α_B0 = α_B(z=0)

as functions of x = ln a, all three positivity conditions become input
constraints satisfied by construction, and the braiding follows by integrating
the Riccati equation that the exact c_s² expression (B&S eq. 3.13, α_T = 0)

    c_s² = −[(2−α_B)(Ḣ − ½H²α_B − H²α_M) − H·α̇_B + ρ̃_m + p̃_m] / (H² D_kin)

defines for α_B once c_s² is prescribed (overdots: coordinate time; tildes:
divided by M*²; units 3H² = ρ_tot). In d/d ln a form, integrated backward
from α_B(0) = α_B0:

    dα_B/dx = c_s² D_kin + (2−α_B)(d ln H/dx − ½α_B − α_M) + (ρ̃_m+p̃_m)/H²

The kineticity is then α_K = D_kin − (3/2)α_B². The background is Cosmic's
own (the designer approach: H(a) is whatever the species bag gives, and the
scalar's energy is defined by the Friedmann closure) — matching mochi_class's
default of fixing the expansion history and modifying only the perturbations.
"""

using OrdinaryDiffEq: ODEProblem, solve, Vern9
using Interpolations: linear_interpolation, Line

export HorndeskiEFT, HorndeskiFunctions, stable_basis_solve

"""
    HorndeskiEFT(; D_kin, M_star2, c_s2, α_B0 = 0.0)

Stable-basis specification of a linear Horndeski model (α_T = 0). `D_kin`,
`M_star2` and `c_s2` are callables of x = ln a returning the de-mixed kinetic
term D = α_K + 3α_B²/2, the effective Planck mass squared, and the scalar
sound speed squared; all three must be positive at every x (this IS the
no-ghost + no-gradient stability requirement, imposed by construction).
`α_B0` is the braiding today.
"""
struct HorndeskiEFT{FD,FM,FC}
    D_kin::FD
    M_star2::FM
    c_s2::FC
    α_B0::Float64
end

HorndeskiEFT(; D_kin, M_star2=x -> 1.0, c_s2=x -> 1.0, α_B0=0.0) =
    HorndeskiEFT(D_kin, M_star2, c_s2, float(α_B0))

"""
    HorndeskiFunctions

Solved α-functions on an x = ln a grid: callables `α_B(x)`, `α_K(x)`,
`α_M(x)`, `M_star2(x)`, their x-derivatives `α_Bx`, `α_Kx`, `α_Mx` (needed by
the scalar equation of motion, which carries α̇ terms), plus the grid and the
spec that generated them. Produced by [`stable_basis_solve`](@ref); consumed
by the perturbation layer. `α_Bx` is exact — it is the Riccati right-hand
side evaluated on the solution, not a numerical derivative.
"""
struct HorndeskiFunctions{I1,I2,I3,I4,I5,I6,FM,S}
    xs::Vector{Float64}
    α_B::I1
    α_K::I2
    α_M::I3
    α_Bx::I4
    α_Kx::I5
    α_Mx::I6
    M_star2::FM
    spec::S
    # activation point: the scalar is frozen (v_X = v̇_X = 0, GR equations)
    # for x < x_on. Before x_on the braiding dominates the kinetic term
    # (α_B² > D_kin, D_kin at its numerical floor): the true D_kin of a
    # late-time spec dies like Ω_DE(a) while the Riccati decaying mode of
    # α_B only falls like a, so the early-time system is singularly braided
    # and integrating it amplifies noise at rate ≫ H — the same reason
    # mochi_class activates modifications only below a threshold redshift.
    x_on::Float64
end

# d ln f/dx by 5-point central stencil. The background and the user's basis
# functions are smooth; h = 1e-3 puts the O(h⁴) truncation error near 1e-13
# without complex-step or AD requirements on the callables.
function _dln_dx(f, x; h=1e-3)
    (-log(f(x + 2h)) + 8log(f(x + h)) - 8log(f(x - h)) + log(f(x - 2h))) / (12h)
end

# (ρ_m + p_m)/H² in units 3H² = ρ_tot, matter = every species except dark
# energy: 3 Σ_m (1 + w_s) Ω_s(a). Curvature is geometry, not a fluid — the
# Newtonian-gauge system here is derived for flat space, so its presence is
# rejected in stable_basis_solve rather than summed.
function _matter_enthalpy_over_H2(c::Cosmology, a)
    E2 = zero(a)
    hm = zero(a)
    for s in c.species
        ρ = ρ_over_ρc0(s, a)
        E2 += ρ
        s isa AbstractDarkEnergy && continue
        hm += (1 + w(s, a)) * ρ
    end
    3hm / E2
end

"""
    stable_basis_solve(c, spec; x_min = log(1e-6), nx = 1200) -> HorndeskiFunctions

Integrate the braiding Riccati equation backward from `spec.α_B0` at x = 0 and
assemble the α-functions. Errors if the cosmology is curved, or if any of the
three stability inputs fails positivity on the grid (which can only happen
through a user callable that changes sign — the point of the basis is that
stability is an input).

Early times: for any spec whose D_kin, ΔM*² and α_M vanish as x → x_min the
source of the Riccati equation vanishes and α_B → 0 dynamically — the ΛCDM
limit is reached by construction rather than imposed by hand.
"""
function stable_basis_solve(c::Cosmology, spec::HorndeskiEFT;
    x_min=log(1e-6), nx=1200)

    any(s -> s isa Curvature, c.species) &&
        throw(ArgumentError("stable-basis Horndeski requires a flat cosmology"))

    dlnE_dx(x) = (-log(E(c, exp(x + 2e-3))) + 8log(E(c, exp(x + 1e-3))) -
                  8log(E(c, exp(x - 1e-3))) + log(E(c, exp(x - 2e-3)))) / 12e-3

    αM(x) = _dln_dx(spec.M_star2, x)

    function rhs(u, p, x)
        α_B = u
        a = exp(x)
        cs2 = spec.c_s2(x)
        D = spec.D_kin(x)
        M2 = spec.M_star2(x)
        (cs2 > 0 && D > 0 && M2 > 0) || error(
            "stability input violated at x=$x: c_s²=$cs2, D_kin=$D, M*²=$M2")
        cs2 * D + (2 - α_B) * (dlnE_dx(x) - α_B / 2 - αM(x)) +
        _matter_enthalpy_over_H2(c, a) / M2
    end
    sol = solve(ODEProblem(rhs, spec.α_B0, (0.0, x_min)), Vern9();
        reltol=1e-12, abstol=1e-14)

    xs = collect(range(x_min, 0.0; length=nx))
    αBv = [sol(x) for x in xs]
    αMv = [αM(x) for x in xs]
    αKv = [spec.D_kin(x) - 1.5 * sol(x)^2 for x in xs]
    # α_B' exactly from the ODE; D_kin' and α_M' by stencil of the (smooth)
    # user callables, α_K' = D_kin' − 3 α_B α_B' by the defining relation
    αBxv = [rhs(sol(x), nothing, x) for x in xs]
    Dx(x) = (-spec.D_kin(x + 2e-3) + 8 * spec.D_kin(x + 1e-3) -
             8 * spec.D_kin(x - 1e-3) + spec.D_kin(x - 2e-3)) / 12e-3
    αKxv = [Dx(x) - 3 * sol(x) * rhs(sol(x), nothing, x) for x in xs]
    αMxv = [(-αM(x + 2e-3) + 8αM(x + 1e-3) - 8αM(x - 1e-3) + αM(x - 2e-3)) /
            12e-3 for x in xs]
    itp(v) = linear_interpolation(xs, v, extrapolation_bc=Line())
    # Activation floor is RELATIVE to the spec's own peak kinetic term: the
    # scalar stays frozen (GR equations) until D_kin exceeds 1e-4 of its
    # maximum. An absolute floor (1e-8) let the scalar switch on while D_kin
    # was still ~3e-8, and the 1/D_kin in the second-order EOM then amplified
    # any α-function noise into a blow-up (2026-07-18 testset). The frozen
    # window carries no scalar dynamics, so raising the floor shifts nothing
    # physical — validated: the hi_class propto_omega φ ratio is unchanged.
    Dv = [spec.D_kin(x) for x in xs]
    Dfloor = max(1e-8, 1e-4 * maximum(Dv))
    i_on = findfirst(i -> Dv[i] > Dfloor && αBv[i]^2 <= Dv[i], eachindex(xs))
    i_on === nothing && error(
        "stable-basis spec never reaches a dynamically meaningful regime " *
        "(needs D_kin > $(Dfloor) with α_B² ≤ D_kin somewhere)")
    HorndeskiFunctions(xs, itp(αBv), itp(αKv), itp(αMv),
        itp(αBxv), itp(αKxv), itp(αMxv), spec.M_star2, spec, xs[i_on])
end
