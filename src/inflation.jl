"""
Single-field inflation: the primordial spectra are *computed* from a potential
V(φ) instead of being posited as a power law. The background inflaton trajectory
and the Mukhanov–Sasaki mode functions are integrated exactly (no slow-roll
truncation), so any feature the potential imprints — a step, a bump, a brief
departure from slow roll — shows up in Δ²_ℛ(k) with the right shape.

Everything is in reduced-Planck units M_p = 1 and in e-folds N (dN = H dt).

Background (state φ, φ_N ≡ dφ/dN):

    H² = V / (3 − ε₁),   ε₁ = ½ φ_N²,   φ_NN = −(3 − ε₁) φ_N − V_φ/H²,

integrated from the initial field value to the end of inflation ε₁ = 1.

Curvature perturbation ℛ_k (comoving), exact Mukhanov–Sasaki in e-folds:

    ℛ_k'' + (3 − ε₁ + d ln ε₁/dN) ℛ_k' + (k/aH)² ℛ_k = 0,

with d ln ε₁/dN = 2 φ_NN/φ_N. Each mode starts in the Bunch–Davies vacuum deep
inside the horizon (k/aH = `sub_horizon`), v = e^{−ikτ}/√(2k), ℛ = v/z with
z = a√(2ε₁), and is evolved until it freezes on super-horizon scales. Then

    Δ²_ℛ(k) = k³/(2π²) |ℛ_k|²,   Δ²_t(k) = 4 · k³/(2π²) |h_k|² (both polarisations)

where the tensor mode obeys the same equation with z_t = a. The comoving k of
inflation is mapped to physical k via horizon crossing k = aH, with the pivot
placed `N_pivot` e-folds before the end of inflation; only the *shape* and
amplitude are physical, so this fixes the k-axis the way any inflationary
prediction does.

Feed the result straight into a cosmology: `cosmology(primordial = spectrum)`.
"""

using OrdinaryDiffEq: ODEProblem, solve, Vern9, Tsit5, ContinuousCallback, terminate!
using Interpolations: linear_interpolation, Line

# --- background trajectory ---------------------------------------------------

struct InflationBackground{IφN,Iε,IφNN,IaH}
    N_end::Float64
    φ_N::IφN          # φ_N(N)
    ε1::Iε            # ε₁(N)
    φNN_over_φN::IφNN # (φ_NN/φ_N)(N)
    aH::IaH           # a(N) H(N) = e^N H(N)
end

"""
    inflation_background(V, dV; φ0, φ0_dot=nothing, N_max=400, npts=6000)

Integrate the exact single-field background for potential `V` (with derivative
`dV`), starting from field value `φ0` on the slow-roll attractor (initial
φ_N = −V_φ/V unless `φ0_dot` is given), until inflation ends at ε₁ = 1.
"""
function inflation_background(V, dV; φ0, φ0_dot=nothing, N_max=400.0, npts=6000)
    φN0 = φ0_dot === nothing ? -dV(φ0) / V(φ0) : φ0_dot
    function rhs!(du, u, p, N)
        φ, φN = u
        ε1 = 0.5 * φN^2
        H2 = V(φ) / (3 - ε1)
        du[1] = φN
        du[2] = -(3 - ε1) * φN - dV(φ) / H2
    end
    # stop when ε₁ = ½φ_N² reaches 1
    cb = ContinuousCallback((u, N, integ) -> 0.5 * u[2]^2 - 1.0,
        (integ) -> terminate!(integ))
    prob = ODEProblem(rhs!, [float(φ0), φN0], (0.0, N_max))
    sol = solve(prob, Tsit5(); reltol=1e-10, abstol=1e-12, callback=cb)
    N_end = sol.t[end]

    Ns = collect(range(0.0, N_end; length=npts))
    φNv = similar(Ns)
    ε1v = similar(Ns)
    ratio = similar(Ns)
    aHv = similar(Ns)
    for (i, N) in enumerate(Ns)
        φ, φN = sol(N)
        ε1 = 0.5 * φN^2
        H2 = V(φ) / (3 - ε1)
        H = sqrt(H2)
        φNN = -(3 - ε1) * φN - dV(φ) / H2
        φNv[i] = φN
        ε1v[i] = ε1
        ratio[i] = φNN / φN
        aHv[i] = exp(N) * H
    end
    itp(y) = linear_interpolation(Ns, y; extrapolation_bc=Line())
    InflationBackground(N_end, itp(φNv), itp(ε1v), itp(ratio), itp(aHv))
end

# --- mode evolution ----------------------------------------------------------

# Evolve one comoving k. `scalar=true` uses z = a√(2ε₁) (curvature ℛ); false uses
# z_t = a (tensor h). State is [Re, Im, Re', Im'] of the perturbation in e-folds.
function _evolve_mode(bg::InflationBackground, k; scalar=true, sub_horizon=50.0,
    reltol=1e-9, abstol=1e-12, cs=nothing, dlncs=nothing)
    # start a few e-folds inside the horizon: c_s k/aH = sub_horizon (the sound
    # horizon is what matters when c_s ≠ 1)
    cs_of = cs === nothing ? (N -> 1.0) : cs
    # d ln c_s/dN: user-supplied, or a relative-step finite difference
    dlncs_of = dlncs !== nothing ? dlncs :
               cs === nothing ? (N -> 0.0) :
               (N -> (log(cs_of(N + 1e-4)) - log(cs_of(max(N - 1e-4, 0.0)))) /
                     (N < 1e-4 ? 1e-4 : 2e-4))
    f_ic(N) = cs_of(N) * bg.aH(N) - k / sub_horizon
    Nic = _bisect(f_ic, 0.0, bg.N_end)
    aH = bg.aH(Nic)
    ε1 = bg.ε1(Nic)
    a = exp(Nic)
    cs0 = cs_of(Nic)
    if scalar
        # EFT quadratic action: S = ∫ z²(ζ'² − c_s²(∂ζ)²) with z = a√(2ε)/c_s,
        # so z'/z = 1 + ε_N/(2ε) − c_s,N/c_s and the Bunch–Davies amplitude is
        # |ζ| = 1/(z√(2 c_s k)).
        z = a * sqrt(2 * ε1) / cs0
        zratio = 1 + bg.φNN_over_φN(Nic) - dlncs_of(Nic)
        amp = 1 / (z * sqrt(2 * cs0 * k))
    else
        amp = 1 / (a * sqrt(k))                   # z_t = a; h = √2 v/(a) with M_p=1
        zratio = 1.0
    end
    # Bunch–Davies: X = amp, X' = X(−(z'/z) − i c_s k/aH)
    u0 = [amp, 0.0, -amp * zratio, -amp * (cs0 * k / aH)]

    function rhs!(du, u, p, N)
        ε = bg.ε1(N)
        coef = scalar ? (3 - ε + 2 * bg.φNN_over_φN(N) - 2 * dlncs_of(N)) :
               (3 - ε)
        csN = scalar ? cs_of(N) : 1.0
        w2 = (csN * k / bg.aH(N))^2
        du[1] = u[3]
        du[2] = u[4]
        du[3] = -coef * u[3] - w2 * u[1]
        du[4] = -coef * u[4] - w2 * u[2]
    end
    prob = ODEProblem(rhs!, u0, (Nic, bg.N_end))
    sol = solve(prob, Vern9(); reltol, abstol)
    ue = sol.u[end]
    mag2 = ue[1]^2 + ue[2]^2
    pref = scalar ? 1.0 : 4.0                      # tensor: 4 = 2 polarisations × norm
    pref * k^3 / (2π^2) * mag2
end

# tiny bisection on a monotonic increasing f
function _bisect(f, lo, hi; tol=1e-10, maxit=200)
    flo = f(lo)
    fhi = f(hi)
    flo * fhi > 0 && return flo < 0 ? hi : lo     # clamp if out of bracket
    for _ in 1:maxit
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        (abs(fm) < tol || hi - lo < tol) && return mid
        (flo * fm <= 0) ? (hi = mid; fhi = fm) : (lo = mid; flo = fm)
    end
    0.5 * (lo + hi)
end

# --- the spectrum object -----------------------------------------------------

"""
    InflatonSpectrum

Computed primordial spectrum. Callable: `s(k)` returns the dimensionless curvature
power Δ²_ℛ(k) at physical wavenumber `k` (1/Mpc); `tensor_power(s, k)` returns
Δ²_t(k). Derived summary numbers: `n_s`, `r`, `A_s` at the pivot.
"""
struct InflatonSpectrum{IR,IT}
    k_pivot::Float64
    lnPR::IR          # ln Δ²_ℛ as a function of ln(k/k_pivot)
    lnPT::IT
    n_s::Float64
    r::Float64
    A_s::Float64
end

(s::InflatonSpectrum)(k) = exp(s.lnPR(log(k / s.k_pivot)))
tensor_power(s::InflatonSpectrum, k) = exp(s.lnPT(log(k / s.k_pivot)))

"""
    inflaton_spectrum(V, dV; φ0, N_pivot=55, k_pivot=0.05,
                      decades_below=4, decades_above=4, per_decade=20)

Compute the primordial curvature and tensor spectra from potential `V` (Planck
units, `dV` its derivative), starting the background at field value `φ0`. The
pivot scale `k_pivot` (1/Mpc) is placed `N_pivot` e-folds before the end of
inflation. Returns an [`InflatonSpectrum`](@ref) suitable for
`cosmology(primordial = ...)`.
"""
function inflaton_spectrum(V, dV; φ0, N_pivot=55.0, k_pivot=0.05,
    decades_below=4.0, decades_above=4.0, per_decade=20, sub_horizon=50.0,
    cs=nothing, dlncs=nothing)
    bg = inflation_background(V, dV; φ0)

    # EFT sound speed: c_s(N) on the background's own e-fold axis (N = 0 at
    # φ0, inflation ends at bg.N_end). Superluminal propagation is a causality
    # red flag referees will catch — refuse it outright.
    if cs !== nothing
        cs_max = maximum(cs(N) for N in range(0.0, bg.N_end; length=2000))
        cs_max ≤ 1 + 1e-12 || error(
            "c_s(N) reaches $(round(cs_max, digits=4)) > 1: superluminal sound " *
            "speed violates causality (a causality violation)")
        cs_min = minimum(cs(N) for N in range(0.0, bg.N_end; length=2000))
        cs_min > 0 || error("c_s(N) must be positive (min = $(cs_min))")
        cs_min < 0.02 && @warn "c_s(N) drops to $(round(cs_min, sigdigits=2)) — " *
                              "strong-coupling risk in the EFT" maxlog = 1
    end
    # The largest scale (smallest k) exits the horizon decades_below·ln10 e-folds
    # before the pivot and must still start sub-horizon (a further ln(sub_horizon)
    # e-folds earlier), so inflation must last at least this long before its end.
    need = N_pivot + decades_below * log(10) + log(sub_horizon) + 1
    bg.N_end ≥ need || error(
        "inflation lasts $(round(bg.N_end,digits=1)) e-folds, but the requested " *
        "k-range needs ≥ $(round(need,digits=1)) e-folds before the end. " *
        "Start at a larger |φ0|, or reduce decades_below / N_pivot.")

    # comoving k that exits N_pivot e-folds before the end == physical k_pivot
    # comoving pivot: the mode whose SOUND horizon exit (c_s k = aH) happens
    # N_pivot e-folds before the end — reduces to k = aH for c_s = 1
    N_piv = bg.N_end - N_pivot
    k_piv_inf = (cs === nothing ? 1.0 : cs(N_piv)) \ bg.aH(N_piv)
    nlo = round(Int, decades_below * per_decade)
    nhi = round(Int, decades_above * per_decade)
    ln_rel = collect(range(-decades_below * log(10), decades_above * log(10);
        length=nlo + nhi + 1))
    lnPR = similar(ln_rel)
    lnPT = similar(ln_rel)
    for (i, lr) in enumerate(ln_rel)
        k_inf = k_piv_inf * exp(lr)
        lnPR[i] = log(_evolve_mode(bg, k_inf; scalar=true, sub_horizon, cs, dlncs))
        lnPT[i] = log(_evolve_mode(bg, k_inf; scalar=false, sub_horizon))
    end
    PRitp = linear_interpolation(ln_rel, lnPR; extrapolation_bc=Line())
    PTitp = linear_interpolation(ln_rel, lnPT; extrapolation_bc=Line())

    # summary numbers at the pivot (central finite difference in ln k)
    h = 0.05
    ns = 1 + (PRitp(h) - PRitp(-h)) / (2h)
    A_s = exp(PRitp(0.0))
    r = exp(PTitp(0.0)) / A_s
    InflatonSpectrum(k_pivot, PRitp, PTitp, ns, r, A_s)
end
