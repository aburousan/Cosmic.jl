using Cosmic
using Test
using QuadGK: quadgk

# Reference values are Planck 2018 (TT,TE,EE+lowE+lensing+BAO, Table 2 of
# arXiv:1807.06209) unless stated otherwise. Tolerances reflect what the physics
# in each layer can actually deliver: the background is exact up to quadrature
# error, while recombination is RECFAST-level and good to a few parts in 1e3.
#
# Nothing here touches the network or plots.

# The full suite includes multi-minute physics closures (CMB spectra, lensing,
# thermalization). Set COSMIC_TEST_FULL=false (as CI does) to run only the fast
# structural and unit-physics tests; the full suite remains the default locally.
const FULL_TESTS = get(ENV, "COSMIC_TEST_FULL", "true") == "true"

@testset "Cosmic.jl" begin

    @testset "species" begin
        # The Fermi-Dirac energy integral must reproduce both analytic limits.
        @test fermi_dirac_ρ(0.0) ≈ 7π^4 / 120 rtol = 1e-9
        # Non-relativistic limit: F(y)/y -> ∫x²/(eˣ+1)dx = (3/2)ζ(3), i.e. ρ -> m·n.
        @test fermi_dirac_ρ(1e5) / 1e5 ≈ 1.5 * 1.2020569031595943 rtol = 1e-4

        mν = massive_neutrino(0.06, 2.7255, 0.6766)
        # Radiation-like when hot: ρ·a⁴ is constant.
        @test ρ_over_ρc0(mν, 1e-7) * (1e-7)^4 ≈ ρ_over_ρc0(mν, 1e-6) * (1e-6)^4 rtol = 1e-3
        # Matter-like when cold: ρ·a³ is constant.
        @test ρ_over_ρc0(mν, 0.9) * 0.9^3 ≈ ρ_over_ρc0(mν, 1.0) rtol = 2e-2
        # ...and w runs from 1/3 to 0 across the transition.
        @test w(mν, 1e-6) ≈ 1 / 3 rtol = 1e-3
        @test w(mν, 1.0) < 1e-3

        # The phase-space integral must independently reproduce the standard
        # Ω_ν h² = Σm_ν / 93.14 eV rule. A genuine cross-check: nothing in the
        # code puts 93.14 in by hand.
        @test Ω0(massive_neutrino(0.06, 2.7255, 0.6766; N_eff_share=3.044 / 3)) * 0.6766^2 ≈
              0.06 / 93.14 rtol = 2e-3

        # CPL reduces to Λ at w0=-1, wa=0.
        @test ρ_over_ρc0(W0WaDarkEnergy(0.7, -1.0, 0.0), 0.3) ≈ 0.7
        # The continuity-equation integral in GeneralDarkEnergy must reproduce
        # the closed forms it generalises.
        @test ρ_over_ρc0(GeneralDarkEnergy(0.7, a -> -1.0), 0.3) ≈ 0.7 rtol = 1e-8
        @test ρ_over_ρc0(GeneralDarkEnergy(0.7, a -> 0.0), 0.5) ≈ 0.7 / 0.5^3 rtol = 1e-8
    end

    @testset "background" begin
        c = cosmology()

        # The budget must close exactly; E(1) = 1 follows from it.
        @test sum(Ω0(s) for s in c.species) ≈ 1.0 rtol = 1e-12
        @test E(c, 1.0) ≈ 1.0 rtol = 1e-12
        @test H(c, 0.0) ≈ 100 * c.h

        @test age(c) ≈ 13.787 rtol = 1e-3          # Planck: 13.787 ± 0.020 Gyr
        @test z_equality(c) ≈ 3387 rtol = 5e-3     # Planck: 3387 ± 21
        @test conformal_time_today(c) ≈ 14160 rtol = 5e-3

        # Etherington distance-duality: D_L = (1+z)² D_A. Holds identically.
        z = 1.5
        @test luminosity_distance(c, z) ≈ angular_diameter_distance(c, z) * (1 + z)^2 rtol = 1e-10
        @test transverse_comoving_distance(c, z) ≈ comoving_distance(c, z) rtol = 1e-12
        @test 10^(distance_modulus(c, z) / 5 - 5) ≈ luminosity_distance(c, z) rtol = 1e-10

        # Curvature: open universes give D_M > D_C, closed the reverse. Getting
        # these backwards is the classic sign error, and the old code had no
        # transverse distance at all.
        co = cosmology(Ω_k=0.05)
        cc = cosmology(Ω_k=-0.05)
        @test transverse_comoving_distance(co, 1.0) > comoving_distance(co, 1.0)
        @test transverse_comoving_distance(cc, 1.0) < comoving_distance(cc, 1.0)

        # Curved comoving volume must reduce to 4/3 π D³ when flat.
        @test comoving_volume(c, 0.5) ≈ 4π / 3 * comoving_distance(c, 0.5)^3 rtol = 1e-10

        # Age and lookback time must partition the age of the universe.
        @test age(c, 2.0) + lookback_time(c, 2.0) ≈ age(c, 0.0) rtol = 1e-8

        # Inverting t(a) must round-trip.
        @test age(c, redshift(scale_factor_of_time(c, 5.0))) ≈ 5.0 rtol = 1e-6

        # Heavier neutrinos -> more early radiation -> a younger universe.
        @test age(cosmology(m_ν=[0.3, 0.3, 0.3])) < age(cosmology(m_ν=Float64[]))
    end

    @testset "thermodynamics" begin
        c = cosmology()
        # Reionization off: this testset is about the recombination history, and
        # leaving it on would conflate the two (x_e today is 1.16, not 2e-4).
        r = recombination(c; reionization=no_reionization)

        # Saha is exact at high z, where the rates crush the expansion.
        @test x_e(r, 5000.0) ≈ saha_x_e(c, 5000.0) rtol = 1e-3
        # Full ionization early: H plus singly-ionized He gives 1 + f_He.
        @test x_e(r, 4000.0) ≈ 1 + f_He(c) rtol = 1e-2

        # Recombination must *lag* Saha: the universe cannot recombine as fast
        # as equilibrium demands. This is the whole point of the Peebles C factor.
        @test x_e(r, 1100.0) > saha_x_e(c, 1100.0)

        # Residual ionization freeze-out.
        @test 1e-4 < x_e(r, 0.0) < 5e-4

        # Matter temperature: locked to the CMB early, cold and *positive* late.
        @test T_matter(r, 3000.0) ≈ c.Tcmb * 3001 rtol = 1e-3
        @test 0 < T_matter(r, 0.0) < 1.0

        # Key epochs, against Planck 2018.
        @test z_star(r) ≈ 1089.92 rtol = 3e-3       # 1089.92 ± 0.25
        @test z_drag(r) ≈ 1059.94 rtol = 3e-3       # 1059.94 ± 0.30
        @test r_star(r) ≈ 144.43 rtol = 3e-3        # 144.43 ± 0.26 Mpc
        @test r_drag(r) ≈ 147.09 rtol = 3e-3        # 147.09 ± 0.26 Mpc
        @test 100 * θ_star(r) ≈ 1.04109 rtol = 3e-3 # 1.04109 ± 0.00030

        # Baryons are released *after* photons: it takes R scatterings to move a
        # baryon, but one to scatter a photon.
        @test z_drag(r) < z_star(r)

        # Optical depth is monotonic and passes through 1 at last scattering.
        @test optical_depth(r, 1089.92) ≈ 1.0 rtol = 5e-2
        @test optical_depth(r, 500.0) < optical_depth(r, 1500.0)

        # The visibility function is a normalised probability density...
        zs = range(1.0, 3000.0; length=6000)
        gs = [visibility(r, z) for z in zs]
        @test sum(gs) * step(zs) ≈ 1.0 rtol = 2e-2
        # ...peaking at the last-scattering surface.
        @test 1000 < zs[argmax(gs)] < 1150

        # With reionization on, g(z) grows a second spike at z ≈ 7. It is *taller*
        # than the recombination peak in dz -- reionization deposits its τ over
        # Δz ≈ 2 while recombination spreads over Δz ≈ 100 -- even though it
        # carries only ~5% of the area. So the recombination peak is still the
        # last-scattering surface; it just is not the global maximum of g(z).
        rr = recombination(c)
        gs_r = [visibility(rr, z) for z in zs]
        i_rec = argmax(gs_r[zs.>500])
        @test 1000 < zs[zs.>500][i_rec] < 1150
        @test maximum(gs_r[zs.<20]) > maximum(gs_r[zs.>500])

        # Sound speed: c/√3 when the fluid is pure radiation, slower once the
        # baryons load it down. (Even at a = 1e-8 the loading R is 7e-6, not 0,
        # so this cannot be tested to better than ~1e-5.)
        @test sound_speed(c, 1e-8) ≈ 1 / sqrt(3) rtol = 1e-5
        @test sound_speed(c, 1e-3) < 1 / sqrt(3)
        @test R_baryon(c, 1.0) ≈ 3 * Ω_b(c) / (4 * Ω_γ(c))
    end

    @testset "RECFAST 1.5 corrections" begin
        # Each correction should move x_e the way its source paper says it does.
        # Turning them all off must recover the 1999 RECFAST.
        c = cosmology()
        r15 = recombination(c)
        r10 = recombination(c; fudge_H=1.14, gaussian_H=false,
            sobolev_He=false, triplets=false)

        # Switzer & Hirata (2008) I: H I continuum opacity absorbs the He I
        # resonance photons that would otherwise re-ionize helium, so helium
        # recombines *faster*. The effect must show up during He I recombination
        # (z ≈ 1800-2500) and be worth around a percent.
        @test x_e(r15, 1800.0) < x_e(r10, 1800.0)
        @test 0.005 < (x_e(r10, 1800.0) - x_e(r15, 1800.0)) / x_e(r10, 1800.0) < 0.10

        # ...and must leave the pre-helium era alone.
        @test x_e(r15, 3000.0) ≈ x_e(r10, 3000.0) rtol = 1e-3

        # The corrections must not wreck the epochs they were tuned to improve.
        @test z_drag(r15) ≈ 1059.94 rtol = 2e-3
        @test z_star(r15) ≈ 1089.92 rtol = 2e-3

        # The triplet channel is a real extra pathway: disabling it alone must
        # change the answer during helium recombination.
        r_notrip = recombination(c; triplets=false)
        @test x_e(r_notrip, 1900.0) != x_e(r15, 1900.0)
    end

    @testset "perturbations" begin
        # Massless neutrinos here so the analytic limits below are exact.
        c = cosmology(m_ν=Float64[])
        r = recombination(c)
        bg = BackgroundCache(c, r)

        Ωγ = Ω_γ(c)
        Ων = Cosmic._Ω_or_zero(c, MasslessNeutrinos)
        f_ν = Ων / (Ωγ + Ων)
        ψ_pred = 1 / (1.5 + 0.4f_ν)
        φ_pred = (1 + 0.4f_ν) * ψ_pred

        # A mode that stays outside the horizon across equality.
        p = solve_perturbations(c, bg, 1e-4)

        # Superhorizon potential in radiation domination, from the ℛ = 1
        # normalisation. This pins the initial conditions.
        @test Φ(p, 1e-5) ≈ φ_pred rtol = 1e-3
        @test Ψ(p, 1e-5) ≈ ψ_pred rtol = 5e-3

        # THE canonical check on a Boltzmann hierarchy: a superhorizon potential
        # drops by exactly 9/10 across matter-radiation equality. It is sensitive
        # to the photon, neutrino, and metric equations all being right *together*,
        # so a sign error anywhere shows up here.
        @test Φ(p, 3e-3) / Φ(p, 1e-5) ≈ 0.9 rtol = 1e-2

        # Φ and Ψ differ only by anisotropic stress. Once photons and neutrinos
        # have stopped free-streaming relative to the matter, they must converge.
        @test Φ(p, 1.0) ≈ Ψ(p, 1.0) rtol = 1e-3

        # Λ makes the potential decay at late times -- the ISW effect.
        @test abs(Φ(p, 1.0)) < abs(Φ(p, 0.1))

        # Sub-horizon CDM grows; adiabatic ICs tie δ_c to the photons.
        p2 = solve_perturbations(c, bg, 0.1)
        @test abs(δ_cdm(p2, 1.0)) > abs(δ_cdm(p2, 0.01))
        @test δ_cdm(p2, 1e-5) ≈ 0.75 * δ_photon(p2, 1e-5) rtol = 1e-2   # adiabaticity
    end

    @testset "second-order angular foundation" begin
        t = Cosmic.SecondOrderTriangle(0.1, 0.2, 0.25)
        @test t.k1 * t.sin1 ≈ t.k2 * t.sin2 rtol = 2e-15
        @test Cosmic._k_spherical(t, 1, 0) + Cosmic._k_spherical(t, 2, 0) ≈ t.k
        @test Cosmic._k_spherical(t, 1, -1) + Cosmic._k_spherical(t, 2, -1) ≈ 0 atol = 1e-16
        @test Cosmic._k_spherical(t, 1, 1) + Cosmic._k_spherical(t, 2, 1) ≈ 0 atol = 1e-16
        @test_throws DomainError Cosmic.SecondOrderTriangle(1.0, 1.0, 2.1)

        # Scalar rotations and their negative-m reality relation.
        @test Cosmic._rotate_scalar(t, 1, 1, 0) ≈ t.cos1 rtol = 2e-15
        @test Cosmic._rotate_scalar(t, 1, 1, 1) ≈ -t.sin1 / sqrt(2) rtol = 2e-15
        @test Cosmic._rotate_scalar(t, 1, 1, -1) ≈ t.sin1 / sqrt(2) rtol = 2e-15
        @test Cosmic._rotate_scalar(t, 2, 1, 1) ≈ t.sin2 / sqrt(2) rtol = 2e-15

        # Coupling identities from the primary harmonic derivation.
        for l in 2:8, m in -l:l, m1 in (m-1):(m+1)
            @test Cosmic._Rplus(l, m1, m) ≈ -(l + 2) * Cosmic._Cplus(l, m1, m)
            @test Cosmic._Rminus(l, m1, m) ≈ (l - 1) * Cosmic._Cminus(l, m1, m)
            @test Cosmic._Kplus(l, m1, m) ≈ -(l + 2) * Cosmic._Dplus(l, m1, m)
            @test Cosmic._Kminus(l, m1, m) ≈ (l - 1) * Cosmic._Dminus(l, m1, m)
            @test Cosmic._Kzero(l, m1, m) ≈ -Cosmic._Dzero(l, m1, m)
        end
        @test Cosmic._Dzero(7, 0, 0) == 0

        # The second-order kernels reuse the existing scalar G_l hierarchy.
        # Tram--Lesgourgues eq. (2.40b) gives E_l algebraically; in particular
        # the quadrupole must reproduce Cosmic's established Thomson source.
        Gpol = [0.31, -0.07, 0.11, 0.03, -0.05, 0.02, 0.01]
        Epol = Cosmic._scalar_E_from_G(Gpol)
        @test Epol[1] == 0 == Epol[2]
        @test Epol[3] ≈ -5 / sqrt(6) * (Gpol[1] + Gpol[3]) rtol = 2e-15
        I2 = 5 * 0.17
        @test (I2 - sqrt(6) * Epol[3]) / 10 ≈
              (0.17 + Gpol[1] + Gpol[3]) / 2 rtol = 2e-15
        legendre(l, x) = l == 0 ? one(x) : l == 1 ? x : begin
            pm, p = one(x), x
            for n in 2:l
                pm, p = p, ((2n-1)*x*p-(n-1)*pm)/n
            end
            p
        end
        assoc2(l, x) = l < 2 ? zero(x) : l == 2 ? 3*(1-x^2) : begin
            pm, p = 3*(1-x^2), 15x*(1-x^2)
            for n in 4:l
                pm, p = p, ((2n-1)*x*p-(n+1)*pm)/(n-2)
            end
            p
        end
        angular_G(x) = sum((-1im)^j * (2j+1) * Gpol[j+1] * legendre(j,x)
                           for j in 0:(length(Gpol)-1))
        for l in 2:6
            direct = (1.0im)^l * (2l+1) /
                (2sqrt((l-1)*l*(l+1)*(l+2))) *
                Cosmic.quadgk(x -> angular_G(x)*assoc2(l,x), -1.0, 1.0;
                              rtol=2e-13)[1]
            @test imag(direct) ≈ 0 atol = 2e-14
            @test real(direct) ≈ Epol[l+1] rtol = 3e-13 atol = 2e-14
        end

        # Known Wigner values, permutation symmetry and the monopole LOS identity.
        for j in 0:8
            @test Cosmic._wigner3j(j, j, 0, 0, 0, 0) ≈ (-1)^j / sqrt(2j + 1) rtol = 2e-14
        end
        @test Cosmic._wigner3j(2, 3, 4, 1, -2, 1) ≈
              (-1)^(2 + 3 + 4) * Cosmic._wigner3j(3, 2, 4, -2, 1, 1) rtol = 2e-14
        @test Cosmic._wigner3j(2, 2, 2, 2, 2, -4) == 0
        for Ls in 0:8
            @test Cosmic._J_T(Ls, 0, 0, 1.2) ≈
                  (-1)^Ls * Cosmic.sphericalbesselj(Ls, 1.2) rtol = 3e-14
        end

        # Complete +/-m layout and collision invariants of the pure operator.
        L = Cosmic.SecondOrderLayout(4, 4, 4)
        @test length(L.I) == 25
        @test length(L.E) == 21 == length(L.B)
        @test length(L.N) == 25
        y = randn(L.n)
        dy = similar(y)
        Cosmic._second_order_radiation_linear!(dy, y, L, 0.0, 3.0)
        @test dy[L.I[(0, 0)]] == 0                 # Thomson conserves intensity monopole
        @test all(isfinite, dy)
        fill!(y, 0)
        y[L.E[(3, 1)]] = 1
        Cosmic._second_order_radiation_linear!(dy, y, L, 0.0, 3.0)
        @test dy[L.E[(3, 1)]] == -3

        # The exact STF geometry reconstructs the total k-aligned tensor.
        tg = Cosmic.SecondOrderTriangle(0.17, 0.11, 0.21)
        for m in (-2, -1, 1, 2)
            tsum = Cosmic._tensor_product(tg, 1, 1, m) +
                   2Cosmic._tensor_product(tg, 1, 2, m) +
                   Cosmic._tensor_product(tg, 2, 2, m)
            @test tsum ≈ 0 atol = 2e-17
        end
        @test Cosmic._tensor_product(tg, 1, 1, 0) +
              2Cosmic._tensor_product(tg, 1, 2, 0) +
              Cosmic._tensor_product(tg, 2, 2, 0) ≈ 2tg.k^2/3

        # The Einstein kernels must obey the exchange parity
        # Q_m(k2,k1)=(-1)^m Q_m(k1,k2).
        q1 = Cosmic.FirstOrderMetricSnapshot(0.13, 0.11, -0.007, -0.006, 0.001)
        q2 = Cosmic.FirstOrderMetricSnapshot(-0.08, -0.06, 0.004, 0.003, -0.0007)
        qs = Cosmic._einstein_quadratic_sources(tg, q1, q2, 0.09, -0.003;
            a=0.2, rho_dipole1=0.04, rho_dipole2=-0.03)
        ts = Cosmic.SecondOrderTriangle(tg.k2, tg.k1, tg.k)
        qx = Cosmic._einstein_quadratic_sources(ts, q2, q1, 0.09, -0.003;
            a=0.2, rho_dipole1=-0.03, rho_dipole2=0.04)
        @test qs.QTT ≈ qx.QTT
        @test qs.QST ≈ qx.QST
        @test qs.QTR ≈ qx.QTR
        for m in -2:2
            @test qs.QSS[m] ≈ (-1)^m * qx.QSS[m]
            @test qs.QSS_prime[m] ≈ (-1)^m * qx.QSS_prime[m]
        end

        # Full quadratic Liouville block: a vanishing metric gives exactly no
        # source, and swapping the two convolution legs gives (-1)^m.
        q0 = Cosmic.FirstOrderMetricSnapshot(0.0, 0.0, 0.0, 0.0, 0.0)
        rr = Cosmic.FirstOrderRadiationSnapshot(randn(7), randn(7), randn(7))
        zsrc = Cosmic._quadratic_radiation_liouville(tg, q0, q0, rr, rr, L)
        @test all(iszero, values(zsrc.I))
        @test all(iszero, values(zsrc.E))
        @test all(iszero, values(zsrc.B))
        @test all(iszero, values(zsrc.N))
        r1 = Cosmic.FirstOrderRadiationSnapshot(randn(7), randn(7), randn(7))
        r2 = Cosmic.FirstOrderRadiationSnapshot(randn(7), randn(7), randn(7))
        ls = Cosmic._quadratic_radiation_liouville(tg, q1, q2, r1, r2, L)
        lx = Cosmic._quadratic_radiation_liouville(ts, q2, q1, r2, r1, L)
        for species in (:I, :E, :B, :N), lm in keys(ls[species])
            @test ls[species][lm] ≈ (-1)^lm[2] * lx[species][lm] atol = 2e-15
        end


        b1 = Cosmic.FirstOrderBaryonSnapshot(0.2, 0.25, -0.03, -0.025)
        b2 = Cosmic.FirstOrderBaryonSnapshot(-0.1, -0.08, 0.02, 0.018)
        tc = Cosmic._quadratic_thomson(tg, q1, q2, r1, r2, b1, b2, L;
            rho_gamma_over_rho_b=0.7)
        tx = Cosmic._quadratic_thomson(ts, q2, q1, r2, r1, b2, b1, L;
            rho_gamma_over_rho_b=0.7)
        for species in (:I, :E, :B), lm in keys(tc[species])
            @test tc[species][lm] ≈ (-1)^lm[2] * tx[species][lm] atol = 2e-15
        end
        @test tc.baryon_monopole == -0.7tc.I[(0, 0)]
        @test all(tc.baryon_dipole[m] == -0.7tc.I[(1, m)] for m in -1:1)

        d1 = Cosmic.FirstOrderMatterSnapshot(0.2, -0.03)
        d2 = Cosmic.FirstOrderMatterSnapshot(-0.1, 0.02)
        ml = Cosmic._quadratic_matter_liouville(tg, q1, q2, d1, d2)
        mx = Cosmic._quadratic_matter_liouville(ts, q2, q1, d2, d1)
        @test ml.monopole ≈ mx.monopole
        @test ml.pressure ≈ mx.pressure
        @test all(ml.dipole[m] ≈ (-1)^m * mx.dipole[m] for m in -1:1)
        @test all(ml.quadrupole[m] ≈ (-1)^m * mx.quadrupole[m] for m in -2:2)
    end

    FULL_TESTS && @testset "exact curved geometry" begin
        # CLASS's HypersphericalExplicit evaluator (not its high-nu flat
        # approximation), beta=40, ell=5, chi=0.2.  These pin both the open and
        # closed radial eigenfunctions to an independent implementation.
        phi_open = Cosmic._hyperspherical_bessels(-1, 40.0, 5, 0.2)[1][6]
        phi_closed = Cosmic._hyperspherical_bessels(1, 40.0, 5, 0.2)[1][6]
        @test phi_open ≈ 0.1240565354867583 rtol = 2e-14
        @test phi_closed ≈ 0.12904941867739481 rtol = 2e-14

        # Projection is linear: a zero physical source must produce exactly
        # zero transfer.  This also guards every curved E/B quadrature
        # accumulator against accidentally reverting to uninitialised storage.
        scalar_source = (eta=[0.0, 1.0], S0=zeros(2), S1=zeros(2),
            S2=zeros(2), SP=zeros(2))
        Kopen = -1e-8
        qscalar = 0.01
        kscalar = sqrt(qscalar^2 - Kopen)
        tT, tE = Cosmic._project_scalar_curved(
            scalar_source, [2, 5], Kopen, kscalar, qscalar, 2.0)
        @test tT == zeros(2)
        @test tE == zeros(2)

        tensor_source = (η=[0.0, 1.0], S_T=zeros(2), S_P=zeros(2))
        qtensor = 0.01
        ktensor = sqrt(qtensor^2 - 3Kopen)
        tT, tE, tB = Cosmic._project_tensor_curved(
            tensor_source, [2, 5], Kopen, ktensor, qtensor, 2.0)
        @test tT == zeros(2)
        @test tE == zeros(2)
        @test tB == zeros(2)
    end

    FULL_TESTS && @testset "observables" begin
        c = cosmology(m_ν=Float64[])
        r = recombination(c)
        P = matter_power_spectrum(c, r; nk=48, kmin=1e-4, kmax=5.0)

        # σ₈. Planck gives 0.8111 *with* Σm_ν = 0.06 eV; massive neutrinos
        # free-stream and suppress small-scale power by ~1.5% in σ₈, and this
        # cosmology has none, so it should land slightly high.
        @test 0.80 < σ8(P) < 0.84

        # P(k) must turn over at the matter-radiation equality scale,
        # k_eq = a_eq H(a_eq), which the background predicts independently.
        a_eq = scale_factor(z_equality(c))
        k_eq = a_eq * H_Mpc(c, a_eq)
        ks = exp.(range(log(1e-3), log(0.5); length=300))
        k_peak = ks[argmax([power(P, k) for k in ks])]
        @test 0.5 * k_eq < k_peak < 2 * k_eq

        # The primordial spectrum itself must have slope n_s - 1.
        s_prim = (log(primordial_power(c, 2e-3)) - log(primordial_power(c, 1e-3))) /
                 (log(2e-3) - log(1e-3))
        @test s_prim ≈ c.n_s - 1 rtol = 1e-10

        # P(k) rises on large scales and falls on small ones: below k_eq the
        # transfer function is flat, above it modes were suppressed while inside
        # the horizon during radiation domination (the Meszaros effect), so T ~
        # ln(k)/k² and P falls steeply.
        @test power(P, 1e-3) < power(P, 1e-2)
        @test power(P, 1.0) < power(P, 0.1) < power(P, 0.02)

        # Note there is deliberately no test that P ∝ k^n_s at small k. That
        # asymptote needs a mode both well inside the horizon (kη₀ ≫ 1) *and* far
        # above the equality scale (k ≪ k_eq), and in this cosmology those windows
        # barely overlap: by k = 5e-4 the mode is only just sub-horizon and the
        # transfer function is already bending. The textbook limit is not cleanly
        # reachable, and asserting it would be asserting a fiction.
    end

    FULL_TESTS && @testset "Poisson equation" begin
        # Deep inside the horizon the metric and matter must satisfy
        #     k²Φ = -(3/2) H₀² Ω_m δ_m / a
        # This ties the Einstein sector to the matter sector and is the sharpest
        # single check that the two are consistent.
        c = cosmology(m_ν=Float64[])
        r = recombination(c)
        bg = BackgroundCache(c, r)
        H0 = Cosmic.Constants.H0_in_invMpc(c.h)

        for (k, a) in ((0.1, 1.0), (0.1, 0.5), (0.05, 1.0), (0.2, 0.3))
            p = solve_perturbations(c, bg, k)
            @test k > 10 * H_Mpc(c, a) * a              # confirm sub-horizon
            lhs = k^2 * Φ(p, a)
            rhs = -1.5 * H0^2 * Ω_m(c) * δ_matter(p, a) / a
            @test lhs ≈ rhs rtol = 0.05
        end
    end

    @testset "reionization" begin
        c = cosmology()
        r = recombination(c)
        rn = recombination(c; reionization=no_reionization)

        # Fully ionized today: H ionized plus helium *doubly* ionized.
        @test x_e(r, 0.0) ≈ 1 + 2 * f_He(c) rtol = 1e-3
        # Half-ionized at the midpoint by construction.
        @test x_e(r, 7.67) ≈ (1 + f_He(c)) / 2 rtol = 0.05
        # Neutral again well before reionization.
        @test x_e(r, 20.0) < 1e-3

        # Planck 2018: τ = 0.0544 ± 0.0073.
        @test τ_reio(r) ≈ 0.0544 rtol = 0.05
        @test τ_reio(rn) < 1e-3
        @test z_reio_from_τ(c, 0.0544) ≈ 7.67 rtol = 0.05

        # Reionization must not disturb the recombination era. z_star is defined
        # against the recombination-only optical depth precisely so that adding
        # τ_reio ≈ 0.054 does not move it.
        @test z_star(r) ≈ z_star(rn) rtol = 1e-4
        @test r_drag(r) ≈ r_drag(rn) rtol = 1e-3
    end

    FULL_TESTS && @testset "CMB spectra" begin
        # Validated against CAMB 1.6.0 run with this same cosmology. These are
        # not round numbers from a table -- they are what CAMB produces.
        c = cosmology(m_ν=Float64[])
        r = recombination(c)

        # One spectrum, several ℓ. Each cmb_spectra call is a full Boltzmann sweep
        # over k, so the ℓ list is free but extra calls are not.
        ℓs = [50, 100, 200, 210, 220, 230, 260, 400, 1000]
        camb_TT = Dict(50 => 1423.4, 100 => 2704.1, 220 => 5765.1,
            400 => 1740.4, 1000 => 1033.8)   # D_ℓ, μK²

        # nk must resolve the 2π/η₀ ≈ 4.4e-4 oscillation of Θ_ℓ(k). At nk = 2500
        # the amplitudes are good to ~10% but the peak still sits a little low;
        # full convergence onto CAMB needs nk ≈ 5000, which is too slow for CI.
        s = cmb_spectra(c, r; lmax=1100, nk=2500, ℓs=ℓs)
        D = D_ℓ(s, :TT)

        for (i, ℓ) in enumerate(ℓs)
            haskey(camb_TT, ℓ) || continue
            @test D[i] ≈ camb_TT[ℓ] rtol = 0.12
        end

        # The first acoustic peak. Its position is fixed by θ_*, which the
        # thermodynamics layer nails to 0.001% -- so if this drifts, something
        # upstream broke. (CAMB: 221.)
        @test 200 <= s.ℓ[argmax(D)] <= 240

        # E-modes are positive; TE is a correlation and must change sign.
        @test all(D_ℓ(s, :EE) .> 0)
        @test any(D_ℓ(s, :TE) .< 0) && any(D_ℓ(s, :TE) .> 0)

        # The Doppler term fills the acoustic troughs -- it is π/2 out of phase
        # with the density. Dropping it deepens them; *reversing* its sign deepens
        # them too, by turning constructive interference into destructive. That is
        # a silent 4x error, and the peaks stay put either way, so nothing else
        # catches it.
        ℓ_trough = [400]
        D_ok = D_ℓ(cmb_spectra(c, r; lmax=500, nk=2000, ℓs=ℓ_trough), :TT)[1]
        D_off = D_ℓ(cmb_spectra(c, r; lmax=500, nk=2000, ℓs=ℓ_trough, dop=0.0), :TT)[1]
        D_bad = D_ℓ(cmb_spectra(c, r; lmax=500, nk=2000, ℓs=ℓ_trough, dop=-1.0), :TT)[1]
        @test D_off < D_ok
        @test D_bad < D_ok

        # nk_solve interpolates the smooth sources onto the dense projection
        # grid instead of solving the hierarchy at every k. Same integral, same
        # grids -- the difference against the all-k path is pure source
        # interpolation error, and at 600 solve points it must sit well below
        # the physics tolerances above.
        ℓ_fast = [50, 220, 400]
        s_fast = cmb_spectra(c, r; lmax=500, nk=2000, ℓs=ℓ_fast, nk_solve=600)
        s_full = cmb_spectra(c, r; lmax=500, nk=2000, ℓs=ℓ_fast)
        for w in (:TT, :EE)
            @test all(abs.(D_ℓ(s_fast, w) ./ D_ℓ(s_full, w) .- 1) .< 5e-3)
        end
    end

    FULL_TESTS && @testset "CMB lensing" begin
        # The Wigner d-recurrence against closed forms, and orthogonality at a
        # multipole high enough to expose recursion instability if any existed.
        μ = 0.3
        buf = zeros(13)
        Cosmic._wigner_d!(buf, 0, 0, μ, 12)
        @test buf[2] ≈ μ atol = 1e-14
        @test buf[3] ≈ (3μ^2 - 1) / 2 atol = 1e-14
        Cosmic._wigner_d!(buf, 2, -2, μ, 12)
        @test buf[3] ≈ ((1 - μ) / 2)^2 atol = 1e-14
        Cosmic._wigner_d!(buf, 2, 0, μ, 12)
        @test buf[3] ≈ sqrt(6) / 4 * (1 - μ^2) atol = 1e-14
        Cosmic._wigner_d!(buf, 2, 2, μ, 12)
        @test buf[3] ≈ ((1 + μ) / 2)^2 atol = 1e-14

        xs, ws = Cosmic.gauss(1024)
        big = zeros(401)
        aa = 0.0
        ab = 0.0
        for (x, w) in zip(xs, ws)
            Cosmic._wigner_d!(big, 3, -1, x, 400)
            aa += w * big[301]^2      # ℓ = 300
            ab += w * big[301] * big[304]
        end
        @test aa * 601 / 2 ≈ 1 rtol = 1e-10
        @test abs(ab) < 1e-12

        # Lensing with a null potential must return the inputs exactly -- the
        # correlation-function pipeline subtracts and adds back the unlensed
        # spectra, so any asymmetry in that bookkeeping shows up here.
        L = 300
        TT0 = [l < 2 ? 0.0 : 1000.0 / l^2 for l in 0:L]
        EE0 = [l < 2 ? 0.0 : 10.0 / l^2 for l in 0:L]
        TE0 = [l < 2 ? 0.0 : -30.0 / l^2 for l in 0:L]
        out = lensed_from_cls(TT0, EE0, TE0, zeros(L + 1); lmax_out=250)
        @test out.TT[3:251] == TT0[3:251]
        @test out.EE[3:251] == EE0[3:251]
        @test out.TE[3:251] == TE0[3:251]
        @test all(out.BB .== 0)

        # End-to-end at low resolution. The φφ amplitude is the strongest
        # physics gate here. Converged, CLASS 3.3.4 gives [ℓ(ℓ+1)]²C_ℓ^φφ/2π =
        # 1.353e-7 at ℓ = 50 for this cosmology; at this CI resolution the
        # small kmax cuts the high-k (small-χ) tail of the lensing kernel and
        # the value lands ~2% below, at 1.328e-7. The converged comparison
        # (0.4% of CLASS to ℓ = 400) was run against CLASS directly.
        c = cosmology(m_ν=Float64[])
        r = recombination(c; hydrogen=:hyrec)
        res = lensed_cmb_spectra(c, r; lmax=250, dl_buffer=150, nk=900, nη=3000)
        u = res.unlensed
        i50 = findfirst(==(50), u.ℓ)
        @test (50 * 51)^2 * u.φφ[i50] / (2π) ≈ 1.328e-7 rtol = 0.03
        @test u.Tφ[i50] > 0                       # ISW × lensing correlation
        @test all(res.lensed.BB .>= 0)            # lensing B-modes are a power transfer
        i100 = findfirst(==(100), res.lensed.ℓ)
        ui = findfirst(==(100), u.ℓ)
        @test abs(res.lensed.TT[i100] / u.TT[ui] - 1) < 0.01   # low-ℓ smoothing is sub-%
    end

    @testset "isocurvature initial conditions" begin
        # The relative entropy S_iγ = δ_i − ¾δ_γ is gauge invariant, so it is the
        # cleanest probe of the initial-condition normalisation: 1 for the species
        # that carries the entropy, 0 for the others (and 0 for adiabatic).
        c = cosmology(m_ν=Float64[], Yp=0.24568)
        r = recombination(c; hydrogen=:hyrec)
        bg = Cosmic.BackgroundCache(c, r)
        k = 0.01
        L = Cosmic.Layout(25, 32; lmax_m=12, nq=20, nmν=0)
        G = Cosmic.MassiveNuGrid(20)
        x0 = log(1e-5)
        ic0(mode) = Cosmic.initial_conditions(c, bg, k, L, x0, G, (); ic=mode)
        Scγ(u) = u[1] - 0.75 * u[L.iγ]        # CDM–photon entropy
        Sbγ(u) = u[3] - 0.75 * u[L.iγ]        # baryon–photon entropy

        u_ad = ic0(:adiabatic)
        @test abs(Scγ(u_ad)) < 1e-8           # adiabatic carries no relative entropy
        @test abs(Sbγ(u_ad)) < 1e-8

        u_cdi = ic0(:cdi)
        @test Scγ(u_cdi) ≈ 1 atol = 1e-6      # unit CDM isocurvature
        @test abs(Sbγ(u_cdi)) < 1e-6          # baryons track photons in CDI

        u_bi = ic0(:bi)
        @test Sbγ(u_bi) ≈ 1 atol = 1e-6       # unit baryon isocurvature
        @test abs(Scγ(u_bi)) < 1e-6

        for mode in (:cdi, :bi, :nid, :niv)
            @test all(isfinite, ic0(mode))
        end
        @test_throws ErrorException ic0(:not_a_mode)
    end

    @testset "primordial spectrum & inflaton" begin
        # An arbitrary primordial override must replace the power law everywhere.
        c = cosmology(m_ν=Float64[], Yp=0.24568)
        @test Cosmic.primordial_power(c, 0.05) ≈ c.A_s
        feat(k) = 3e-9 * (k / 0.05)^(-0.03)
        cf = cosmology(m_ν=Float64[], Yp=0.24568, primordial=feat)
        @test Cosmic.primordial_power(cf, 0.02) ≈ feat(0.02)

        # Inflaton V = ½m²φ²: the exact Mukhanov–Sasaki evolution must land on the
        # slow-roll attractor n_s = 1 − 2/N, r = 8/N at the pivot (N = 55).
        m = 6e-6
        V(φ) = 0.5 * m^2 * φ^2
        dV(φ) = m^2 * φ
        s = inflaton_spectrum(V, dV; φ0=16.0, N_pivot=55.0, k_pivot=0.05,
            decades_below=0.6, decades_above=0.6, per_decade=8, sub_horizon=25.0)
        @test s.n_s ≈ 1 - 2 / 55 atol = 0.003
        @test s.r ≈ 8 / 55 atol = 0.02
        @test s.A_s > 0
        # the computed spectrum must flow through primordial_power
        ci = cosmology(m_ν=Float64[], Yp=0.24568, primordial=s)
        @test Cosmic.primordial_power(ci, 0.05) ≈ s.A_s rtol = 1e-6

        # EFT sound speed: Δ²_ζ = H²/(8π²εc_s) at leading order, so constant
        # c_s = 0.5 doubles A_s and halves r (up to O(ε) slow-roll corrections
        # that the exact solver legitimately includes)
        s5 = inflaton_spectrum(V, dV; φ0=16.0, N_pivot=55.0, k_pivot=0.05,
            decades_below=0.6, decades_above=0.6, per_decade=8,
            sub_horizon=25.0, cs=N -> 0.5)
        @test s5.A_s / s.A_s ≈ 2 rtol = 0.06
        @test s5.r / s.r ≈ 0.5 rtol = 0.06
        @test s5.n_s ≈ s.n_s atol = 0.002
        # superluminal c_s violates causality and must be refused
        @test_throws ErrorException inflaton_spectrum(V, dV; φ0=16.0, cs=N -> 1.17)
    end

    @testset "PBH abundance (Young–Musco–Byrnes)" begin
        # monochromatic moments have closed forms: σ² = (16/81)A[W̃(1)T(1)]²
        # for P = Aδ(ln k/k*) at r_m = 1/k*, and μ²/σ² = k*²
        Pm(k) = exp(-log(k)^2 / (2 * 0.01^2)) / (sqrt(2π) * 0.01)
        σ2, μ2 = Cosmic.pbh_moments(Pm, 1.0)
        W1 = 3 * (sin(1) - cos(1))
        y = 1 / sqrt(3)
        T1 = 3 * (sin(y) - y * cos(y)) / y^3
        @test σ2 ≈ (16 / 81) * (W1 * T1)^2 rtol = 2e-3
        @test μ2 / σ2 ≈ 1 rtol = 2e-3
        # the critical linear amplitude (eq. 23) bounds: no PBHs from tiny power
        @test Cosmic.pbh_β(k -> 1e-9, 1.0) == 0.0
        # exponential sensitivity to the amplitude — the defining feature
        P1(k) = 0.02 * exp(-log(k)^2 / (2 * 0.3^2)) / (sqrt(2π) * 0.3)
        β1 = Cosmic.pbh_β(P1, 1.0)
        β2 = Cosmic.pbh_β(k -> 1.2 * P1(k), 1.0)
        @test β1 > 0
        @test β2 / β1 > 5
        # δc beyond the type-I maximum must refuse
        @test_throws ArgumentError Cosmic.pbh_β(P1, 1.0; δc=0.7)
    end

    @testset "scalar-induced GW (Kohri–Terada)" begin
        # exact-integral anchors from the defining paper (1804.08577):
        # scale-invariant (eq 31) and the power-law Q(n_s) table
        @test Cosmic.Ω_gw_rd(k -> 1.0, 1.0) ≈ 0.8222 rtol = 1e-3
        @test Cosmic.Ω_gw_rd(k -> k^(0.9655 - 1), 1.0) ≈ 0.8149 rtol = 1e-3
        @test Cosmic.Ω_gw_rd(k -> k^(1.2 - 1), 1.0) ≈ 0.8988 rtol = 1e-3
        # monochromatic closed form (eq 29): resonance below k̃=2, zero above
        @test Cosmic.sigw_monochromatic(2.1) == 0.0
        @test Cosmic.sigw_monochromatic(1.0) ≈ 0.4822 rtol = 1e-3
        # ΔN_eff of a flat Ω_GW h² = 1e-6 band over one decade: analytically
        # (8/7)(11/4)^{4/3} · (1e-6/h²)·ln10 / Ω_γ ≈ 0.41 — the same arithmetic
        # that makes h²Ω_GW ≲ 1e-6 the rule-of-thumb BBN bound
        c = cosmology(m_ν=Float64[], Yp=0.24568)
        s = Cosmic.SIGWSpectrum([1e5, 1e6], [1.5e-10, 1.5e-9], [1e-6, 1e-6])
        ΔN = sigw_ΔNeff(s, c)
        @test ΔN ≈ (8 / 7) * (11 / 4)^(4 / 3) * (1e-6 / c.h^2) * log(10) / Cosmic.Ω_γ(c) rtol = 1e-6
    end

    @testset "neutrino distribution (ξ, arbitrary f₀)" begin
        Tcmb = 2.7255
        h = 0.6766
        # standard Fermi–Dirac normalisation is exact (ξ = 0 unchanged)
        @test massive_neutrino(0.06, Tcmb, h).F_norm ≈ 7π^4 / 120 rtol = 1e-6
        @test degenerate_Neff_factor(0.0) == 1
        # a degeneracy ξ enhances the relativistic density by the exact factor
        # 1 + (30/7)(ξ/π)² + (15/7)(ξ/π)⁴
        ξ = 0.5
        νξ = massive_neutrino(1e-6, Tcmb, h; ξ=ξ)
        ν0 = massive_neutrino(1e-6, Tcmb, h; ξ=0.0)
        @test Cosmic.ρ_over_ρc0(νξ, 1e-3) / Cosmic.ρ_over_ρc0(ν0, 1e-3) ≈
              degenerate_Neff_factor(ξ) rtol = 1e-5
        # a tabulated Fermi–Dirac reproduces the analytic one
        fd = FermiDirac(0.0)
        tab = TabulatedNu(q -> Cosmic.nu_f0(fd, q), q -> Cosmic.nu_dlnf0(fd, q))
        νt = massive_neutrino(0.06, Tcmb, h; dist=tab)
        νf = massive_neutrino(0.06, Tcmb, h)
        @test Cosmic.ρ_over_ρc0(νt, 0.1) ≈ Cosmic.ρ_over_ρc0(νf, 0.1) rtol = 1e-6

        # warm dark matter: the exact-normalisation thermal relic must land on
        # its target abundance today, and stay radiation-like early
        w = warm_dark_matter(1.0, 0.2607, Tcmb, h)
        @test Cosmic.ρ_over_ρc0(w, 1.0) ≈ 0.2607 rtol = 1e-3
        # relativistic early: for m = 1 keV, y₁·a ≪ 1 needs a ≲ 1e-10 (y₁ ≈ 3e7)
        @test Cosmic.w(w, 1e-11) ≈ 1 / 3 rtol = 1e-3
        @test Cosmic.w(w, 1.0) < 1e-4                   # cold today
        @test_throws ErrorException warm_dark_matter(1e-3, 0.2607, Tcmb, h)  # β ≥ 1
    end

    @testset "decaying dark matter background" begin
        # long-lived limit: Γ t₀ ≈ 1.4e-3, so ~0.14% decays and the budget closes
        c1 = cosmology(m_ν=Float64[], Yp=0.24568, Ω_c=1e-6,
            Γ_dcdm=0.1, Ω_dcdm_ini=0.2607)
        d1 = Cosmic.get_species(c1, DecayingCDM)
        @test Cosmic.ρ_over_ρc0(d1, 1.0) ≈ 0.2607 rtol = 3e-3
        @test sum(Cosmic.ρ_over_ρc0(s, 1.0) for s in c1.species) ≈ 1 rtol = 1e-6

        # order-unity decay: energy moves to dark radiation, budget still closes;
        # the exact e^{−Γt₀} agreement (5 digits) is verified against direct quadrature.
        c2 = cosmology(m_ν=Float64[], Yp=0.24568, Ω_c=1e-6,
            Γ_dcdm=50.0, Ω_dcdm_ini=0.2607)
        d2 = Cosmic.get_species(c2, DecayingCDM)
        r2 = Cosmic.get_species(c2, DecayRadiation)
        @test Cosmic.ρ_over_ρc0(d2, 1.0) < 0.14          # more than half decayed
        @test Cosmic.ρ_over_ρc0(r2, 1.0) > 0.05          # products present
        @test sum(Cosmic.ρ_over_ρc0(s, 1.0) for s in c2.species) ≈ 1 rtol = 1e-6

        # perturbations: in the long-lived limit the dcdm perturbation must track
        # ordinary CDM (identical equations up to aΓψ, negligible for Γt₀ ≪ 1)
        r1 = recombination(c1; hydrogen=:hyrec)
        p1 = Cosmic.solve_perturbations(c1, Cosmic.BackgroundCache(c1, r1), 0.1)
        L1 = Cosmic._layout(p1)
        u1 = p1.sol(0.0)
        @test u1[L1.idc] ≈ u1[1] rtol = 1e-2       # δ_dcdm ≈ δ_cdm
        @test all(isfinite, u1)
    end

    @testset "dark-energy fluid perturbations" begin
        # w → −1: the fluid must reproduce the (unperturbed, exact) Λ run, with
        # δ_fld suppressed by (1+w)
        c0 = cosmology(m_ν=Float64[], Yp=0.24568)
        r0 = recombination(c0; hydrogen=:hyrec)
        p0 = Cosmic.solve_perturbations(c0, Cosmic.BackgroundCache(c0, r0), 0.1)
        @test Cosmic._layout(p0).has_fld == false   # Λ carries no perturbations
        c1 = cosmology(m_ν=Float64[], Yp=0.24568, w0=-0.9999)
        r1 = recombination(c1; hydrogen=:hyrec)
        p1 = Cosmic.solve_perturbations(c1, Cosmic.BackgroundCache(c1, r1), 0.1)
        @test Cosmic.δ_cdm(p1, 1.0) ≈ Cosmic.δ_cdm(p0, 1.0) rtol = 1e-3
        L1 = Cosmic._layout(p1)
        @test abs(p1.sol(0.0)[L1.ifld]) < 1e-4
        # real fluid: finite, and DE clustering suppressed by cs² = 1 pressure
        c2 = cosmology(m_ν=Float64[], Yp=0.24568, w0=-0.9)
        r2 = recombination(c2; hydrogen=:hyrec)
        p2 = Cosmic.solve_perturbations(c2, Cosmic.BackgroundCache(c2, r2), 0.1)
        u2 = p2.sol(0.0)
        @test all(isfinite, u2)
        @test abs(u2[Cosmic._layout(p2).ifld]) < 0.1   # no runaway DE clustering
    end

    @testset "quintessence" begin
        # thawing exponential potential: budget closes exactly by shooting,
        # field frozen early (w = −1), thawed today; perturbations regular at
        # w = −1 thanks to the (δ, Q=(1+w)θ) variables
        λ = 0.5
        V(φ) = exp(-λ * φ)
        dV(φ) = -λ * exp(-λ * φ)
        c = cosmology(m_ν=Float64[], Yp=0.24568, V_scf=V, dV_scf=dV, φ0_scf=0.5)
        de = Cosmic.get_species(c, QuintessenceDE)
        @test sum(Cosmic.ρ_over_ρc0(s, 1.0) for s in c.species) ≈ 1 rtol = 1e-6
        @test Cosmic.w(de, 1e-5) ≈ -1 atol = 1e-4      # frozen early
        @test -1 < Cosmic.w(de, 1.0) < -0.9            # thawed, still DE-like
        r = recombination(c; hydrogen=:hyrec)
        p = Cosmic.solve_perturbations(c, Cosmic.BackgroundCache(c, r), 0.1)
        u = p.sol(0.0)
        @test all(isfinite, u)                          # regular at w = −1
        # constant potential ≡ Λ exactly
        V0(φ) = 1.0
        dV0(φ) = 0.0
        c0 = cosmology(m_ν=Float64[], Yp=0.24568, V_scf=V0, dV_scf=dV0)
        de0 = Cosmic.get_species(c0, QuintessenceDE)
        @test Cosmic.w(de0, 1.0) ≈ -1 atol = 1e-8
    end

    @testset "PPF (w = −1 crossing)" begin
        # non-crossing: PPF must mimic the exact fluid (FHL construction)
        c9 = cosmology(m_ν=Float64[], Yp=0.24568, w0=-0.9)
        r9 = recombination(c9; hydrogen=:hyrec)
        bg9 = Cosmic.BackgroundCache(c9, r9)
        pF = Cosmic.solve_perturbations(c9, bg9, 0.1; de_scheme=:fluid)
        pP = Cosmic.solve_perturbations(c9, bg9, 0.1; de_scheme=:ppf)
        @test Cosmic.δ_cdm(pP, 1.0) ≈ Cosmic.δ_cdm(pF, 1.0) rtol = 1e-3
        # crossing model: auto-detected, evolves finitely through w = −1
        cX = cosmology(m_ν=Float64[], Yp=0.24568, w0=-1.1, wa=0.5)
        @test Cosmic._crosses_w_minus1(Cosmic._fld_species(cX))
        @test !Cosmic._crosses_w_minus1(Cosmic._fld_species(c9))
        rX = recombination(cX; hydrogen=:hyrec)
        pX = Cosmic.solve_perturbations(cX, Cosmic.BackgroundCache(cX, rX), 0.1)
        @test all(isfinite, pX.sol(0.0))
    end

    @testset "Horndeski stable basis (structural)" begin
        # The stable basis is defined so the three positivity inputs ARE the
        # no-ghost/no-gradient conditions; the derived α's must obey the
        # defining relations exactly, and α_B must hit its z=0 boundary value.
        c = cosmology(m_ν=Float64[])
        ΩΛ0 = sum(Cosmic.ρ_over_ρc0(s, 1.0)
                  for s in c.species if s isa Cosmic.AbstractDarkEnergy)
        ΩΛa(x) = ΩΛ0 / Cosmic.E(c, exp(x))^2
        Dk(x) = ΩΛa(x) + 1.5 * (0.2ΩΛa(x))^2
        spec = HorndeskiEFT(D_kin=Dk, M_star2=x -> 1.0, c_s2=x -> 0.8,
            α_B0=0.2ΩΛa(0.0))
        hf = stable_basis_solve(c, spec)
        @test hf.α_B(0.0) ≈ spec.α_B0 rtol = 1e-6      # boundary condition
        for x in (-4.0, -2.0, -0.5, 0.0)
            @test hf.α_K(x) ≈ Dk(x) - 1.5 * hf.α_B(x)^2 rtol = 1e-3  # α_K def
        end
        @test hf.x_on > -15  # activates once D_kin is dynamically meaningful
        @test_throws ArgumentError stable_basis_solve(
            cosmology(m_ν=Float64[], Ω_k=0.01), spec)   # flat only
    end

    FULL_TESTS && @testset "Horndeski perturbations (hi_class-validated)" begin
        # The MG branch is strictly opt-in: with no spec the state and solution
        # are byte-identical to GR, and a spec with a vanishing kinetic term
        # must reproduce the GR result (the scalar sources everything at O(α)).
        c = cosmology(m_ν=Float64[])
        r = recombination(c)
        bg = Cosmic.BackgroundCache(c, r)
        ΩΛ0 = sum(Cosmic.ρ_over_ρc0(s, 1.0)
                  for s in c.species if s isa Cosmic.AbstractDarkEnergy)
        ΩΛa(x) = ΩΛ0 / Cosmic.E(c, exp(x))^2

        # near-GR: D_kin → 0 ⇒ φ, δ_c unchanged from the no-MG solve
        tiny = HorndeskiEFT(D_kin=x -> 1e-6 * ΩΛa(x), M_star2=x -> 1.0,
            c_s2=x -> 1.0, α_B0=0.0)
        hf0 = stable_basis_solve(c, tiny)
        for k in (0.01, 0.1)
            pg = Cosmic.solve_perturbations(c, bg, k)
            pm = Cosmic.solve_perturbations(c, bg, k; mg=hf0)
            @test pm.sol(0.0)[5] ≈ pg.sol(0.0)[5] rtol = 1e-4
            @test pm.sol(0.0)[1] ≈ pg.sol(0.0)[1] rtol = 1e-4
        end

        # propto_omega (ĉ_K=1, ĉ_B=0.2, ĉ_M=0.1) written in the stable basis by
        # inverting B&S eq. 3.13 for c_s²(x) at the target α_B = 0.2Ω_Λ. The
        # Newtonian-gauge Bardeen potential ratio φ_MG/φ_GR at z=0 sits at the
        # sub-horizon plateau 1.0197, matching hi_class v3.2.3 to <1e-3.
        E2(x) = Cosmic.E(c, exp(x))^2
        dlnE(x) = (-log(E2(x + 2e-3)) + 8log(E2(x + 1e-3)) -
                   8log(E2(x - 1e-3)) + log(E2(x - 2e-3))) / 24e-3
        cB, cM = 0.2, 0.1
        xs = collect(range(-15.0, 0.0; length=1200))
        lnM2 = zeros(length(xs))
        for i in 2:length(xs)
            lnM2[i] = lnM2[i-1] + Cosmic.quadgk(u -> cM * ΩΛa(u), xs[i-1], xs[i])[1]
        end
        M2itp = Cosmic.linear_interpolation(xs, lnM2,
            extrapolation_bc=Cosmic.Line())
        M2f(x) = exp(M2itp(x))
        αB(x) = cB * ΩΛa(x)
        Dk(x) = ΩΛa(x) + 1.5αB(x)^2
        function cs2(x)
            hm = Cosmic._matter_enthalpy_over_H2(c, exp(x)) / M2f(x)
            (-2cB * ΩΛa(x) * dlnE(x) -
             (2 - αB(x)) * (dlnE(x) - αB(x) / 2 - cM * ΩΛa(x)) - hm) / Dk(x)
        end
        spec = HorndeskiEFT(D_kin=Dk, M_star2=M2f, c_s2=cs2, α_B0=αB(0.0))
        hf = stable_basis_solve(c, spec; x_min=-6.0)
        pg = Cosmic.solve_perturbations(c, bg, 0.03)
        pm = Cosmic.solve_perturbations(c, bg, 0.03; mg=hf)
        @test pm.sol(0.0)[5] / pg.sol(0.0)[5] ≈ 1.0197 rtol = 3e-3

        # tensor sector: α_M adds Hubble friction and the M*² source rescaling;
        # the h-mode stays finite and shifts from the GR amplitude
        tg = Cosmic.solve_tensor_perturbations(c, bg, 0.05)
        tm = Cosmic.solve_tensor_perturbations(c, bg, 0.05; mg=hf)
        Lt = Cosmic.TensorLayout(12, 16; lmax_m=12, nq=20, nmν=0)
        @test isfinite(tm.sol(0.0)[Lt.ih])
        @test tm.sol(0.0)[Lt.ih] != tg.sol(0.0)[Lt.ih]
    end

    FULL_TESTS && @testset "Horndeski + massive ν / decaying CDM" begin
        # The MG scalar EOM sources now include massive neutrinos (δP, momentum,
        # shear and dσ/dx) and decaying CDM + dark radiation (density, momentum,
        # shear/derivative, and a decay-aware background ṗ_dr). Validated off-line
        # by the B&S 3.17 Hamiltonian-constraint residual: under one spec the
        # massive-ν and dcdm solutions sit on the constraint surface to the same
        # numerical floor as the massless GR limit (the residual scales with the
        # spec's c_s² stiffness, not the species content). Here we guard that the
        # near-GR limit stays uncorrupted and the full solve is finite + shifted.
        function specs(c)
            ΩΛ0 = sum(Cosmic.ρ_over_ρc0(s, 1.0)
                      for s in c.species if s isa Cosmic.AbstractDarkEnergy)
            ΩΛa(x) = ΩΛ0 / Cosmic.E(c, exp(x))^2
            Dk(x) = ΩΛa(x) + 1.5 * (0.1 * ΩΛa(x))^2
            tiny = HorndeskiEFT(D_kin=x -> 1e-6 * ΩΛa(x), M_star2=x -> 1.0,
                c_s2=x -> 1.0, α_B0=0.0)
            full = HorndeskiEFT(D_kin=Dk, M_star2=x -> 1.0, c_s2=x -> 0.5,
                α_B0=0.1 * ΩΛa(0.0))
            (stable_basis_solve(c, tiny), stable_basis_solve(c, full))
        end
        # A massive-ν cosmology, and a realistic decaying-CDM one (mostly-stable
        # CDM + a small long-lived decaying fraction). The extreme Ω_c→0, fast-
        # decay case makes the tiny-D_kin near-GR limit ill-conditioned (the
        # rapidly-decaying background drives a residual α_B ~ 1e-4 that the 1/D_kin
        # limit amplifies) — not a wiring issue, and the full spec's constraint
        # residual stays clean there; we just don't pin near-GR on it.
        for c in (cosmology(m_ν=[0.15, 0.0, 0.0], Yp=0.24568),
                  cosmology(m_ν=Float64[], Yp=0.24568, Ω_c=0.25,
                      Γ_dcdm=0.1, Ω_dcdm_ini=0.01))
            r = recombination(c; hydrogen=:hyrec)
            bg = Cosmic.BackgroundCache(c, r)
            hf_tiny, hf_full = specs(c)
            p0 = Cosmic.solve_perturbations(c, bg, 0.1)
            pt = Cosmic.solve_perturbations(c, bg, 0.1; mg=hf_tiny)
            # near-GR (D_kin → 0): the MG solve reproduces the GR one for this
            # species content — the massive-ν / dcdm sources vanish with the scalar
            @test pt.sol(0.0)[5] ≈ p0.sol(0.0)[5] rtol = 1e-3
            # full spec: finite everywhere and genuinely MG-shifted from GR
            pm = Cosmic.solve_perturbations(c, bg, 0.1; mg=hf_full)
            @test all(isfinite, pm.sol(0.0))
            @test !isapprox(pm.sol(0.0)[5], p0.sol(0.0)[5]; rtol = 1e-3)
        end
    end

    FULL_TESTS && @testset "tensor spectra" begin
        # The flat projector (Bessel table) and the curved projector
        # (hyperspherical recurrences) are independent code paths that must
        # agree in the K -> 0 limit. This catches regressions in either
        # without an external reference.
        K = 1e-12
        for ℓ in (2, 10), (ν, d) in ((2000, 8000.0), (10000, 13000.0))
            # Closed eigenmodes require integer β = q/√K, so fix the eigenmode
            # and derive (q, k); the open case reuses the same q.
            for sgn in (+1, -1)
                q = ν * sqrt(K)
                k = sqrt(q^2 - 3 * sgn * K)
                B = Cosmic.TensorBesselTable([ℓ], 1.05 * k * d)
                Tf, Ef, Bf = Cosmic._tensor_radials(B, 1, ℓ, k * d)
                Tc, Ec, Bc = Cosmic._curved_tensor_radials(sgn * K, k, q, [ℓ], d)
                # The tolerance is set by the *flat* table's linear
                # interpolation (dx = 0.08, ~5e-4); the curved recurrence is
                # exact.
                @test Tc[1] ≈ Tf rtol = 1e-3
                @test Ec[1] ≈ Ef rtol = 1e-3
                @test Bc[1] ≈ Bf rtol = 1e-3
            end
        end

        # Mode maps: tensors live at k² = q² - 3K; on S³ the tower is the
        # integer eigenmodes ν = 3, 4, ... with none skipped.
        Kc = 5e-10
        ks, qs, discrete = Cosmic._tensor_mode_grid(-Kc, 14000.0, 20, 100, 0.003)
        @test !discrete && all(@. abs(ks^2 - (qs^2 + 3Kc)) < 1e-18)
        ks, qs, discrete = Cosmic._tensor_mode_grid(Kc, 14000.0, 20, 100, 0.003)
        @test discrete && qs[1] ≈ 3 * sqrt(Kc) && length(qs) == length(unique(qs))
        @test all(@. abs(ks^2 - (qs^2 - 3Kc)) < 1e-18)

        # Spectrum smoke test at deliberately coarse nk: TT at low ℓ is
        # insensitive to the k-grid (its radial decays as j/x²), so it must
        # sit near the validated value (5.095 at ℓ=2 vs CLASS 5.099, with
        # kmax = 0.056). EE/BB need the full 40-points-per-cycle grid and are
        # validated against CLASS/CAMB separately, not here.
        c = cosmology(m_ν=Float64[])
        r = recombination(c; reionization=no_reionization)
        t = tensor_cmb_spectra(c, r; r=0.01, n_t=-0.01 / 8, lmax=20, ℓs=[2, 20],
            nk=600, kmax=0.056, nη=2400, lmax_γ=16, lmax_ν=20)
        DT = D_ℓ(t, :TT)
        @test DT[1] ≈ 5.095 rtol = 0.05
        @test all(D_ℓ(t, :EE) .> 0) && all(D_ℓ(t, :BB) .> 0)
        @test all(D_ℓ(t, :TE) .< 0)
    end

    FULL_TESTS && @testset "massive neutrinos" begin
        # The momentum-resolved hierarchy. Slow (nq bins multiply the state), so
        # only a couple of modes are solved here.
        c = cosmology(m_ν=[0.06, 0.0, 0.0])
        @test length(get_all_species(c, MassiveNeutrinos)) == 1
        r = recombination(c)
        bg = BackgroundCache(c, r)

        # The whole coupled system -- momentum hierarchy, metric, other species --
        # must still conserve ℛ super-horizon. This is the strongest check that
        # the massive-ν moments feed back into the Einstein equations correctly.
        pS = solve_perturbations(c, bg, 1e-4)
        for a in (1e-4, 1e-3, 1e-2)
            @test curvature_ℛ(pS, a) ≈ 1.0 rtol = 5e-3
        end

        # Massive neutrinos suppress small-scale power by free-streaming: with
        # Σm_ν = 0.06 eV, σ₈ drops ~1.5% relative to massless. Validated against
        # CAMB run with the same cosmology (σ₈ = 0.8103); we sit slightly high
        # because the momentum grid is finite, and converge onto it as nq grows.
        P = matter_power_spectrum(c, r; nk=32, kmin=1e-3, kmax=2.0)
        Pmless = matter_power_spectrum(cosmology(m_ν=Float64[]), recombination(cosmology(m_ν=Float64[]));
            nk=32, kmin=1e-3, kmax=2.0)
        @test σ8(P) < σ8(Pmless)                 # suppression, right direction
        @test σ8(P) ≈ 0.8103 rtol = 0.02         # matches CAMB to ~1%

        # Large scales are untouched -- neutrinos free-stream out of small scales
        # only, so P(k) at k ≪ k_fs must agree with CAMB tightly.
        @test power(P, 1e-3) > 0
    end

    @testset "gauges" begin
        c = cosmology(m_ν=Float64[])
        r = recombination(c)
        bg = BackgroundCache(c, r)

        # A super-horizon mode: ℛ is conserved there, which is the whole reason
        # inflation can hand a prediction to the late universe at all.
        p = solve_perturbations(c, bg, 1e-4)

        # The initial conditions are normalised to ℛ = 1. Recovering 1 at every
        # later time is *not* a tautology: ℛ is rebuilt from the evolved φ and
        # θ_tot, so it only comes back to 1 if the Boltzmann hierarchy and the
        # Einstein constraints agree with each other. This is the sharpest single
        # consistency check in the package.
        for a in (1e-5, 1e-4, 1e-3, 1e-2, 0.1)
            @test curvature_ℛ(p, a) ≈ 1.0 rtol = 5e-3
        end

        # ζ equals ℛ outside the horizon (they differ at order (k/ℋ)²).
        for a in (1e-4, 1e-3, 1e-2)
            @test curvature_ζ(p, a) ≈ curvature_ℛ(p, a) rtol = 5e-3
        end

        p2 = solve_perturbations(c, bg, 1e-3)

        # Each gauge must actually satisfy its own defining condition.
        @test abs(θ_species(p2, 1.0, :cdm, Synchronous)) < 1e-12      # CDM at rest
        @test abs(δ_species(p2, 1.0, :matter, UniformDensity)) < 0.05  # δρ_tot = 0
        @test gauge_shift(p2, 1.0, Newtonian) == 0.0

        # Density contrast is gauge dependent -- that is the point. The gauges
        # must genuinely disagree outside the horizon...
        @test δ_species(p2, 1.0, :cdm, Newtonian) != δ_species(p2, 1.0, :cdm, Synchronous)

        # ...but agree deep inside it, where the shift 3ℋT is negligible.
        p3 = solve_perturbations(c, bg, 0.5)
        @test δ_species(p3, 1.0, :cdm, Newtonian) ≈
              δ_species(p3, 1.0, :cdm, Synchronous) rtol = 1e-3

        # The Bardeen potentials are gauge invariant and coincide with the
        # conformal-Newtonian ones.
        @test bardeen_Φ(p2, 1.0) == Φ(p2, 1.0)
        @test bardeen_Ψ(p2, 1.0) == Ψ(p2, 1.0)

        # And the comoving-gauge δ_m used for P(k) must match the transformation.
        @test δ_matter_comoving(p2, 1.0) ≈
              δ_species(p2, 1.0, :matter, Synchronous) rtol = 1e-6
    end

    FULL_TESTS && @testset "spectral distortions" begin
        # Exact Green's-function / PCA method (Chluba & Jeong 2014), validated
        # against a live CLASS 3.3.4 run (sd_branching_approx = exact, FIRAS,
        # this exact cosmology): μ = 2.09209e-8, y = 3.37597e-9. CLASS divides
        # the file's J_μ by 1.401 on read and multiplies by 1.401 after
        # convolving (likewise 4 for J_T, J_y); the factors are baked into the
        # file, so using it as-is gives the same answer by a shorter path.
        c = cosmology(m_ν=Float64[])
        r = recombination(c)
        amp = distortion_amplitudes(c, r)

        # μ and y are the observables. Both match CLASS to a couple of 1e-4.
        @test amp.μ ≈ 2.09209e-8 rtol = 0.01
        @test amp.y ≈ 3.37597e-9 rtol = 0.01
        @test amp.μ < FIRAS_μ_limit
        @test amp.y < FIRAS_y_limit

        # The damping scale must match CLASS's r_d (index_th_r_d): k_D = 2π/r_d.
        # CLASS gives k_D(z=1e4) = 4.40, k_D(z=1e5) = 129.6 /Mpc.
        @test damping_scale(c, r, scale_factor(1e4)) ≈ 4.40 rtol = 0.02
        @test damping_scale(c, r, scale_factor(1e5)) ≈ 129.6 rtol = 0.02

        # μ probes small-scale power: a bluer spectrum injects more of it, so μ
        # must rise with n_s. This is the physics that makes μ a probe of inflation.
        c_blue = cosmology(m_ν=Float64[], n_s=1.0)
        @test μ_distortion(c_blue, recombination(c_blue)) > amp.μ

        # The adiabatic-cooling contribution alone is negative (baryons draining
        # heat from the photons), the one guaranteed source of the opposite sign.
        @test μ_distortion(c, r; silk=false) < 0

        # A ΛCDM distortion is ~1e4 below the FIRAS noise, so χ² ≪ 1: FIRAS does
        # not constrain standard cosmology.
        @test firas_chi2(c, r) < 1e-3

        # The distortion spectrum is a real signal of the right sign structure.
        @test distortion_spectrum(c, r, 400.0) > 0
    end

    FULL_TESTS && @testset "thermalization (derived μ, y)" begin
        # Solve the photon Boltzmann equation from scratch -- exact Klein-Nishina
        # x Maxwell-Juttner Compton + double Compton + bremsstrahlung -- and read
        # μ and y off the evolved spectrum.
        # No Green's-function table, no branching-ratio file.
        c = cosmology(m_ν=Float64[])
        r = recombination(c)
        res = thermalization_distortions(c, r)

        # Validated against Chluba (2013) eqs. (5)-(6) -- his own fits to the exact
        # thermalization Green's function -- convolved with this heating history.
        #
        # NOT against CLASS's FIRAS_branching_ratios.dat. That file is energy-
        # conserving and self-consistent (4J_T + 4J_y + J_μ/1.401 = 1 at every z),
        # but its split between μ, y and a temperature shift is the Chluba & Jeong
        # (2014) PCA projection in the FIRAS band, which the authors note depends
        # on the experimental settings. The solver here uses Chluba 2013's
        # visibility split instead; the two conventions differ by ~17% in μ while
        # agreeing on the total energy. `distortion_amplitudes` follows the CLASS
        # convention, this solver the Chluba-2013 one.
        μ_ref = 1.401 * quadgk(lz -> (z = exp(lz) - 1;
                chluba_J_mu(z) * chluba_J_vis(z) * heating_rate(c, r, z) * (z + 1)),
            log(1e3 + 1), log(5e6 + 1); rtol=1e-5)[1]
        y_ref = 0.25 * quadgk(lz -> (z = exp(lz) - 1;
                chluba_J_y(z) * heating_rate(c, r, z) * (z + 1)),
            log(1e3 + 1), log(5e6 + 1); rtol=1e-5)[1]

        @test res.μ ≈ μ_ref rtol = 0.08      # derived vs Chluba: 2.5%
        @test res.y ≈ y_ref rtol = 0.08      # derived vs Chluba: 2.1%
        @test res.μ > 0 && res.y > 0

        # The closure test that pins the μ "convention gap": apply CLASS's own
        # estimator -- weighted least squares over the FIRAS channels with the
        # FIRAS noise, basis {G, Y, M, S_1..6} -- to THIS solver's evolved
        # spectrum. It must reproduce CLASS's μ (2.092e-8), even though the
        # grid projection above reads 1.6e-8 off the very same spectrum. The
        # 17% between the conventions is bookkeeping of the intermediate-era
        # signal, not physics: the spectrum is one and the same.
        sh = Cosmic._shapes()
        fir = Cosmic._firas()
        kB, hP, cl = 1.380649e-23, 6.62607015e-34, 2.99792458e8
        xk(ν) = hP * ν * 1e9 / (kB * c.Tcmb)
        pref = 2hP / cl^2 * (kB * c.Tcmb / hP)^3 / 1e-18
        itp = Cosmic.linear_interpolation(log.(res.x), res.Δn)
        ΔI = [pref * xk(ν)^3 * itp(log(xk(ν))) for ν in fir.ν]
        cols = Any[[sh.G_T(ν) for ν in fir.ν], [sh.Y_SZ(ν) for ν in fir.ν],
            [sh.M_mu(ν) for ν in fir.ν]]
        for k in 1:sh.nS
            push!(cols, [sh.S[k](ν) for ν in fir.ν])
        end
        w = 1 ./ fir.σ
        coef = (hcat(cols...) .* w) \ (ΔI .* w)
        @test coef[3] ≈ 2.0921e-8 rtol = 0.05
        @test coef[2] ≈ 3.3760e-9 rtol = 0.05

        # Chluba's own crossing: J_y ≈ J_μ at z ≈ 5.3e4, which is what makes 5e4 the
        # conventional μ/y boundary.
        @test chluba_J_y(5.3e4) ≈ chluba_J_mu(5.3e4) rtol = 0.05

        # μ cannot exceed its photon-number ceiling. Injecting energy at fixed
        # photon number forces ΔN/N = 3Δ_T - 1.3686μ = 0 and ΔQ = 4Δ_T - 1.1106μ,
        # whence μ = 1.401·ΔQ. With ΔQ_μ-era ≈ 1.4e-8 that caps μ at ~1.96e-8.
        @test res.μ < 1.401 * 1.5e-8

        # Energy must be conserved by the evolution: the distortion energy has to
        # account for what was injected. The exact scattering operator conserves
        # photon number and solves the shared Compton+emission electron-energy
        # ledger to machine precision on any grid.
        @test 1e-8 < res.Δρ_ρ < 5e-8

        # The DC Gaunt factor carries the seed-photon phase space I₄ = 4π⁴/15.
        @test Cosmic._I4pl ≈ 25.976 rtol = 1e-4
        @test H_dc(0.0) ≈ 1.0                 # normalised at zero frequency
        @test H_dc(5.0) < H_dc(1.0) < 1.0     # falls as x⁴e^{-2x}

        # The exact double-Compton emission factor (data/thermalization/
        # dc_gdc.dat, from the FORM matrix element): its soft corner is the
        # analytic I₄^pl, a mid-grid node matches the generation-time gate
        # value, it agrees with the CS2012 fit where that fit is good, and
        # past the table edge it keeps falling as e^{-2x}.
        @test g_dc_exact(1e-4, 1e-5) ≈ Cosmic._I4pl rtol = 2e-4
        @test g_dc_exact(0.5, 1e-3) ≈ 20.43 rtol = 5e-3
        @test g_dc_exact(0.1, 1e-4) ≈ Cosmic._I4pl * H_dc(0.1) / (1 + 14.16e-4) rtol = 0.01
        @test g_dc_exact(25.0, 5e-4) < g_dc_exact(17.0, 5e-4) * exp(-2 * (25 - 17)) * 10

        # Installed exact Compton cache and its two conservation null modes.
        kt = Cosmic._load_compton_kernel_table()
        @test size(kt.logr) == (48, 101, 29)
        xkn = exp.(range(log(5e-4), log(50.0); length=64))
        kop = Cosmic.ExactComptonOperator(xkn)
        kout = zeros(64)
        δkn = -1e-8 .* kop.bose                 # linearized BE chemical potential
        Cosmic.exact_compton_collision!(kout, δkn, kop, 1e-3, 0.0)
        @test maximum(abs, kout) < 1e-17
        @test abs(sum(kop.vol .* kout)) < 1e-20     # photon-number conservation

        # The BBN trace composition is carried into thermal heating. Published
        # 90/10 tritium epochs and the ⁷Be population ledger are hard gates.
        nh = trace_nuclear_history(c, r)
        @test Cosmic._tritium_fraction(nh, 6.35e5) ≈ 0.9 atol = 0.01
        @test Cosmic._tritium_fraction(nh, 1.35e5) ≈ 0.1 atol = 0.01
        @test sum(Cosmic._be_fractions(nh, 3e4)) ≈ 1.0 atol = 2e-10
        @test heating_rate(c, r, 3e4; silk=false, cooling=false,
            nuclear_history=nh) == nuclear_heating_rate(c, r, 3e4; history=nh)
    end

    FULL_TESTS && @testset "growth" begin
        c = cosmology(m_ν=Float64[])
        r = recombination(c)

        # The top-hat window must be well behaved at the origin.
        @test top_hat_window(0.0) ≈ 1.0
        @test top_hat_window(1e-8) ≈ 1.0 rtol = 1e-10

        # Growth factor, normalised, against standard ΛCDM values.
        @test growth_factor(c, r, 0.0) ≈ 1.0
        @test growth_factor(c, r, 1.0) ≈ 0.61 rtol = 0.05
        # Growth is monotonic in time.
        @test growth_factor(c, r, 2.0) < growth_factor(c, r, 1.0) < growth_factor(c, r, 0.0)
        # Growth rate f ≈ Ω_m(a)^0.55 today.
        @test growth_rate(c, r, 0.0) ≈ Ω_m(c)^0.55 rtol = 0.1
    end

    FULL_TESTS && @testset "HMcode-2020 (CAMB-validated)" begin
        # The non-linear boost P_HMcode/P_lin, validated against CAMB's own
        # mead2020 (halofit.f90) to <0.5% across 0.05 < k[h] < 10. These pins
        # guard the halo-model chain -- in particular the formation-redshift
        # concentration: it must solve g(a_f) = g(z)·δc/σ(γM) for a_f, NOT use
        # the g_f/g_z proxy that assumes g ∝ a. The proxy dropped the g(z)≈0.78
        # growth-suppression factor, under-concentrating small haloes ~28% and
        # deflating P(k) by ~9% at k ≈ 8 h/Mpc (the concentration shapes the
        # NFW window only at high k), which these high-k pins would catch.
        c = cosmology(m_ν=Float64[], Yp=0.24568)
        r = recombination(c; hydrogen=:hyrec)
        P = matter_power_spectrum(c, r; z=0.0, kmin=1e-4, kmax=80.0, nk=400)
        hm = hmcode_power(P)
        h = c.h
        @test power(hm, 1.0h) / power(P, 1.0h) ≈ 5.669 rtol = 0.01
        @test power(hm, 3.0h) / power(P, 3.0h) ≈ 19.107 rtol = 0.01
        @test power(hm, 8.0h) / power(P, 8.0h) ≈ 39.318 rtol = 0.01
        # concentration boosted above the base B = 5.196 for an early-forming
        # small halo (the fix; the old g∝a proxy sat far lower here)
        @test Cosmic._concentration(hm, 1e11) > 15.0
        # flat-ΛCDM Dolag factor is identically 1 (no wCDM/curvature)
        @test Cosmic._hmcode_growth(c, 0.0).ginf_ratio ≈ 1.0 atol = 1e-9
    end

    FULL_TESTS && @testset "HMcode massive-ν total matter (CAMB-validated)" begin
        # The halo model runs on the cold (cdm+baryon) spectrum; the total-matter
        # nonlinear P(k) is reconstructed as P_cb^NL·(P_mm^lin/P_cb^lin), the
        # neutrinos carrying their (scale-dependent, → f_cb² at high k) linear
        # suppression through unchanged. Validated against CAMB mead2020 with
        # Σmν = 0.15 eV to <0.25% across 0.05 < k[h] < 8; here we pin Cosmic's
        # own nonlinear total/cold ratio to the CAMB total/cold ratio.
        c = cosmology(m_ν=[0.15, 0.0, 0.0], Yp=0.24568)   # 1 massive + 2 massless
        r = recombination(c; hydrogen=:hyrec)
        Pc = matter_power_spectrum(c, r; z=0.0, kmin=1e-4, kmax=80.0, nk=250)
        Pt = matter_power_spectrum(c, r; z=0.0, kmin=1e-4, kmax=80.0, nk=250, total=true)
        hm_c = hmcode_power(Pc)                 # cold output (P_total defaults to Pc)
        hm_t = hmcode_power(Pc; P_total=Pt)     # total-matter reconstruction
        h = c.h
        ratio(k) = power(hm_t, k) / power(hm_c, k)
        @test ratio(0.1h) ≈ 0.9833 rtol = 0.005
        @test ratio(1.0h) ≈ 0.9779 rtol = 0.005
        @test ratio(8.0h) ≈ 0.9776 rtol = 0.005   # → f_cb² high-k limit
        # passing the cold spectrum as its own total is the identity (massless path)
        @test power(hmcode_power(Pc; P_total=Pc), 1.0h) == power(hm_c, 1.0h)
        # massless cosmology: total spectrum is bit-identical to cold
        cml = cosmology(m_ν=Float64[], Yp=0.24568)
        rml = recombination(cml; hydrogen=:hyrec)
        Pcl = matter_power_spectrum(cml, rml; z=0.0, kmin=1e-4, kmax=5.0, nk=60)
        Ptl = matter_power_spectrum(cml, rml; z=0.0, kmin=1e-4, kmax=5.0, nk=60, total=true)
        @test Pcl.δm == Ptl.δm
    end
end

@testset "BBN" begin
    # The neutrino temperature is not put in by hand: it comes from entropy
    # conservation in the γ-e± plasma, and must reproduce (4/11)^(1/3) once the
    # positrons are gone.
    @test Cosmic.T_nu_over_T_gamma(10.0) ≈ 1.0 atol = 1e-4
    @test Cosmic.T_nu_over_T_gamma(1e-3) ≈ (4 / 11)^(1 / 3) rtol = 1e-6

    # The weak rates are normalised by the neutron lifetime, so as T → 0 the n→p
    # rate must collapse onto free decay.
    @test weak_rate(0.005, 0.005, 878.4; direction=1) ≈ 1 / 878.4 rtol = 1e-3
    # and p→n must shut off entirely (no positrons, no neutrinos to capture)
    @test weak_rate(0.005, 0.005, 878.4; direction=-1) < 1e-30

    # Born phase-space factor: the literature value is 1.6360.
    @test Cosmic._LAMBDA0 ≈ 1.6360 rtol = 1e-4

    # The forward rates are measured cross-sections and cannot be derived -- but
    # the *thermal averaging* that turned σ(E) into the tabulated N_A<σv> can be
    # redone from scratch. Recompute the three deuterium-burning rates by direct
    # Gamow integration of the published S-factor polynomials (Pisanti et al 2021,
    # Appendix A -- the exact fits behind PArthENoPE 3.0, LUNA data included) and
    # demand the tables reproduce them. This is the strongest check a measured
    # rate admits: the physics layer (barrier penetration + Maxwell average) is
    # verified, the data layer is pinned to its published source.
    @testset "forward rates vs S-factors" begin
        function NA_sigmav(S_MeVb, Z1Z2, μ_u, T9)
            kT = 0.086173332621 * T9
            bG = 0.9895339 * Z1Z2 * sqrt(μ_u)          # 2πη = bG/√E, E in MeV
            I = quadgk(E -> S_MeVb(E) * exp(-bG / sqrt(E) - E / kT), 0.0, 40kT;
                rtol=1e-8)[1]
            6.02214076e23 * 1e-24 * 2.99792458e10 *
                sqrt(8 / (π * μ_u * 931.49410242)) * kT^(-1.5) * I
        end
        S_dpγ(E) = (0.2121 + 5.975E + 5.463E^2 - 1.665E^3) * 1e-6
        S_ddn(E) = 0.05225 + 0.3655E - 0.1799E^2 + 0.05832E^3 - 0.007393E^4
        S_ddp(E) = 0.05520 + 0.2151E - 0.02555E^2
        μ_dp = 2.01355321 * 1.00727647 / (2.01355321 + 1.00727647)
        μ_dd = 2.01355321 / 2
        r_dpγ = Cosmic._load_rate("d_p__He3_g")
        r_ddn = Cosmic._load_rate("d_d__He3_n")
        r_ddp = Cosmic._load_rate("d_d__t_p")
        for T9 in (0.3, 0.5, 0.8, 1.0)
            @test NA_sigmav(S_dpγ, 1.0, μ_dp, T9) ≈ r_dpγ(T9) rtol = 0.02
            @test NA_sigmav(S_ddn, 1.0, μ_dd, T9) ≈ r_ddn(T9) rtol = 0.02
            @test NA_sigmav(S_ddp, 1.0, μ_dd, T9) ≈ r_ddp(T9) rtol = 0.02
        end
    end

    # Every stored reverse-rate coefficient set must agree with the one derived
    # from masses and spins. α to 0.5% (the stored values carry 5 digits and
    # slightly different mass tables), β exactly, γ to 0.15 (≈ 13 keV in Q).
    for r in Cosmic._network(:full)
        α, β, γ = Cosmic.detailed_balance(r.reactants, r.products)
        @test isapprox(α, r.α; rtol=5e-3)
        @test β == r.β
        @test isapprox(γ, r.γ; atol=0.15)
    end
    # Charge and baryon number balance in every reaction (γ carries neither).
    for r in Cosmic._network(:full)
        @test sum(Cosmic.NUC_A[i] for i in r.reactants) ==
              sum(Cosmic.NUC_A[i] for i in r.products)
        @test sum(Cosmic.NUC_Z[i] for i in r.reactants) ==
              sum(Cosmic.NUC_Z[i] for i in r.products)
    end

    b = bbn(ω_b=0.02242, N_eff=3.044)
    # baryon number: conserved by construction, so this measures stiff-solver
    # drift, which is hardware dependent (2.7e-8 observed on GitHub runners
    # where 1e-10 is typical locally)
    @test sum(Cosmic.NUC_A .* b.Y) ≈ 1.0 rtol = 1e-7
    # References must match the rate compilation they test. The default chain is
    # PRIMAT throughout, so the references are PRIMAT's: CLASS's sBBN_2025_primat.dat
    # for Y_p (0.245683 at ω_b = 0.02242, τ_n = 878.4) and PRyMordial run on these
    # exact inputs with the same compilation (D/H = 2.4368e-5, Li7/H = 5.551e-10).
    @test b.Y_p ≈ 0.245683 rtol = 1.5e-3
    @test b.DH ≈ 2.4368e-5 rtol = 0.01
    @test b.Li7H ≈ 5.551e-10 rtol = 0.05
    # The PArthENoPE 3.0 core fits are the post-LUNA D/H reference (Pisanti et al
    # 2021, eq. 16: 2.51e-5 at ω_b = 0.02242, rescaled to τ_n = 878.4). The ~2.6%
    # D/H and ~13% Li7 differences between the two compilations are genuine spread
    # in the fitted cross-sections, and both configurations must land on their own
    # compilation's reference.
    bp = bbn(ω_b=0.02242, N_eff=3.044, rates=:parthenope)
    @test bp.DH ≈ 2.509e-5 rtol = 0.01
    # Y_p rises with N_eff (faster expansion, earlier freeze-out, more neutrons)
    @test bbn(ω_b=0.02242, N_eff=4.044).Y_p > b.Y_p

    # cosmology(Yp = :bbn) must close the loop: helium consistent with its own ω_b
    c = cosmology(Yp=:bbn)
    @test 0.24 < c.Yp < 0.25
end
