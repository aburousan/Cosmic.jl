#
# Reproduces the perturbation-chapter figures of Baumann, *Cosmology* (2022),
# straight out of Cosmic.jl.
#
#   Fig. 1  Φ(a)/Φᵢ for k = {0.001, 0.05, 1} h/Mpc   (Baumann Fig. 6.5)
#   Fig. 2  δ_γ and θ_b for k = 1 h/Mpc              (Baumann Fig. 6.9)
#   Fig. 3  the matter power spectrum and its k⁴ → (ln k)² shape
#
# Baumann's figures use an Einstein-de Sitter background: Ω_b = 0.05, Ω_c = 0.95,
# no cosmological constant. That is *not* our universe, but it is what makes the
# potential exactly constant in matter domination, which is the point of the plot.
#
# Run with:  julia --project=examples examples/baumann_figures.jl
#
using Cosmic
using Plots
using Printf

gr(size=(760, 520), dpi=150, legend=:topright, framestyle=:box,
    guidefontsize=11, tickfontsize=9, linewidth=2)

# --- Einstein-de Sitter, as Baumann uses --------------------------------------
# Build the species list by hand: the default `cosmology()` closes the budget with
# dark energy, and here we want none.
function eds_cosmology(; h=0.67, Ω_b=0.05, Ω_c=0.95, Tcmb=2.7255, N_eff=3.044)
    Ωγ = Cosmic.Constants.Ω_γ(Tcmb, h)
    ν = massless_neutrinos(N_eff, Ωγ)
    species = (Photons(Ωγ), Baryons(Ω_b), ColdDarkMatter(Ω_c), ν,
        CosmologicalConstant(0.0))          # EdS: no Λ
    Cosmology(h, Tcmb, 0.2454, 2.1e-9, 0.9665, 0.05, species)
end

c = eds_cosmology()
r = recombination(c)
bg = BackgroundCache(c, r)
h = c.h

@printf("Einstein-de Sitter: Ω_m = %.3f, Ω_γ = %.2e, Ω_ν = %.2e\n", Ω_m(c), Ω_γ(c), Ω_ν(c))
a_eq = scale_factor(z_equality(c))
k_eq = a_eq * H_Mpc(c, a_eq)
@printf("z_eq = %.0f,  a_eq = %.3e,  k_eq = %.4f /Mpc = %.4f h/Mpc\n\n",
    z_equality(c), a_eq, k_eq, k_eq / h)

# ---------------------------------------------------------------------------
# Fig. 1 — evolution of the gravitational potential
# ---------------------------------------------------------------------------
# Three modes bracketing k_eq: one that never enters the horizon before equality
# (Φ constant, bar the 9/10 drop), one entering right at equality, and one deep
# inside the horizon during radiation domination (Φ decays and rings).
println("Fig 1: Φ(a)/Φᵢ  ...")
ks_hMpc = [0.001, 0.05, 1.0]
labels = ["k < k_eq  (0.001 h/Mpc)", "k ≈ k_eq  (0.05 h/Mpc)", "k > k_eq  (1 h/Mpc)"]

p1 = plot(xscale=:log10, xlabel="a", ylabel="Φ / Φᵢ",
    xlims=(1e-7, 1), ylims=(-0.05, 1.05), legend=:left,
    title="Baumann Fig. 6.5 — potential evolution (EdS)")

for (kh, lab) in zip(ks_hMpc, labels)
    k = kh * h                                   # h/Mpc -> 1/Mpc
    p = solve_perturbations(c, bg, k; lmax_γ=24, lmax_ν=24)
    a0 = exp(p.sol.t[1])                         # start of this mode's solve
    as = exp.(range(log(max(a0, 1e-7)), 0.0; length=600))
    Φi = Φ(p, as[1])
    plot!(p1, as, [Φ(p, a) / Φi for a in as], label=lab)
end
savefig(p1, "examples/fig1_potential.png")
println("  -> examples/fig1_potential.png")

# ---------------------------------------------------------------------------
# Fig. 2 — acoustic oscillations and Silk damping
# ---------------------------------------------------------------------------
# The photon density oscillates as a sound wave in the photon-baryon fluid; the
# baryon velocity oscillates a quarter-period out of phase. Both damp away as
# photon diffusion erases the wave on scales below the mean free path.
println("Fig 2: δ_γ and θ_b for k = 1 h/Mpc  ...")
k = 1.0 * h
p = solve_perturbations(c, bg, k; lmax_γ=32, lmax_ν=32)
as = exp.(range(log(1e-6), log(1e-3); length=1500))
δγ = [δ_photon(p, a) for a in as]
θb = [θ_baryon(p, a) for a in as]
# Baumann plots the baryon *velocity*, v_b = θ_b/k, which is the quantity that sits
# naturally on the same axis as δ_γ: in tight coupling v_b ≈ (√3/4) δ_γ.
θb_scaled = θb ./ k

p2 = plot(xscale=:log10, xlabel="a", ylabel="δ_γ , θ_b",
    xlims=(1e-6, 1e-3), title="Baumann Fig. 6.9 — acoustic oscillations, k = 1 h/Mpc")
plot!(p2, as, δγ, label="photons  δ_γ", color=:gray, linewidth=1.6)
plot!(p2, as, θb_scaled, label="baryons  θ_b", color=:black, linewidth=2)
savefig(p2, "examples/fig2_acoustic.png")
println("  -> examples/fig2_acoustic.png")

# ---------------------------------------------------------------------------
# Fig. 3 — the matter power spectrum
# ---------------------------------------------------------------------------
# The shape Baumann sketches: P ∝ k on large scales (the primordial k^{n_s}), a
# turnover at k_eq, and P ∝ k^{-3}(ln k)² beyond it, because modes entering during
# radiation domination only grow logarithmically (the Meszaros effect).
println("Fig 3: P(k) and its asymptotic shape  ...")
cΛ = cosmology(m_ν=Float64[])                    # a realistic ΛCDM for this one
rΛ = recombination(cΛ)
P = matter_power_spectrum(cΛ, rΛ; nk=140, kmin=1e-4, kmax=10.0)

ks = exp.(range(log(1e-4), log(10.0); length=400))
Pk = [power(P, k) for k in ks]
a_eqΛ = scale_factor(z_equality(cΛ))
k_eqΛ = a_eqΛ * H_Mpc(cΛ, a_eqΛ)

p3 = plot(ks, Pk, xscale=:log10, yscale=:log10,
    xlabel="k  [1/Mpc]", ylabel="P(k)  [Mpc³]", label="Cosmic.jl",
    title="Matter power spectrum (ΛCDM)", legend=:bottomleft)
vline!(p3, [k_eqΛ], ls=:dash, color=:gray, label="k_eq = $(round(k_eqΛ,sigdigits=3))")
# asymptotes
kl = ks[ks.<k_eqΛ/2]
kh_ = ks[ks.>3k_eqΛ]
plot!(p3, kl, Pk[argmin(abs.(ks .- k_eqΛ / 4))] .* (kl ./ (k_eqΛ / 4)) .^ cΛ.n_s,
    ls=:dot, color=:red, label="∝ k^n_s")
plot!(p3, kh_, Pk[argmin(abs.(ks .- 3k_eqΛ))] .* (kh_ ./ (3k_eqΛ)) .^ (cΛ.n_s - 4) .*
           (log.(kh_ ./ (3k_eqΛ)) .+ 1) .^ 2,
    ls=:dot, color=:blue, label="∝ k^(n_s-4) (ln k)²")
savefig(p3, "examples/fig3_pk.png")
println("  -> examples/fig3_pk.png")

println("\ndone.")
