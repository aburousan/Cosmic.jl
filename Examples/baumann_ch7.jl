#
# Baumann, *Cosmology* (2022), chapter 7 -- the CMB figures, from Cosmic.jl.
#
#   Fig. 7.16   baryon loading of the Sachs-Wolfe transfer function |F(k)|²
#   Fig. 7.10   the photon transfer function |δ_γ(η,k)| over the (η, 1/k) plane
#
# Run with:  julia --project=examples examples/baumann_ch7.jl
#
using Cosmic
using Plots
using Printf

gr(size=(900, 420), dpi=150, framestyle=:box, guidefontsize=10, tickfontsize=8)

c = cosmology(m_ν=Float64[])
r = recombination(c)
bg = BackgroundCache(c, r)
h = c.h

z_rec = z_star(r)
η_rec = conformal_time(c, scale_factor(z_rec))
rs_rec = sound_horizon(c, z_rec)
R_rec = R_baryon(c, scale_factor(z_rec))
kD_rec = damping_scale(c, r, scale_factor(z_rec))

@printf("z_* = %.1f   η_* = %.1f Mpc   r_s = %.2f Mpc   R_* = %.3f   k_D = %.4f /Mpc\n",
    z_rec, η_rec, rs_rec, R_rec, kD_rec)

# ===========================================================================
# Fig. 7.16 -- how baryons load the acoustic oscillator
# ===========================================================================
# The photon-baryon fluid oscillates, but the baryons weigh it down: the zero
# point of the oscillation is displaced from 0 to -R·ψ. So the transfer function
# is not cos(kr_s) but cos(kr_s) - shift, and once you *square* it the compressions
# (which are pushed away from zero) beat the rarefactions (which are pulled toward
# it). That is the odd/even peak-height alternation in the CMB, and it is what
# makes the CMB able to weigh the baryons at all.
#
# Left panel: the idea in its purest form, no damping.
# Right panel: the same thing inside Weinberg's solution (Baumann 7.113), where
# Silk damping and the WKB factor (1+R)^{-1/4} partly hide it.

x = range(0.0, 6.5; length=1200)          # x = k r_s / π
kr = x .* π

p_left = plot(x, cos.(kr) .^ 2, ls=:dot, color=:black, lw=1.4,
    label="[cos(kr_s)]²", xlabel="k r_s,* / π", ylabel="|F_*(k)|²",
    xlims=(0, 6.5), ylims=(0, 2), legend=:topright,
    title="idealised: shifting the zero point")
plot!(p_left, x, (cos.(kr) .- 0.15) .^ 2, color=:black, lw=2,
    label="[cos(kr_s) - 0.15]²")

# Weinberg's ansatz, Baumann (7.113):
#
#   F_*(k) = (1/5) [ S(κ)/(1+R_*)^{1/4} · cos(k r_s,* + θ(κ))  -  3 R_* T(κ) ]
#
# with κ = √2 k/k_eq, and S, T, θ the fitting functions (7.115) obtained by solving
# the coupled photon-baryon + CDM equations numerically. S is the radiation-driving
# boost (S → 5 well inside the horizon), T the Meszaros suppression, and θ the phase
# shift, which is what fixes the peak positions -- Weinberg's agreement with them is,
# in his words, "almost embarrassingly good".
#
# Silk damping is not in (7.113); Baumann adds it (below 7.133) by replacing
# S(k) → e^{-k²/k_S²} S(k). Note it multiplies *only* the oscillating term, never the
# -3R_*T offset -- diffusion erases the wave, not the baryon loading. That is exactly
# why the effect survives at all: damping kills the oscillation but leaves the shifted
# zero point behind.
S_w(κ) = ((1 + (1.209κ)^2 + (0.5116κ)^4 + sqrt(5) * (0.1657κ)^6) /
          (1 + (0.9459κ)^2 + (0.4249κ)^4 + (0.167κ)^6))^2
T_w(κ) = log(1 + (0.124κ)^2) / (0.124κ)^2 *
         sqrt((1 + (1.257κ)^2 + (0.4452κ)^4 + (0.2197κ)^6) /
              (1 + (1.606κ)^2 + (0.8568κ)^4 + (0.3927κ)^6))
θ_w(κ) = sqrt(((1.1547κ)^2 + (0.5986κ)^4 + (0.2578κ)^6) /
              (1 + (1.723κ)^2 + (0.8707κ)^4 + (0.4581κ)^6 + (0.2204κ)^8))

a_eq = scale_factor(z_equality(c))
k_eq = a_eq * H_Mpc(c, a_eq)

function F_weinberg(k, R; damp=true)
    κ = sqrt(2) * k / k_eq
    D = damp ? exp(-(k / kD_rec)^2) : 1.0
    (D * S_w(κ) / (1 + R)^0.25 * cos(k * rs_rec + θ_w(κ)) - 3R * T_w(κ)) / 5
end

ks = kr ./ rs_rec
p_right = plot(x, [F_weinberg(k, 0.0)^2 for k in ks], ls=:dot, color=:black, lw=1.4,
    label="R_* = 0", xlabel="k r_s,* / π", ylabel="|F_*(k)|²",
    xlims=(0, 6.5), ylims=(0, 0.6), legend=:topright,
    title="Weinberg (7.113), R_* = $(round(R_rec, digits=2)) from Cosmic")
plot!(p_right, x, [F_weinberg(k, R_rec)^2 for k in ks], color=:black, lw=2,
    label="R_* = $(round(R_rec, digits=2))")

savefig(plot(p_left, p_right, layout=(1, 2)), "examples/fig7_16_baryon_loading.png")
println("  -> examples/fig7_16_baryon_loading.png")

# ---------------------------------------------------------------------------
# The same thing, but from the actual Boltzmann solve rather than the WKB formula.
# F_* = Θ₀ + ψ evaluated at last scattering: no fitted envelope, no assumed phase.
# ---------------------------------------------------------------------------
println("computing the true |Θ₀+ψ|² at last scattering ...")
kk = range(1e-3, 6.5π / rs_rec; length=500)
Feff = zeros(length(kk))
Threads.@threads for i in eachindex(kk)
    p = solve_perturbations(c, bg, kk[i]; lmax_γ=25, lmax_ν=32)
    a = scale_factor(z_rec)
    Feff[i] = δ_photon(p, a) / 4 + Ψ(p, a)          # Θ₀ + ψ
end

p_true = plot(kk .* rs_rec ./ π, Feff .^ 2, color=:crimson, lw=1.6,
    xlabel="k r_s,* / π", ylabel="|Θ₀ + ψ|²  at z = z_*",
    xlims=(0, 6.5), legend=false, size=(900, 380),
    title="the real thing: Cosmic.jl Boltzmann solve, no WKB, no fitted envelope")
savefig(p_true, "examples/fig7_16_boltzmann.png")
println("  -> examples/fig7_16_boltzmann.png")

# ===========================================================================
# Fig. 7.10 -- the photon transfer function across (η, 1/k)
# ===========================================================================
# One panel that contains the whole story. Above the sound horizon, δ_γ is frozen
# at its initial value (the dark band, top left). Below it, the mode has had time
# to oscillate, and the diagonal stripes are the sound waves. Below the damping
# scale, diffusion has erased them (the white wedge, bottom). The dashed line is
# recombination, after which the photons free-stream and the pattern is frozen into
# the sky we see.
println("Fig 7.10: |δ_γ(η,k)| map (this one takes a couple of minutes) ...")

nk, nη = 220, 400
kinv = exp.(range(log(0.7), log(120.0); length=nk))     # 1/k in Mpc/h, as Baumann plots
ks_map = h ./ kinv                                       # -> 1/Mpc
ηs = exp.(range(log(10.0), log(3000.0); length=nη))

# η -> a, by inverting the (monotonic) conformal time
as_grid = exp.(range(log(1e-8), 0.0; length=4000))
ηs_grid = [conformal_time(c, a) for a in as_grid]
a_of_η(η) = as_grid[searchsortedfirst(ηs_grid, clamp(η, ηs_grid[1], ηs_grid[end]))]

M = zeros(nk, nη)
Threads.@threads for i in 1:nk
    p = solve_perturbations(c, bg, ks_map[i]; lmax_γ=25, lmax_ν=32)
    for j in 1:nη
        M[i, j] = abs(δ_photon(p, a_of_η(ηs[j])))
    end
end

p10 = heatmap(ηs, kinv, clamp.(M, 0, 4.0),
    xscale=:log10, yscale=:log10, color=cgrad(:grays, rev=true),
    xlabel="η  [Mpc]", ylabel="k⁻¹  [Mpc/h]",
    xlims=(10, 3000), ylims=(0.7, 120), clims=(0, 4),
    colorbar_title="|δ_γ(η,k)|", size=(820, 560),
    title="Baumann Fig. 7.10 — photon transfer function")

# The three scales that carve the plot up. All are computed, none are drawn by hand.
ηl = exp.(range(log(10.0), log(3000.0); length=300))
al = [a_of_η(η) for η in ηl]
plot!(p10, ηl, [h / (a * H_Mpc(c, a)) for a in al], color=:black, lw=2.5, label="Hubble scale")
plot!(p10, ηl, [h * sound_horizon(c, redshift(a)) for a in al], color=:black, lw=2.5,
    ls=:solid, label="sound horizon")
plot!(p10, ηl, [h / damping_scale(c, r, a) for a in al], color=:black, lw=2.5,
    label="damping scale")
vline!(p10, [η_rec], ls=:dash, color=:black, lw=1.5, label="recombination")

savefig(p10, "examples/fig7_10_transfer_map.png")
println("  -> examples/fig7_10_transfer_map.png")
println("\ndone.")
