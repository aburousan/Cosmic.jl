# The Einstein-Boltzmann system as implemented

Conventions follow Ma & Bertschinger (1995) in the conformal Newtonian gauge,

    ds² = a²(η) [ -(1 + 2ψ) dη² + (1 - 2φ) δᵢⱼ dxⁱ dxʲ ]

with ψ the Newtonian potential and φ the spatial curvature perturbation. Overdots
are derivatives with respect to conformal time η (Mpc). The code integrates in
`x = ln a`, so every equation below is divided through by ℋ = aH.

The Thomson scattering rate is κ̇ = a nₑ σ_T > 0 (`-τ_dot` in the thermodynamics
layer).

## Fluid variables

δ = δρ/ρ̄, θ = ikʲvⱼ (velocity divergence), σ = anisotropic stress.
For the photon hierarchy, δ_γ = F_γ0, θ_γ = (3k/4) F_γ1, σ_γ = F_γ2/2.

## Metric: solved algebraically, not evolved

Combining the energy and momentum constraints (M&B eq. 23a, 23b) eliminates φ̇ and
gives φ directly, which is both cheaper and more stable than carrying φ as a
dynamical variable:

    φ  = -(3/2) (H₀²a²/k²) [ Δ + 3ℋ Θ_tot / k² ]
    ψ  = φ - (9/2) (H₀²a²/k²) Σ_σ
    φ̇  = -ℋψ + (3/2) (H₀²a²/k²) Θ_tot

where, summing over species i,

    Δ      = Σ (ρᵢ/ρ_c0) δᵢ
    Θ_tot  = Σ (1 + wᵢ)(ρᵢ/ρ_c0) θᵢ
    Σ_σ    = Σ (1 + wᵢ)(ρᵢ/ρ_c0) σᵢ

The energy constraint is then redundant and is used as an error monitor.

## Species

CDM (pressureless, no anisotropic stress):

    δ̇_c = -θ_c + 3φ̇
    θ̇_c = -ℋ θ_c + k²ψ

Baryons (pressureless but for the sound speed c_s², dragged by the photons):

    δ̇_b = -θ_b + 3φ̇
    θ̇_b = -ℋ θ_b + c_s² k² δ_b + k²ψ + (4ρ̄_γ)/(3ρ̄_b) κ̇ (θ_γ - θ_b)

Photons, with polarization. Π ≡ F_γ2 + G_γ0 + G_γ2:

    δ̇_γ  = -(4/3) θ_γ + 4φ̇
    θ̇_γ  = k²(δ_γ/4 - σ_γ) + k²ψ + κ̇(θ_b - θ_γ)
    Ḟ_γ2 = (8/15) θ_γ - (3/5) k F_γ3 + κ̇[ -F_γ2 + Π/10 ]
    Ḟ_γℓ = k/(2ℓ+1) [ ℓ F_γ,ℓ₋₁ - (ℓ+1) F_γ,ℓ₊₁ ] - κ̇ F_γℓ            (ℓ ≥ 3)
    Ġ_γℓ = k/(2ℓ+1) [ ℓ G_γ,ℓ₋₁ - (ℓ+1) G_γ,ℓ₊₁ ]
           + κ̇[ -G_γℓ + (Π/2)(δ_ℓ0 + δ_ℓ2/5) ]

Massless neutrinos: the photon hierarchy with κ̇ = 0 (they free-stream).

    δ̇_ν  = -(4/3) θ_ν + 4φ̇
    θ̇_ν  = k²(δ_ν/4 - σ_ν) + k²ψ
    Ḟ_νℓ = k/(2ℓ+1) [ ℓ F_ν,ℓ₋₁ - (ℓ+1) F_ν,ℓ₊₁ ]                     (ℓ ≥ 2)

Hierarchies are truncated (M&B eq. 65) by

    Ḟ_ℓmax = k F_ℓmax₋₁ - (ℓmax + 1)/η · F_ℓmax  [ - κ̇ F_ℓmax for photons ]

which absorbs outgoing power instead of reflecting it back down the hierarchy.

## Adiabatic initial conditions

Derived, not quoted, so they can be checked. Deep in radiation domination and far
outside the horizon (kη ≪ 1), write f_ν = ρ_ν/(ρ_γ + ρ_ν) and normalize to ψ:

Photon/neutrino density: from δ̇_γ = -(4/3)θ_γ + 4φ̇ with φ̇ ≈ 0 and θ small, δ_γ is
constant. The adiabatic condition ties it to ψ:

    δ_γ = δ_ν = -2ψ,      δ_c = δ_b = (3/4) δ_γ = -(3/2) ψ

Velocities: in tight coupling θ_b = θ_γ, and

    θ̇_γ ≈ k²(δ_γ/4 + ψ) = k²(-ψ/2 + ψ) = k²ψ/2   ⟹   θ_γ = (1/2) k² η ψ

For CDM, θ̇_c = -ℋθ_c + k²ψ with ℋ = 1/η in RD gives (ηθ_c)' = k²ψη, hence the same
θ_c = (1/2) k² η ψ. All four velocities agree at leading order, as adiabaticity
requires.

Neutrino anisotropic stress: Ḟ_ν2 ≈ (8/15)θ_ν gives

    σ_ν = F_ν2/2 = (1/15) k² η² ψ

Consistency check on φ: the anisotropic-stress constraint k²(φ-ψ) = 12πG a²(ρ̄+P̄)σ,
with 4πGa²ρ_tot = (3/2)ℋ² = (3/2)/η² in RD, gives

    k²(φ - ψ) = (4/3)(9/2)(f_ν/η²) σ_ν = 6 f_ν σ_ν / η² = (2/5) f_ν k² ψ

    ⟹   φ = (1 + 2f_ν/5) ψ

which reproduces the standard result and confirms the set is self-consistent.

Normalization to the primordial curvature perturbation ℛ. Superhorizon in RD, with
w = 1/3 and φ̇ ≈ 0,

    ℛ = φ + (2/3)(ℋ⁻¹φ̇ + ψ)/(1+w) = φ + ψ/2 = ψ (3/2 + 2f_ν/5)

    ⟹   ψ = ℛ / (3/2 + 2 f_ν / 5)

The code sets ℛ = 1 per mode and applies the primordial spectrum
P_ℛ(k) = A_s (k/k_pivot)^{n_s-1} · 2π²/k³ afterwards, since the system is linear.

## Massive neutrinos

Massive neutrinos cannot be reduced to a fluid. A fluid has one velocity at each
point; a free-streaming gas of massive particles has particles of *different
momenta moving at different speeds* through the same point, and it is exactly that
spread that erases small-scale power. So the perturbation must be tracked as a
function of momentum, not just position: the state is Ψ_ℓ(q), a hierarchy per
momentum bin.

Write f(q, n̂, η) = f₀(q)[1 + Ψ(q, n̂, η)] with f₀(q) = 1/(e^q + 1), where q is the
comoving momentum in units of T_ν,0 and

    ε(q, a) = √(q² + a²m̃²),    m̃ = m c²/(k_B T_ν,0)

The hierarchy (Ma & Bertschinger eq. 100, conformal Newtonian):

    Ψ̇₀ = -(qk/ε) Ψ₁ - φ̇ · dlnf₀/dlnq
    Ψ̇₁ = (qk/3ε)(Ψ₀ - 2Ψ₂) - (εk/3q) ψ · dlnf₀/dlnq
    Ψ̇_ℓ = (qk/(2ℓ+1)ε) [ ℓΨ_{ℓ-1} - (ℓ+1)Ψ_{ℓ+1} ]        (ℓ ≥ 2)

with, for the Fermi-Dirac f₀,

    dlnf₀/dlnq = -q / (1 + e^{-q})

The factor q/ε is the particle's velocity. When q ≫ am̃ it is 1 and the equations
reduce to the massless neutrino hierarchy; when q ≪ am̃ it goes to zero and the
free-streaming term switches off, which is the neutrinos becoming non-relativistic
and starting to cluster.

### Back-reaction on the metric

The moments are recovered by integrating over q, weighted by f₀:

    δρ_ν     ∝ (1/a⁴) ∫ dq q² ε   f₀ Ψ₀
    δP_ν     ∝ (1/3a⁴) ∫ dq q² (q²/ε) f₀ Ψ₀
    (ρ̄+P̄)θ_ν ∝ (k/a⁴) ∫ dq q³     f₀ Ψ₁
    (ρ̄+P̄)σ_ν ∝ (2/3a⁴) ∫ dq q⁴/ε  f₀ Ψ₂

normalised against ρ̄_ν ∝ (1/a⁴)∫dq q² ε f₀ — which is exactly the F(y) integral the
species layer already uses for the background, so the two are guaranteed consistent.

### Initial conditions

Deep in radiation domination the neutrinos are ultra-relativistic and adiabatic,
so their distribution perturbation is just the massless answer reshaped by f₀:

    Ψ₀ = -(1/4) δ_ν      · dlnf₀/dlnq
    Ψ₁ = -(ε/3qk) θ_ν    · dlnf₀/dlnq
    Ψ₂ = -(1/2) σ_ν      · dlnf₀/dlnq

Check of the monopole: in the relativistic limit ε → q, integrating by parts gives
∫dq q²ε f₀ (dlnf₀/dlnq) = -4∫dq q³ f₀, while ρ̄_ν ∝ ∫dq q³ f₀. So
δρ_ν = -(1/4)δ_ν · (-4ρ̄_ν) = δ_ν ρ̄_ν, as required.
