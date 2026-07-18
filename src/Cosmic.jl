"""
    Cosmic

A cosmology library in pure Julia, following Baumann's *Cosmology* and
Dodelson & Schmidt's *Modern Cosmology*.

The package is built in layers, each resting on the one below:

1. `constants.jl`      -- CODATA constants, internal Mpc/c=1 unit system
2. `species.jl`        -- energy components: γ, ν (massless and massive), CDM,
                          baryons, curvature, dark energy with arbitrary w(a)
3. `background.jl`     -- E(a), distances, times, conformal time
4. `thermodynamics.jl` -- recombination, x_e(z), optical depth, visibility
5. `perturbations.jl`  -- the Boltzmann hierarchy
6. `observables.jl`    -- P(k), σ₈, growth, CMB C_ℓ

The organising idea is that a cosmology is a *list of species*, each of which
knows only its own ρ(a). Everything above the species layer is generic in the
composition, so a new component costs one method rather than a new struct in
every file.
"""
module Cosmic

include("constants.jl")
using .Constants

include("species.jl")
include("background.jl")
include("curvature.jl")
include("modified_gravity.jl")
include("bbn.jl")
include("bbn_thermal.jl")
include("neff.jl")
include("hyrec.jl")
include("thermodynamics.jl")
include("nuclear_deposition.jl")
include("perturbations.jl")
include("gauges.jl")
include("observables.jl")
include("inflation.jl")
include("sigw.jl")
include("pbh.jl")
include("second_order.jl")
include("cmb.jl")
include("lensing.jl")
include("nonlinear.jl")
include("hmcode.jl")
include("eft.jl")
include("tensors.jl")
include("distortions.jl")
include("compton_kernel.jl")
include("thermalization.jl")

# Species and composition
export Species, ρ_over_ρc0, w, Ω0
export Photons, MasslessNeutrinos, MassiveNeutrinos, ColdDarkMatter, Baryons,
    Curvature, AbstractDarkEnergy, CosmologicalConstant, W0WaDarkEnergy,
    GeneralDarkEnergy
export photons_from_Tcmb, massless_neutrinos, massive_neutrino, fermi_dirac_ρ
export NeutrinoDistribution, FermiDirac, TabulatedNu, degenerate_Neff_factor, nu_ρ_integral
export warm_dark_matter, DecayingCDM, DecayRadiation, QuintessenceDE

# Background
export Cosmology, cosmology
export E, E_z, H, H_Mpc, ℋ, scale_factor, redshift
export get_species, get_all_species
export Ω_b, Ω_c, Ω_γ, Ω_k, Ω_m, Ω_de, Ω_ν, Ω_r
export hubble_distance, hubble_time
export comoving_distance, transverse_comoving_distance, angular_diameter_distance,
    luminosity_distance, distance_modulus, comoving_volume
export age, lookback_time, conformal_time, conformal_time_today
export z_equality, scale_factor_of_time, ρ_matter, ρ_radiation

# Thermodynamics
export RecombinationSolution, recombination, saha_x_e, saha_x_H, saha_state
export Reionization, no_reionization, reionization_x_e, τ_reio, z_reio_from_τ
export Injection, no_injection, injection_rate
export TraceNuclearHistory, trace_nuclear_history, nuclear_deposition_rate,
    nuclear_heating_rate
export x_e, T_matter, n_e, n_H_of_z, f_He, α_H, α_He
export optical_depth, drag_optical_depth, recombination_optical_depth, visibility, dτ_dz, τ_dot, z_star, z_drag
export R_baryon, sound_speed, sound_horizon, r_drag, r_star, θ_star
export bbn, bbn_mc, BBNSolution, Y_He_bbn, primordial_bbn, weak_rate, T_nu_over_T_gamma, NUCLIDES, detailed_balance
export thermal_weak_correction, generate_thermal_cache
export qed_pressure, qed_energy_density, F_ann, F_scat, neff, N_eff_sm, nu_temperatures, decoupling_history, scale_factor_at

# Perturbations
export PerturbationSolution, BackgroundCache, solve_perturbations, initial_conditions
export δ_cdm, θ_cdm, δ_baryon, θ_baryon, δ_photon, δ_neutrino, δ_matter, δ_matter_comoving, δ_matter_total_comoving, Θ, Φ, Ψ
export MassiveNuGrid

# Observables
export MatterPowerSpectrum, matter_power_spectrum, bao_k_grid, transfer, power,
    dimensionless_power, primordial_power, top_hat_window, σ_R, σ8,
    growth_factor, growth_rate
export HalofitSpectrum, halofit_power, nonlinear_scale
export HMcodeSpectrum, hmcode_power
export EFTSpectrum, eft_power, p22, p13

# CMB
export CMBSpectra, cmb_spectra, D_ℓ, SourceFunction, source_function, PrimordialMatrix
export InflatonSpectrum, inflaton_spectrum, inflation_background, tensor_power
export SIGWSpectrum, sigw_spectrum, Ω_gw_rd, sigw_monochromatic, sigw_ΔNeff
export sigw_spectrum_ng, sigw_hybrid, sigw_gnl_factor, sigw_omega_g_kernel
export sigw_Z_connected, sigw_C_connected, sigw_f4_reducible, sigw_gnl_reducible
export PBHAbundance, pbh_abundance, pbh_β, pbh_moments
export pbh_β_ng, pbh_compaction_variances
export LensedCMBSpectra, lensed_cmb_spectra, lensed_from_cls
export TensorPerturbationSolution, solve_tensor_perturbations, tensor_h, tensor_hdot
export TensorSourceFunction, TensorCMBSpectra, tensor_cmb_spectra

# Spectral distortions
export distortion_amplitudes, μ_distortion, y_distortion, g_distortion,
    distortion_spectrum, firas_chi2, heating_rate, silk_heating_rate,
    damping_scale, adiabatic_cooling_rate,
    FIRAS_μ_limit, FIRAS_y_limit, PIXIE_μ_target
export thermalization_grid, solve_thermalization, distortion_from_spectrum, H_dc, g_dc_exact, K_emission,
    thermalization_distortions, chluba_J_y, chluba_J_mu, chluba_J_vis
export single_electron_compton_kernel, thermal_compton_kernel,
    minimum_compton_momentum, compton_frequency_bounds,
    compton_detailed_balance_residual

# Gauges
export Gauge, Newtonian, Synchronous, Comoving, UniformDensity, SpatiallyFlat
export gauge_shift, δ_species, θ_species
export curvature_ℛ, curvature_ζ, bardeen_Φ, bardeen_Ψ

end # module
