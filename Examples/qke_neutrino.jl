# Running a cosmology on a solved neutrino quantum-kinetic spectrum.
#
# The standard treatment gives every neutrino an equilibrium Fermi-Dirac
# distribution and folds the non-instantaneous decoupling (N_eff = 3.044, not 3)
# into a rescaled temperature. Cosmic can instead use the real, momentum- and
# flavour-resolved spectrum from a solved quantum kinetic equation (FortEPiaNO),
# which carries the actual spectral distortions left by e+e- annihilation.
#
# Run with:  julia --project=Examples Examples/qke_neutrino.jl

using Cosmic

# The bundled Standard-Model three-flavour reference history (zero lepton
# asymmetry). Point load_flavor_density_history at another directory for a run
# with different masses, mixing, or a nonzero asymmetry.
h = standard_qke_history()

println("N_eff from the flavour trace: ", qke_neff(h))          # 3.0435

# Flavour resolution is the point: electron neutrinos decouple later (they keep
# annihilating with e+/e- through the charged current), so they end up hotter
# than the mu/tau flavours. The distortion delta = f/f_FD - 1 at momentum y = 3:
for (flav, name) in zip(1:3, ("e", "mu", "tau"))
    d = qke_distribution(h, flav; basis=:flavour, coordinate=:comoving)
    delta = Cosmic.nu_f0(d, 3.0) / (1 / (exp(3.0) + 1)) - 1
    println("  nu_$name  delta(y=3) = ", round(delta * 100, digits=3), " %")
end

# Feed the solved mass-basis spectra into a full cosmology. One TabulatedNu per
# mass eigenstate; the spectrum already carries the decoupling enhancement, so
# each species is booked at its own density and N_eff becomes a prediction of
# the spectrum rather than an input number.
dists = [qke_distribution(h, i; basis=:mass) for i in 1:3]
c_qke = cosmology(m_ν = [0.05, 0.05, 0.05], ν_dist = dists)
c_std = cosmology()   # equilibrium Fermi-Dirac, N_eff = 3.044

println("N_eff, solved QKE spectrum:  ", Cosmic._N_eff_of(c_qke))  # 3.0436
println("N_eff, equilibrium default:  ", Cosmic._N_eff_of(c_std))  # 3.0440

# From here c_qke behaves like any cosmology: recombination(c_qke),
# power_spectrum(c_qke, k), cmb_cls(c_qke), ... now propagate the true neutrino
# phase space through the background and the Boltzmann hierarchy.
