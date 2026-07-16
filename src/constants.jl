"""
Physical constants and the internal unit system.

Cosmic works internally in "Mpc units": lengths and times are both measured in
Mpc with c = 1, so H has dimensions of 1/Mpc and wavenumbers k are in 1/Mpc.
This is what CLASS and CAMB do, and it keeps every number in the Boltzmann
hierarchy within a few orders of magnitude of unity.

Values are CODATA 2018 / PDG 2022. Constants that are exact by SI definition
are marked as such.
"""

module Constants

# SI base
const c_SI      = 2.99792458e8        # m/s, exact by definition
const G_SI      = 6.67430e-11         # m^3 kg^-1 s^-2
const ħ_SI      = 1.054571817e-34     # J s, exact (from exact h)
const k_B_SI    = 1.380649e-23        # J/K, exact by definition
const eV_SI     = 1.602176634e-19     # J, exact by definition
const σ_T_SI    = 6.6524587321e-29    # m^2, Thomson cross section
const m_e_SI    = 9.1093837015e-31    # kg
const m_p_SI    = 1.67262192369e-27   # kg
const m_H_SI    = 1.67353284e-27      # kg, neutral hydrogen atom
const m_He_SI   = 6.64647907e-27      # kg, neutral helium-4 atom

# Radiation constant a_rad = 4σ_SB/c = π²k_B⁴/(15 ħ³c³)
const a_rad_SI  = 7.565733e-16        # J m^-3 K^-4

# Length / time
const Mpc_SI    = 3.0856775814913673e22  # m
const Gyr_SI    = 3.1556952e16           # s (Julian year × 1e9)
const km_SI     = 1.0e3

# Electron mass in temperature-friendly units
const m_e_eV    = m_e_SI * c_SI^2 / eV_SI     # 510998.95 eV

# --- Atomic data ------------------------------------------------------------
#
# Energy levels are given as inverse wavelengths in m⁻¹, which is how they are
# measured spectroscopically and how RECFAST/CAMB carry them. Converting to eV
# and back would only lose digits. Use `T_of_wavenumber` to turn one into the
# equivalent temperature, which is the form the Boltzmann factors want.

# Hydrogen (Drake 1993 / NIST)
const L_H_ion    = 1.096787737e7      # H ionization (13.598 eV)
const L_H_alpha  = 8.225916453e6      # H Lyman-α, averaged over the 2 levels

# Helium (Drake 1993; Drake & Morton 2007)
const L_He1_ion  = 1.98310772e7       # He I  -> He II
const L_He2_ion  = 4.389088863e7      # He II -> He III
const L_He_2s    = 1.66277434e7       # He I 2¹S
const L_He_2p    = 1.71134891e7       # He I 2¹P
const L_He_2Pt   = 1.690871466e7      # He I 2³P  (triplet)
const L_He_2St   = 1.5985597526e7     # He I 2³S  (triplet)
const L_He2St_ion = 3.8454693845e6    # ionization from 2³S

# Einstein A coefficients (1/s)
const A2P_s = 1.798287e9              # He I 2¹P₁ -> 1¹S₀  (Morton, Wu & Drake 2006)
const A2P_t = 177.58                  # He I 2³P₁ -> 1¹S₀  (Lach & Pachucki 2001)
                                      # the semi-forbidden intercombination line

# Photoionization cross sections at the 2P thresholds, m² (Hummer & Storey 1998)
const σ_He_2Ps = 1.436289e-22
const σ_He_2Pt = 1.484872e-22

# Two-photon decay rates (1/s)
const Λ_2s1s_H   = 8.2245809          # H  2s -> 1s
const Λ_2s1s_He1 = 51.3               # He 2s -> 1s (Dalgarno)

# Helium recombination fit (Verner & Ferland form), triplet channel
const a_trip = 10.0^(-16.306)
const b_trip = 0.761

"""
    T_of_wavenumber(L)

Convert an inverse wavelength L (m⁻¹) into the equivalent temperature hcL/k_B
in K. This is the form every Boltzmann and Saha exponent wants.
"""
T_of_wavenumber(L) = 2π * ħ_SI * c_SI * L / k_B_SI

# Kept for readability where an energy in eV reads more naturally.
const E_ion_H   = 13.605693122994
const E_2s_H    = E_ion_H / 4

# Boltzmann-suppression helper: 1 K expressed in eV
const K_to_eV   = k_B_SI / eV_SI      # 8.617333e-5 eV/K

# ζ(3), used for photon and baryon number densities
const ζ3        = 1.2020569031595943

# --- Derived quantities in the internal (Mpc, c=1) system --------------------

"Hubble constant H0 = 100h km/s/Mpc expressed in 1/Mpc."
H0_in_invMpc(h) = 100.0 * h * km_SI / c_SI          # = h / 2997.92458

"Hubble time 1/H0 in Gyr."
hubble_time_Gyr(h) = Mpc_SI / (100.0 * h * km_SI) / Gyr_SI

"Hubble distance c/H0 in Mpc."
hubble_distance_Mpc(h) = c_SI / (100.0 * h * km_SI)   # = 2997.92458/h

"Critical density today, kg/m^3, for a given h."
ρ_crit_SI(h) = 3.0 * (100.0 * h * km_SI / Mpc_SI)^2 / (8.0 * π * G_SI)

"""
    Ω_γ(Tcmb, h)

Photon density parameter from the CMB temperature. This is
Ω_γ = a_rad T⁴ / (ρ_crit c²), which evaluates to ≈ 2.4728e-5/h² at T = 2.7255 K.
"""
Ω_γ(Tcmb, h) = a_rad_SI * Tcmb^4 / (ρ_crit_SI(h) * c_SI^2)

"Neutrino-to-photon temperature ratio after e± annihilation, (4/11)^(1/3)."
const T_ν_over_T_γ = (4 / 11)^(1 / 3)

"Thomson cross section in Mpc²."
const σ_T_Mpc2 = σ_T_SI / Mpc_SI^2

"""
Conversion factor for the baryon number density: given Ω_b h², the hydrogen
number density today in 1/Mpc³ is `n_H0 = (1 - Yp) * Ω_b * ρ_crit / m_H`.
"""
n_H_today(Ω_b, h, Yp) = (1 - Yp) * Ω_b * ρ_crit_SI(h) / m_H_SI * Mpc_SI^3  # 1/Mpc^3

end # module Constants
