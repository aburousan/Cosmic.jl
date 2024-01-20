using Cosmic
using Roots, ForwardDiff
using DifferentialEquations
using Plots, LaTeXStrings

c = cosmology()
ag = age(c,0)
sec_to_year(ag)
z_vals = range(0, 2.5, length=1000)
H_1z(z) = H(c,z)/(1+z)
H_vals = H_1z.(z_vals)
plot(z_vals,H_vals,lw=2.5, label="")
xlabel!(L"Redshift ($z$)")
ylabel!(L"$\frac{H(z)}{1+z}$")
scale_fact(c,1)
# 1.380743453918788e10

tgyr = range(0.1,30.8,length=10_000)
dgf(t) = scale_fact(c,t)
a_vals = dgf.(tgyr)
plot(tgyr,a_vals,lw=2.5, label="")
dgf_part(t) = scalefact_1by1(c,t)
a_vals_part = dgf_part.(tgyr)
plot!(tgyr,a_vals_part,lw=2.5, label="")
H0 = 67.74#(67.74/(3.26*9.5))*1e-18
Ωm = 0.32
ΩΛ = 0.68
Ωγ = 5.35e-5
Ων = 3.64e-5
Ωr = Ωγ + Ων

scalefact_part(c,100)
# ((67.74/(3.26*9.5))*1e-18)*(3600*24*365)/1e-9
# Mpc_to_km(mp_km) = mp_km * 3.262e6 * 9.461e12