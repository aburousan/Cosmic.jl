using cosmic
using Test

@testset "cosmic.jl" begin
    @test my_f(2,1) == 11
    @test my_g(2) == 4
    @test redshift(2.9e-4) ≈ 3447.27586
end

using Plots, LaTeXStrings

c = cosmology()
z_vals = range(0, 2.5, length=1000)
H_1z(z) = H(c,z)/(1+z)
H_vals = H_1z.(z_vals)
plot(z_vals,H_vals,lw=2.5, label="")
xlabel!(L"Redshift ($z$)")
ylabel!(L"$\frac{H(z)}{1+z}$")