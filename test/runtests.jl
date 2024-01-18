using cosmic
using Test

@testset "cosmic.jl" begin
    @test my_f(2,1) == 11
    @test my_g(2) == 4
    @test z(2.9e-4) ≈ 3447.27586
end

Tcmb = 2.7260
h = 0.740
OmegaG = 4.48131e-7 * Tcmb^4 / h^2
OmegaN = Neff * OmegaG * (7 / 8) * (4 / 11)^(4 / 3)
OmegaR = OmegaG + OmegaN

