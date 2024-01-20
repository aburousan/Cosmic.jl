using Cosmic
using Test

@testset "Cosmic.jl" begin
    @test my_f(2,1) == 11
    @test my_g(2) == 4
    @test redshift(2.9e-4) ≈ 3447.27586
    @test cosmology().h ≈ 0.6774
    @test scalefact_part(cosmology(),100) ≈ 148.9412838881
    @test ageGyr(cosmology(),0) ≈ 13.80743453918
end