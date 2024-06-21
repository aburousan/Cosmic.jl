using Cosmic
using Test
using Plots

@testset "Cosmic.jl" begin
    @test redshift(2.9e-4) ≈ 3447.27586
    @test cosmology().h ≈ 0.6774
    @test scalefact_part(cosmology(),100) ≈ 148.9412838881
    @test ageGyr(cosmology(),0) ≈ 13.80743453918
end



query = """
SELECT source_id,
    phot_g_mean_mag,
    bp_rp,
    phot_g_mean_flux,
    parallax
FROM gaiadr3.gaia_source
WHERE bp_rp IS NOT NULL AND phot_g_mean_mag IS NOT NULL AND parallax > 50
ORDER BY phot_g_mean_mag ASC
"""

gaia_data = fetch_gaia_data(query)
data = filter_missing_values(gaia_data)

parallax = gaia_data[:, "parallax"]
distance = 1000.0 ./ parallax # distance in parsecs
abs_mag = gaia_data[:, "phot_g_mean_mag"] .+ 5 .* log10.(parallax / 100.0)

luminosity = 10 .^ ((4.83 .- abs_mag) ./ 2.5)

bp_rp = gaia_data[:, "bp_rp"]

normalized_bp_rp = (bp_rp .- minimum(bp_rp)) ./ (maximum(bp_rp) .- minimum(bp_rp))

cmap = cgrad(:plasma)

scatter(bp_rp, log10.(luminosity), 
        marker_z=normalized_bp_rp, 
        color=cmap, 
        xlabel="BP-RP", 
        ylabel="Log(Luminosity/Lsun)", 
        title="HR Diagram", 
        legend=false, 
        reverse=true)

ylims!((minimum(log10.(luminosity)) - 0.5, maximum(log10.(luminosity)) + 0.5))

plot!()