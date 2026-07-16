using Cosmic
using Documenter

DocMeta.setdocmeta!(Cosmic, :DocTestSetup, :(using Cosmic); recursive=true)

makedocs(;
    modules=[Cosmic],
    authors="Kazi Abu Rousan",
    repo="https://github.com/aburousan/Cosmic.jl/blob/{commit}{path}#{line}",
    sitename="Cosmic.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aburousan.github.io/Cosmic.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Background & species" => "background.md",
        "Neutrinos" => "neutrinos.md",
        "Thermal history" => "thermal_history.md",
        "BBN" => "bbn.md",
        "Perturbations" => Any[
            "Overview" => "perturbations.md",
            "Equations & derivations" => "perturbation_equations.md",
        ],
        "CMB, lensing & tensors" => "cmb.md",
        "Matter power spectrum" => "matter_power.md",
        "Primordial & inflation" => "primordial.md",
        "Induced gravitational waves" => "sigw.md",
        "Spectral distortions" => "spectral_distortions.md",
        "API index" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/aburousan/Cosmic.jl",
    devbranch="main",
)
