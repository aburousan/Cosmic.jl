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
    ],
)

deploydocs(;
    repo="github.com/aburousan/Cosmic.jl",
    devbranch="main",
)
