using cosmic
using Documenter

DocMeta.setdocmeta!(cosmic, :DocTestSetup, :(using cosmic); recursive=true)

makedocs(;
    modules=[cosmic],
    authors="Kazi Abu Rousan",
    repo="https://github.com/aburousan/cosmic.jl/blob/{commit}{path}#{line}",
    sitename="cosmic.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://aburousan.github.io/cosmic.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/aburousan/cosmic.jl",
    devbranch="main",
)
