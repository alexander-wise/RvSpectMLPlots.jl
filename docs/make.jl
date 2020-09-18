using RvSpectMLPlots
using Documenter

makedocs(;
    modules=[RvSpectMLPlots],
    authors="Eric Ford",
    repo="https://github.com/RvSpectML/RvSpectMLPlots.jl/blob/{commit}{path}#L{line}",
    sitename="RvSpectMLPlots.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://RvSpectML.github.io/RvSpectMLPlots.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/RvSpectML/RvSpectMLPlots.jl",
)
