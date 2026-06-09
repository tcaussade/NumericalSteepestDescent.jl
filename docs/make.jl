using Documenter
using NumericalSteepestDescent

# DocMeta.setdocmeta!(NumericalSteepestDescent, :DocTestSetup, :(using LinearAlgebra, SpecialFunctions))

makedocs(;
    modules = [NumericalSteepestDescent],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://tcaussade.github.io/NumericalSteepestDescent.jl/",
        repolink = "https://github.com/tcaussade/NumericalSteepestDescent.jl",
    ),
    
    repo = "https://github.com/tcaussade/NumericalSteepestDescent.jl/blob/{commit}{path}#L{line}",
    sitename = "NumericalSteepestDescent.jl",
    authors = "Thomas Caussade",
    pages = [
        "Home" => "index.md",
        "Usage" => "use.md",
        "Phases" => "phases.md",
    ],
)

deploydocs(
    repo   = "github.com/tcaussade/NumericalSteepestDescent.jl.git",
)