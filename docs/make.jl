using Documenter
using NumericalSteepestDescent

# Setup for doctests in docstrings
DocMeta.setdocmeta!(NumericalSteepestDescent, :DocTestSetup, :(using LinearAlgebra, SpecialFunctions, FastGaussQuadrature))

makedocs(;
    modules = [NumericalSteepestDescent],
    format = Documenter.HTML(
        canonical = "https://tcaussade.github.io/NumericalSteepestDescent.jl/",
        assets = [""],
        repolink = "https://github.com/tcaussade/NumericalSteepestDescent.jl"
    ),
    pages = [
        "Home" => "intro.md",
        "References" => "bib.md",
    ],
    repo = "https://github.com/tcaussade/NumericalSteepestDescent.jl",
    sitename = "NumericalSteepestDescent.jl",
    authors = "Me",
)

deploydocs(; repo = "https://github.com/tcaussade/NumericalSteepestDescent.jl")
