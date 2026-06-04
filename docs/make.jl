using Documenter, NumericalSteepestDescent

makedocs(;
    modules = [NumericalSteepestDescent],
    format = Documenter.HTML(
        canonical = "https://tcaussade.github.io/NumericalSteepestDescent.jl/",
        repolink = "https://github.com/tcaussade/NumericalSteepestDescent.jl"
    ),
    pages = [
        "Home" => "index.md",
        # "Gaussian Quadrature" => "gaussquadrature.md",
        # "Benchmark" => "benchmark.md",
        # "Roots of Bessel function" => "besselroots.md",
        # "Misc." => "misc.md",
        # "References" => "reference.md",
    ],
    repo = "https://github.com/tcaussade/NumericalSteepestDescent.jl/blob/{commit}{path}#L{line}",
    sitename = "NumericalSteepestDescent.jl",
    authors = "Thomas Caussade",
)

deploydocs(; repo = "github.com/tcaussade/NumericalSteepestDescent.jl", devbranch = "main")