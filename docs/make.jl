using Documenter
using NumericalSteepestDescent

# Setup for doctests in docstrings
DocMeta.setdocmeta!(NumericalSteepestDescent)

makedocs(;
    modules = [NumericalSteepestDescent],
    format = Documenter.HTML(
        canonical = "https://tcaussade.github.io/NumericalSteepestDescent/",
        assets = [""],
        repolink = "https://github.com/tcaussade/NumericalSteepestDescent.jl",
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    pages = [
        "Home" => "intro.md",
        "References" => "bib.md",
    ],
    repo = "https://github.com/tcaussade/NumericalSteepestDescent.jl",
    sitename = "NumericalSteepestDescent.jl",
    authors = "Me",
)

# Only deploy from GitHub Actions CI
if get(ENV, "CI", "false") == "true" && get(ENV, "GITHUB_EVENT_NAME", "") == "push"
    deploydocs(; repo = "https://github.com/tcaussade/NumericalSteepestDescent.jl")
end
