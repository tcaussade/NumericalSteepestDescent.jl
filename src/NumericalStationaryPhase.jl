module NumericalStationaryPhase

using FastGaussQuadrature
using Polynomials
using Roots
using LinearAlgebra

using CairoMakie
using CairoMakie.Colors
using Graphs
using GraphMakie


include("Core.jl")

include("nonoscillatoryregion.jl")

include("contourgraph.jl")
include("PlotGeneration.jl")
include("quadrature.jl")

# include("special_cases.jl")

include("api.jl")

export 
PolynomialPhaseFunction,
integrate_nsp,
plot_quasiSDdeformation,
QuasiSDcontour,
NonOscillatoryRegion
# specialcases
# MonomialPhaseFunction



end # module NumericalStationaryPhase
