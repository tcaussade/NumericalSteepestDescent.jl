module PathFinder

using FastGaussQuadrature
using Polynomials
using Roots
using LinearAlgebra
using CairoMakie
using CairoMakie.Colors
using Graphs
using GraphMakie

include("PhaseFunction.jl")

include("root_finding.jl")

include("Balls.jl")
include("Contours.jl")
include("ContourGraph.jl")
include("quadrature.jl")
include("PlotGeneration.jl")

# include("special_cases.jl")

include("api.jl")

export 
PolynomialPhaseFunction,
integrate_nsp,
plot_quasiSDdeformation,
QuasiSDcontour,
NonOscillatoryRegion,
integrate
# specialcases
# MonomialPhaseFunction



end # module NumericalStationaryPhase
