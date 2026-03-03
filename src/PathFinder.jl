module PathFinder

"""
    This package provides tools for numerical integration of highly oscillatory integrals 
    using contour deformation based on the method of numerical stationary phase.

    For a MATLAB version, see https://github.com/AndrewGibbs/PathFinder
"""

using FastGaussQuadrature
using Polynomials
using Roots
using LinearAlgebra
using ForwardDiff
using Graphs
using QuadGK
# Plots
# using WGLMakie
# using WGLMakie.Colors
using CairoMakie
using CairoMakie.Colors
using GraphMakie

include("PhaseFunction.jl")

include("root_finding.jl")

include("Balls.jl")
include("Contours.jl")
include("ContourGraph.jl")
include("quadrature.jl")
include("PlotGeneration.jl")
include("Residues.jl")

include("api.jl")

export 
PolynomialPhaseFunction,
RationalPhaseFunction,
LinearPhaseFunction,
SquareRootPhaseFunction,
integrate_nsp,
plot_quasiSDdeformation,
QuasiSDcontour,
NonOscillatoryRegion,
integrate

end # module NumericalStationaryPhase
