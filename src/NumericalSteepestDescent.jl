module NumericalSteepestDescent

"""
    This package provides tools for numerical integration of highly oscillatory integrals 
    using contour deformation based on the (regularised) method of numerical steepest descent.

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

include("phases.jl")

include("root_finding.jl")
include("graph_traversal.jl")

include("balls.jl")
include("contours.jl")
include("contourgraph.jl")
include("quadrature.jl")
include("plotgeneration.jl")
include("residues.jl")

include("diagnosistools.jl")
include("api.jl")


export 
PolynomialPhase,
RationalPhase,
LinearPhase,
SquareRootPhase,
nsd,
plot_quasiSDdeformation,
quasiSDdeformation!,
quasiSDdeformation,
QuasiSDcontour,
NonOscillatoryRegion,
integrate

end # module NumericalStationaryPhase
