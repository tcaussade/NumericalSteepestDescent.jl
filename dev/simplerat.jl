using NumericalSteepestDescent
using QuadGK
using FastGaussQuadrature
using CairoMakie


ω = 2.
a,b = -2, 2
γ0 = 

# trying g(z) = z + 1/z
acoefs = [0,1]
ps     = [0.0]
pcoefs = [[1]] #

RatPhase = RationalPhaseFunction(acoefs, ps, pcoefs)

Ω = NonOscillatoryRegion(RatPhase, ω; Cball = 2π, δball = 1e-3,  Nrays = 16)

CG, CtoG, NodesDict, EdgesList = NumericalSteepestDescent.ContourGraph(RatPhase, a, b, Ω; δODE = 0.1, δcoarse = 0.01)
# a,b = NodesDict[:endpoint]
sd_edges = NumericalSteepestDescent.get_deformation(RatPhase, CG, a,b, CtoG, NodesDict, EdgesList, γ0)

# choose quadrature and precompute weights and nodes
if quadtype == :gaussian
    qleg = gausslegendre(N)
    qlag = gausslaguerre(N)
    if G isa SquareRootPhaseFunction 
        qsin = gausslegendre(Int(floor(N/layersnumber(G.a))) + 1) 
    else qsin = nothing end
    quad = (qleg = qleg, qlag = qlag, qsingular = qsin)
elseif quadtype == :adaptive
else
    @error "quadtype is not in (:gaussian, :adaptive)"
end

# Evaluate each contour on the shortest path
S = zero(ComplexF64)
γtot = Vector{ComplexContour}()
for e in sd_edges 
# for i in eachindex(sd_nodes)[2:end]
    # i1 = sd_nodes[i-1]
    # i2 = sd_nodes[i]
    # i1,i2 = e.src, e.dst
    if haskey(EdgesList, e)
        γ = EdgesList[e]; push!(γtot, γ)
        if quadtype == :gaussian
            S += integrate(γ, f, G, ω, quad; δfine, δquad)
        elseif quadtype == :adaptive
            S += integrate(γ, f, G, ω; δfine, δquad, atol)
        end
    else # the contour is traversed in the opposite direction
        γ = EdgesList[reverse(e)]; push!(γtot, γ)
        if quadtype == :gaussian
            S -= integrate(γ, f, G, ω, quad; δfine, δquad)
        elseif quadtype == :adaptive
            S -= integrate(γ, f, G, ω; δfine, δquad, atol)
        end
    end
end