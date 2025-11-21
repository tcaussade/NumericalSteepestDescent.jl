

"""
    Core functions for particular phase function oscillatory quadrature

    The functions are specialised to special cases, where most computation can be 
    carreid out by hand a priori.
"""


"""
    Let g(z) = z^J where J>1 is a real number (not necessarily integer).  
"""

struct MonomialPhaseFunction <: AbstractPhaseFunction # g(z) = z^J
    J :: Float64 # power in monomial phase function
    ξ  :: Vector{ComplexF64} # stationary point(s)
    function MonomialPhaseFunction(J::Real)
        @assert J > 1.0 "Power J must be greater than 1"
        ξ = [0.0]
        new(J, ξ)
    end
end

# MonomialPhaseFunction(J::Real) = MonomialPhaseFunction(Float64(J))

# local_upperbound(::MonomialPhaseFunction) = 1.0
# local_lowerbound(::MonomialPhaseFunction) = 1.0

function evalphase(G::MonomialPhaseFunction, z)
    return z^G.J
end
function evalphase_derivative(G::MonomialPhaseFunction, z)
    return G.J * z^(G.J - 1)
end

function evalinverse(G::MonomialPhaseFunction, η, s)
    # this might need re working for general case...
    return s^(1/G.J)
end

function ballradius(G::MonomialPhaseFunction, ::Number, C)
    @info "Using ball radius computed a priori"
    return C^(G.J) # C = Cball / ω
end

function exitpoints(::MonomialPhaseFunction, Ω)
    c,r = centre_and_radius(Ω[1])
    return [c-r, c+r]
end


function QuasiSDcontour(G::MonomialPhaseFunction, NonOscBall, L) :: Vector{ComplexContour}
    # compute quasi-sd contour deformation for monomial phase function
    # we assume [0,L] with L>1 is the integration domain

    ξ,r = NonOscBall

    if abs(ξ-0.0)<r && abs(ξ-L)<r # no deformation required
        return [ComplexContour(trace_finite(0.0, L) , :finite, :positive)]
    end

    ηc = ξ + r # exit point from non-oscillatory region
    
    f0 = ComplexContour(trace_finite(0.0, ηc), :finite, :positive)
    hc = ComplexContour(trace_infiniteSDpath(G, ηc), :infiniteSD, :positive)
    hL = ComplexContour(trace_infiniteSDpath(G, L), :infiniteSD, :negative)
    γ = [f0, hc, hL]
    return γ
end
