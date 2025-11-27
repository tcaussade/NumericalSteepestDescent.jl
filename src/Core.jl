""" 
    Abstract type for phase functions

    These must have methods:
        - evalphase(G,z) : evaluate phase function at z
        - evalphase_derivative(G,z) : evaluate derivative of phase function at z
        - evalinverse(G, η, s) : evaluate inverse of phase function at level set

    These must have attributes:
        - ξ : stationary points of the phase function
        - other parameters as needed
"""

abstract type AbstractPhaseFunction end


""" 
    Struct for complex contour parametrisations

    Symbol can be :finite, :infiniteSD or :finiteSD
"""

struct ComplexContour
    parametrisation :: Function
    type :: Symbol # Finite / Infinite SD / Finite SD
    orientation :: Symbol # 
    function ComplexContour(parametrisation::Function, type::Symbol, orientation::Symbol)
        @assert type in (:finite, :infiniteSD, :finiteSD)
        @assert orientation in (:positive, :negative)
        new(parametrisation, type, orientation)
    end
end

contour_type(γ::ComplexContour) = γ.type

_check_contour_type(γ::ComplexContour, s::Symbol ) = @assert γ.type == s "Contour type mismatch: expected $s, got $(γ.type)"

contour_orientation(γ::ComplexContour) = γ.orientation

function trace_finite(a,b)
    # parametrisation of finite straight line from a to b
    u -> 0.5*((b+a) + (b-a)*u) # :: Function
end

# function trace_infiniteSDpath(G::AbstractPhaseFunction, η)
#     # parametrisation of infinite SD path at η
#     g(z) = evalphase(G, z)
#     u -> evalinverse(G, η, g(η) + im*u) # :: Function
# end

"""
    Struct for polynomial phase functions
"""

struct PolynomialPhaseFunction{T} <: AbstractPhaseFunction # arbitrary polynomial
    p   :: Polynomial{T} # coefficients of the polynomial
    dp  :: Polynomial{T} # coefficients of the derivative of polynomial
    dp2 :: Polynomial{T} # coefficients of the second derivative of polynomial
    ξ   :: Vector{ComplexF64} # stationary points
    v   :: Vector{Float64} # angles of each valley
    function PolynomialPhaseFunction(coefs::Vector{T}) where T
        p   = Polynomial(coefs)
        dp  = derivative(p, 1)
        ξ   = roots(dp)
        dp2 = derivative(dp, 1)
        J   = length(coefs)-1
        v   = [((2*(m-1)+1/2)*π - angle(coefs[end]))/J for m=1:J]
        new{T}(p,dp,dp2,ξ,v)
    end
end

degree(G::PolynomialPhaseFunction) = length(G.p)-1

stationary_points(G::PolynomialPhaseFunction) = G.ξ

evalphase(G::PolynomialPhaseFunction, z) = G.p(z)