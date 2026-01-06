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
    Struct for polynomial phase functions
"""

struct PolynomialPhaseFunction{T} <: AbstractPhaseFunction # arbitrary polynomial
    p   :: Polynomial{T} # coefficients of the polynomial
    dp  :: Polynomial{T} # coefficients of the derivative of polynomial
    dp2 :: Polynomial{T} # coefficients of the second derivative of polynomial
    ξ   :: Vector{ComplexF64} # stationary points
    v   :: Vector{Float64} # angles of each valley
    rStar :: Float64
    function PolynomialPhaseFunction(coefs::Vector{T}) where T
        @assert coefs[end] != 0 "Leading coefficient must be non-zero"
        p   = Polynomial(coefs)
        dp  = derivative(p, 1)
        ξ   = roots(dp)
        dp2 = derivative(dp, 1)
        J   = length(coefs)-1
        v   = [((2*(m-1)+1/2)*π - angle(coefs[end]))/J for m=1:J]
        r   = rStar(p)
        new{T}(p,dp,dp2,ξ,v,r)
    end
end

degree(G::PolynomialPhaseFunction) = length(G.p)-1
stationary_points(G::PolynomialPhaseFunction) = G.ξ
evalphase(G::PolynomialPhaseFunction, z) = G.p(z)
evalphase_derivative(G::PolynomialPhaseFunction, z) = G.dp(z)

function rStar(p::Polynomial)
    # define threshold distance for valley region
    α = coeffs(p)
    J = length(α)-1
    β = [k*abs(α[k+1]) for k = 1:J-1]
    poly = Polynomial([β; -J*abs(α[J+1])/sqrt(2)]) 
    return maximum(real.(roots(poly))) # solution is the only positive root
end


""" Struct for g(z) = √(z^2+a^2)
"""

struct DistancePhaseFunction{T} <: AbstractPhaseFunction
end

evalphase(G::DistancePhaseFunction, z) = sqrt(z^2 + G.a^2)

