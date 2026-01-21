""" 
    Abstract type for phase functions

    These must have methods:
        - stationary_points()
        - evalphase(G,z) : evaluate phase function at z
        - evalphase_derivative(G,z) : evaluate first derivative of phase function at z
        - evalphase_derivative2(G,z) : evaluate second derivative of phase function at z

    These should have attributes:
        - ξ : stationary points of the phase function
        - other parameters as needed

    To incorporate a new phase function one should also modify the following functions
        - isinValley() in Contours.jl
        - find_zeros_range() in Balls.jl
        - exitpoints() in Balls.jl
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
        @assert coefs != [0,1]  "Use LinearPhaseFunction instead"
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

degree(G::PolynomialPhaseFunction)                   = length(G.p)-1
stationary_points(G::PolynomialPhaseFunction)        = G.ξ
evalphase(G::PolynomialPhaseFunction, z)             = G.p(z)
evalphase_derivative(G::PolynomialPhaseFunction, z)  = G.dp(z)
evalphase_derivative2(G::PolynomialPhaseFunction, z) = G.dp2(z)

function rStar(p::Polynomial)
    # define threshold distance for valley region
    α = coeffs(p)
    J = length(α)-1
    β = [k*abs(α[k+1]) for k = 1:J-1]
    poly = Polynomial([β; -J*abs(α[J+1])/sqrt(2)]) 
    rstar = maximum(real.(roots(poly))) # solution is the only positive root
    return rstar + 1e-6 # PATCH FIX - rstar = 0 for monomials
end

"""
    Struct for Linear phase function
"""

struct LinearPhaseFunction <: AbstractPhaseFunction 
    ξ :: Vector{Float64} # using Vector struct for consistency
    function LinearPhaseFunction()
        new([])
    end
end

stationary_points(G::LinearPhaseFunction)       = G.ξ
evalphase(::LinearPhaseFunction, z)             = z
evalphase_derivative(::LinearPhaseFunction, z)  = 1.0
evalphase_derivative2(::LinearPhaseFunction, z) = 0.0

"""
    Struct for rational phase function
"""

struct RationalPhaseFunction <: AbstractPhaseFunction 
    analytic  :: Polynomial # analytic part of the phase
    principal :: RationalFunction # singular part of the phase
    ξ :: Vector # stationary points
    p :: Vector # poles
    v :: Vector # valleys at infinity
    rat   :: RationalFunction # phase in the form p(z)/q(z)
    drat  :: RationalFunction # first derivative of phase
    ddrat :: RationalFunction # second derivative of phase
    rstar_valley :: Float64
    rstar_pole :: Float64
    function RationalPhaseFunction(analyticpart_coefs::Vector,poles::Vector)
        analytic_part = Polynomial(analyticpart_coefs)
        # den = fromroots(pole_vals) # Polynomial(den_coefs)

        p = unique(poles)
        poles_mult = [(zp, count(==(zp), poles)) for zp in p] # count multiplicities
        id = Polynomial(1.0)
        principal_part = RationalFunction(Polynomial(0.0), id)
        for (zp,mult) in poles_mult
            p = ones(length(mult)) * zp
            principal_part += id // fromroots(p)
        end
        rat = analytic_part + principal_part

        drat  = derivative(rat)
        ddrat = derivative(drat)

        dnum = derivative(rat.num)*rat.den - rat.num*derivative(rat.den)
        ξ = roots(dnum)  
        J   = length(analyticpart_coefs)-1
        v   = [((2*(m-1)+1/2)*π - angle(analyticpart_coefs[end]))/J for m=1:J]

        # compute r⋆ for poles and valleys only once and store the value
        rvalley = rStar_valley()
        rpole   = rStar_pole()
        new(analytic_part, principal_part, ξ ,p, v, rat, drat, ddrat, rvalley, rpole)
    end
end

# fix: degree of rat phase should be without the singular part??

# degree(G::PolynomialPhaseFunction)               = length(G.p)-1
stationary_points(G::RationalPhaseFunction)        = G.ξ
evalphase(G::RationalPhaseFunction, z)             = G.rat(z)
evalphase_derivative(G::RationalPhaseFunction, z)  = G.drat(z)
evalphase_derivative2(G::RationalPhaseFunction, z) = G.ddrat(z)
numerator(G::RationalPhaseFunction)   = G.rat.num
denominator(G::RationalPhaseFunction) = G.rat.den


function rStar_valley()
    return 2.0
    # define threshold distance for valley region
    # α = coeffs(p)
    # J = length(α)-1
    # β = [k*abs(α[k+1]) for k = 1:J-1]
    # poly = Polynomial([β; -J*abs(α[J+1])/sqrt(2)]) 
    # rstar = maximum(real.(roots(poly))) # solution is the only positive root
    # return rstar + 1e-6 # PATH FIX - rstar = 0 for monomials
end

function rStar_pole()
    return 0.1
end

""" 
    Struct for g(z) = √(z^2+a^2) + b*z
"""

struct SquareRootPhaseFunction{T} <: AbstractPhaseFunction
    a :: Float64
    b :: Float64 
    ξ :: Vector{ComplexF64} # using Vector struct for consistency
    function SquareRootPhaseFunction(a::T, b::T) where T
        @assert a>0 "Parameter `a` should be positive"
        @assert abs(b) ≤ 1 "Parameter `b` should be between -1 and 1"    
        # binf = (1-inftol) * sign(b) # ξ goes to ∞ as b tends to ±1
        # ξ = 1-abs(b) > inftol ? -a*b/sqrt(1-b^2) : -a*binf/sqrt(1-binf^2)
        ξ =  abs(b) == 1 ? -Inf*sign(b) : -a*b/sqrt(1-b^2)
        new{T}(a,b,[ξ])
    end
end

stationary_points(G::SquareRootPhaseFunction)        = G.ξ
evalphase(G::SquareRootPhaseFunction, z)             = sqrt(z^2 + G.a^2) + G.b * z
evalphase_derivative(G::SquareRootPhaseFunction, z)  = z/sqrt(z^2 + G.a^2) + G.b
evalphase_derivative2(G::SquareRootPhaseFunction, z) = G.a^2 / (z^2 + G.a^2)^(3/2)
