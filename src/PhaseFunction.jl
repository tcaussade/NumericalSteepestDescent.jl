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
    rstar_valley :: Float64
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
evalphase(z, G::PolynomialPhaseFunction)             = G.p(z)
evalphase_derivative(z, G::PolynomialPhaseFunction)  = G.dp(z)
evalphase_derivative2(z, G::PolynomialPhaseFunction) = G.dp2(z)
rstar_valley(G::PolynomialPhaseFunction) = G.rstar_valley
infvalleys(G::PolynomialPhaseFunction) = G.v

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
evalphase(z, ::LinearPhaseFunction)             = z
evalphase_derivative(z, ::LinearPhaseFunction)  = 1.0
evalphase_derivative2(z,::LinearPhaseFunction)  = 0.0

""" 
    Struct for g(z) = √(z^2+a^2) + b*z
"""

struct SquareRootPhaseFunction{T} <: AbstractPhaseFunction
    a :: T
    b :: T
    ξ :: Vector{ComplexF64} # using Vector struct for consistency
    function SquareRootPhaseFunction(a, b)
        @assert a>0 "Parameter `a` should be positive"
        @assert abs(b) < 1 "Parameter `b` should be in (-1,1)"   
        # @info "Methods assume the integration contour is real" 
        # binf = (1-inftol) * sign(b) # ξ goes to ∞ as b tends to ±1
        # ξ = 1-abs(b) > inftol ? -a*b/sqrt(1-b^2) : -a*binf/sqrt(1-binf^2)
        ξ =  abs(b) == 1 ? -Inf*sign(b) : -a*b/sqrt(1-b^2)
        T = promote_type(typeof(a), typeof(b))
        new{T}(a,b,[ξ])
    end
end

stationary_points(G::SquareRootPhaseFunction)        = G.ξ
evalphase(z, G::SquareRootPhaseFunction)             = sqrt(z^2 + G.a^2) + G.b * z
evalphase_derivative(z, G::SquareRootPhaseFunction)  = z/sqrt(z^2 + G.a^2) + G.b
evalphase_derivative2(z, G::SquareRootPhaseFunction) = G.a^2 / (z^2 + G.a^2)^(3/2)

# we use knowledge of the phase when "a" is small
# parameter below affects tolerance to decide whether the 
# integral is singular, and how to choose contour deformation by hand.
singular_tol(::SquareRootPhaseFunction) = 2. # this is an arbitrary choice

"""
    Struct for rational phase function
"""

struct RationalPhaseFunction <: AbstractPhaseFunction 
    rat   :: RationalFunction # phase in the form p(z)/q(z)
    drat  :: RationalFunction # first derivative of phase
    ddrat :: RationalFunction # second derivative of phase
    ξ :: Vector # stationary points
    p :: Vector{ComplexF64} # poles
    vinf  :: Vector{ComplexF64} # valleys at infinity
    rstar_valley :: Float64
    vpole :: Vector{Vector{ComplexF64}} # valleys at each pole
    rstar_pole :: Vector{Float64}

    coefs_analytic :: Vector
    coefs_singular :: Vector
    function RationalPhaseFunction(analytic_coefs::Vector,poles::Vector, poles_coefs::Vector)
        
        @assert length(poles) == length(poles_coefs) "Coefficients of singular part are not well specified"
        @assert ~iszero(poles_coefs) "Phase has null singular part"
        
        # construct the function
        id = Polynomial(1.0)
        singular_part = Polynomial(0.0)
        for (i,zp) in enumerate(poles)
            for (k,coef) in enumerate(poles_coefs[i])
                pvec = +zp * ones(k) 
                singular_part += coef * id // fromroots(pvec) 
            end
        end
        rat = lowest_terms(Polynomial(analytic_coefs)+ singular_part)
        drat  = derivative(rat)
        ddrat = derivative(drat)

        dnum = derivative(rat.num)*rat.den - rat.num*derivative(rat.den)
        ξ = Vector{ComplexF64}(undef,0)
        for root in roots(dnum)
            minimum(abs.(root.-poles)) > 1e-2 ? push!(ξ,root) : continue
        end
        
        # setdiff(roots(dnum), poles) # dnum may have more solutions than we need  
        J = length(analytic_coefs)-1

        # valleys at infinity
        vinf   = [((2*(m-1)+1/2)*π - angle(analytic_coefs[end]))/J for m=1:J]
        rinf = compute_rstar_valley(analytic_coefs, poles, poles_coefs)

        # valleys at poles
        vpole = Vector{Vector{ComplexF64}}(undef, length(poles))
        for p = 1:length(poles)
            Kp = length(poles_coefs[p])
            vpole[p] = [(-(2*(m-1)+1/2)*π + angle(poles_coefs[p][end]))/Kp for m=1:Kp]
        end
        rpole   = compute_rstar_pole(analytic_coefs, poles, poles_coefs)

        new(rat, drat, ddrat, ξ ,poles, 
            vinf, rinf, vpole, rpole,
            analytic_coefs, poles_coefs)
    end
end

stationary_points(G::RationalPhaseFunction)        = G.ξ
poles(G::RationalPhaseFunction)                    = G.p
evalphase(z, G::RationalPhaseFunction)             = G.rat(z)
evalphase_derivative(z, G::RationalPhaseFunction)  = G.drat(z)
evalphase_derivative2(z, G::RationalPhaseFunction) = G.ddrat(z)

degree(G::RationalPhaseFunction) = length(G.coefs_analytic)-1

rstar_valley(G::RationalPhaseFunction) = G.rstar_valley
rstar_pole(G::RationalPhaseFunction) = G.rstar_pole
infvalleys(G::RationalPhaseFunction) = G.vinf
polevalleys(G::RationalPhaseFunction) = G.vpole


Base.show(io::IO, G::RationalPhaseFunction) = begin
        analytic = Polynomial(G.coefs_analytic)
        singular = G.rat - analytic
        print(io, "RationalPhase: $analytic + $singular")
end

function compute_rstar_valley(αj, poles, αpk)
    J  = length(αj)-1
    id = Polynomial(1.0)
    # Analytic part: J*|αJ|*r^(J-1)/√2 - ∑j*|αj|*r^(j-1)
    analytic_part = Polynomial([[-(j-1)*abs(αj[j]) for j = 2:J]; J*abs(αj[end])/sqrt(2)])

    # Singular part: ∑∑ |k*α_(p,k)| * (z-zp)^(k-1)
    singular_part = RationalFunction(Polynomial(0.0), id)
    reg = id  # used to convert into a polynomial equation

    id = Polynomial(1.0)
    singular_part = Polynomial(0.0)
    for (i,zp) in enumerate(poles)
        for (k,coef) in enumerate(αpk[i])
            pvec = abs(zp) * ones(k+1) 
            singular_part += k * abs(coef) * id // fromroots(pvec) 
        end
        mult = length(αpk[i])
        reg *= fromroots(ones(mult+1) * abs(zp))
    end

    singular_part_regularised = lowest_terms(reg * singular_part) # cancel out singularities
    @assert singular_part_regularised.den(1.0) ≈ 1.0 # if this fails there is a bug
    G = analytic_part * reg - singular_part_regularised.num

    return maximum(real.(roots(G))) + 1e-12 # fix for monomial analytic part
end

function compute_rstar_pole(αj, poles, αpk)
    J  = length(αj)-1
    id = Polynomial(1.0)

    rp = zeros(length(poles))
    for (i,zp) in enumerate(poles)
        Kp = length(αpk[i])

        # Analytic part: Kp*|α_(P,Kp)|*r^(-Kp-1) - ∑j*|αj|*(|zp|+r)^(j-1)
        analytic_part = Polynomial(0.0)
        for j = 1:J
            vec = -abs(zp) * ones(j-1)
            analytic_part += j * abs(αj[j+1]) * fromroots(vec)
        end
        # Singular part
        singular_part = Polynomial(0.0)

        singular_part += Kp * abs(αpk[i][end]) / sqrt(2) * id // fromroots(zeros(Kp+1)) 
        reg = fromroots(zeros(Kp+1)) # used to convert into polynomial equation

        for k=1:Kp-1
            singular_part -= k * abs(αpk[i][k]) * id // fromroots(zeros(k+1))
        end

        for pp in setdiff(1:length(poles), i) # iterate over p' ≠ p
            zpp = poles[pp] 
            Kpp = length(αpk[pp])
            # @show αpk[pp]
            for k = 1:Kpp # from k = 1 to K_p'
                ppvec = ones(k+1) * abs(zp - zpp)
                singular_part -= (-1)^(-k-1) * k * abs(αpk[pp][k]) * id // fromroots(ppvec)
                # singular_part -= k * abs(αpk[pp][k]) * id // fromroots(ppvec)
            end
            reg *= fromroots(ones(Kpp+1) * abs(zpp - zp)) #* (-1)^(Kpp+1)
        end
        """Testing"""
        # @show singular_part
        # @show analytic_part

        singular_part_regularised = lowest_terms(singular_part * reg)
        @assert singular_part_regularised.den(1.0) ≈ 1.0 # if this fails there is a bug
        G = singular_part_regularised.num - analytic_part * reg

        # @show roots(G)
        rp[i] = maximum(real.(roots(G))) + 1e-12 # fix for monomial analytic part
    end
    return rp
end

function evaluate_noreturn_Gpole(r,θ,G::RationalPhaseFunction; pole_idx)
    αj = G.coefs_analytic
    αpk = G.coefs_singular

    i = pole_idx
    zp = poles(G)[i]
    Kp = length(polevalleys(G)[i])
    J = degree(G)

    Gval = zero(Float64)
    Gval += Kp*abs(αpk[i][end]) * r^(-Kp-1) * min(1/sqrt(2), cos(Kp*θ))
    Gval -= sum([j*abs(αj[j+1]) * (abs(zp)+r)^(j-1) for j = 1:J])
    Gval -= sum([k*abs(αpk[i][k]) * r^(-k-1) for k = 1:Kp-1])
    for pp in setdiff(1:length(poles(G)), i)
        Kpp = length(αpk[pp])
        zpp = poles(G)[pp]
        Gval -= sum([k*abs(αpk[pp][k]) * (abs(zpp-zp) - r)^(-k-1) for k = 1:Kpp])
    end
    return Gval
end