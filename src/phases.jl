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

abstract type AbstractPhase end



@doc raw"""
    PolynomialPhase(coefs::Vector)  # polynomial coefficients

Return a `PolynomialPhase` object representing the phase defined by 
```math
g(z) = \sum_{j=0}^J \alpha_j z^j,
```
where `\alpha_j` are the coefficients in `coefs` and `J` is the degree of the polynomial.

It is assumed that `\alpha_J \neq 0` and `J \geq 2`. For linear phases, use `LinearPhase` instead.
"""
struct PolynomialPhase{T} <: AbstractPhase # arbitrary polynomial
    p   :: Polynomial{T} # coefficients of the polynomial
    dp  :: Polynomial{T} # coefficients of the derivative of polynomial
    dp2 :: Polynomial{T} # coefficients of the second derivative of polynomial
    ξ   :: Vector{ComplexF64} # stationary points
    v   :: Vector{Float64} # angles of each valley
    rstar_valley :: Float64
    function PolynomialPhase(coefs::Vector{T}) where T
        @assert coefs[end] != 0 "Leading coefficient must be non-zero"
        @assert coefs != [0,1]  "Use LinearPhase instead"
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

degree(G::PolynomialPhase)                   = length(G.p)-1
stationary_points(G::PolynomialPhase)        = G.ξ
evalphase(z, G::PolynomialPhase)             = G.p(z)
evalphase_derivative(z, G::PolynomialPhase)  = G.dp(z)
evalphase_derivative2(z, G::PolynomialPhase) = G.dp2(z)
rstar_valley(G::PolynomialPhase) = G.rstar_valley
infvalleys(G::PolynomialPhase) = G.v

function rStar(p::Polynomial)
    # define threshold distance for valley region
    α = coeffs(p)
    J = length(α)-1
    β = [k*abs(α[k+1]) for k = 1:J-1]
    poly = Polynomial([β; -J*abs(α[J+1])/sqrt(2)]) 
    rstar = maximum(real.(roots(poly))) # solution is the only positive root
    return rstar + 1e-6 # PATCH FIX - rstar = 0 for monomials
end

function evaluate_noreturn_Ginf(r,θ,G::PolynomialPhase)
    J = degree(G)
    αj = coeffs(G.p)
    J*αj[end]*r^(J-1) * min(1/sqrt(2), cos(J*θ)) - sum([j*αj[j+1]*r^(j-1) for j=1:J-1])
end


@doc raw"""
    LinearPhase()  

Return a `LinearPhase` object representing the phase defined by 
```math
g(z) = z,
```
"""
struct LinearPhase <: AbstractPhase 
    ξ :: Vector{Float64} # using Vector struct for consistency
    function LinearPhase()
        new([])
    end
end

stationary_points(G::LinearPhase)       = G.ξ
evalphase(z, ::LinearPhase)             = z
evalphase_derivative(z, ::LinearPhase)  = 1.0
evalphase_derivative2(z,::LinearPhase)  = 0.0


@doc raw"""
    SquareRootPhase(a, b)

Return a `SquareRootPhase` object representing the phase defined by 
```math
g(z) = \sqrt{z^2+a^2} + bz
```
where `a>0` and `b` is in `[-1,1]`.
"""
struct SquareRootPhase{T} <: AbstractPhase
    a :: T
    b :: T
    ξ :: Vector{ComplexF64} # using Vector struct for consistency
    function SquareRootPhase(a, b)
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

stationary_points(G::SquareRootPhase)        = G.ξ
evalphase(z, G::SquareRootPhase)             = sqrt(z^2 + G.a^2) + G.b * z
evalphase_derivative(z, G::SquareRootPhase)  = z/sqrt(z^2 + G.a^2) + G.b
evalphase_derivative2(z, G::SquareRootPhase) = G.a^2 / (z^2 + G.a^2)^(3/2)

# we use knowledge of the phase when "a" is small
# parameter below affects tolerance to decide whether the 
# integral is singular, and how to choose contour deformation by hand.
singular_tol(::SquareRootPhase) = 2. # this is an arbitrary choice


@doc raw"""
    RationalPhase(cpoly::Vector, poles::Vector, cs::Vector)

Return a `RationalPhase` object representing the phase defined by 
```math
g(z) = \sum_{j=0}^J \alpha_j z^j + \sum_{p=1}^P \sum_{k=1}^{K_p} \frac{\alpha_{p,k}}{(z-z_p)^k},
```
where `\alpha_j` are the coefficients in `cpoly` and `J` is the degree of the polynomial part,
`z_p` are the poles in `poles`, `\alpha_{p,k}` are the coefficients in `cs[p]` and `K_p` is the order of the pole `z_p`.

It is assumed that the singular part is non-zero (i.e. `cs` is not a vector of zero vectors).
"""
struct RationalPhase <: AbstractPhase 
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
    function RationalPhase(analytic_coefs::Vector,poles::Vector, poles_coefs::Vector)
        
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

stationary_points(G::RationalPhase)        = G.ξ
poles(G::RationalPhase)                    = G.p
evalphase(z, G::RationalPhase)             = G.rat(z)
evalphase_derivative(z, G::RationalPhase)  = G.drat(z)
evalphase_derivative2(z, G::RationalPhase) = G.ddrat(z)

degree(G::RationalPhase) = length(G.coefs_analytic)-1

rstar_valley(G::RationalPhase) = G.rstar_valley
rstar_pole(G::RationalPhase) = G.rstar_pole
infvalleys(G::RationalPhase) = G.vinf
polevalleys(G::RationalPhase) = G.vpole


Base.show(io::IO, G::RationalPhase) = begin
        analytic = Polynomial(G.coefs_analytic)
        singular = G.rat - analytic
        print(io, "RationalPhase: $analytic + $singular")
end

function compute_rstar_valley(αj, poles, αpk)
    J  = length(αj)-1
    

    # Get regularising term ∏(r-|zp|)^(Kp+1)
    # reg = 1.0
    # for (p,zp) in enumerate(poles)
    #     Kp = length(αpk[p])
    #     reg *= fromroots(ones(Kp+1) * abs(zp)) # (r-|zp|)^(Kp+1)
    # end

    # Analytic part: J*|αJ|*r^(J-1)/√2 - ∑j*|αj|*r^(j-1)
    analytic_part = Polynomial([[-j*abs(αj[j+1]) for j = 1:J-1]; J*abs(αj[end])/sqrt(2)])
    # analytic_part *= reg

    # Singular part: ∑∑ |k*α_(p,k)| * (z-zp)^(k-1)
    singular_part = 0.0
    for (p,zp) in enumerate(poles)
        Kp = length(αpk[p])
        for k = 1:Kp
            singular_part += k * abs(αpk[p][k]) * (1.0 // fromroots(abs(zp) * ones(k+1))) 
        end
    end
    # @show singular_part
    # singular_part = lowest_terms(singular_part * reg) # cancel out singularities
    # @assert singular_part.den(1.0) ≈ 1.0 # if this fails there is a bug

    # G = analytic_part * reg - singular_part_regularised.num

    r = roots(analytic_part - singular_part)[1] # works with Polynomial.RationalFunction
    rstar = maximum(real.(r))
    M = (maximum(abs.(poles)))
    # return rstar > M ? rstar : M
    # @assert rstar > M "Simplied criteria for valleys at infinity not verified \n rstar = $rstar < M = $M"
    return maximum(rstar)

    # rstar = maximum(real.(roots(G))) + 1e-12 # fix for monomial analytic part
    # @assert rstar > (maximum(abs.(poles)))
    # return rstar
end

function evaluate_noreturn_Ginf(r,θ,G::RationalPhase)
    αj = G.coefs_analytic
    αpk = G.coefs_singular
    J = degree(G)

    S = zero(Float64)
    S += J*abs(αj[end])*min(1/sqrt(2), cos(J*θ))*r^(J-1)
    S -= sum([j*abs(αj[j+1]*r^(j-1)) for j=1:J-1])
    for (p,zp) in enumerate(poles(G))
        Kp = length(αpk[p])
        S -= sum([k*abs(αpk[p][k])*abs(r-abs(zp))^(-k-1) for k=1:Kp])
    end
    return S
end

function compute_rstar_pole(αj, poles, αpk)
    J  = length(αj)-1
    lenP = length(poles)

    rp = zeros(lenP)
    for (p,zp) in enumerate(poles)
        Kp = length(αpk[p])

        # Dominant term
        singular_part = Kp*abs(αpk[p][end])/sqrt(2) * 1.0//fromroots(zeros(Kp+1)) # * r^(-Kp-1)
        # Lower order terms associated with zp
        if Kp > 1
            singular_part -= sum([k*abs(αpk[p][k]) * 1.0//fromroots(zeros(k+1)) for k=1:Kp-1])
        end
        # Terms associated to other poles
        for pp in setdiff(1:lenP, p)
            Kpp = length(αpk[pp])
            singular_part -= sum([k*abs(αpk[pp][k]) * 1.0//fromroots(ones(k+1)*abs(zp-poles[pp])) for k=1:Kpp])
        end
        # Polynomial part terms
        analytic_part = sum([j*abs(αj[j+1]) * fromroots(-abs(zp)*ones(j-1)) for j=1:J])

        r = roots(singular_part - analytic_part)[1] # works with Polynomial.RationalFunction
        rstar = minimum(real.(r[abs.(imag.(r)) .< 1e-6])) # filter imaginary roots

        if lenP == 1 Mp = Inf
        else Mp = minimum([abs(zp-poles[pp]) for pp in setdiff(1:lenP, p)])
        end
        @assert rstar < Mp "Simplied criteria for valleys at poles not verified"
        rp[p] =  maximum(real.(r))

    end
    return rp
end

function evaluate_noreturn_Gpole(r,θ,G::RationalPhase; pole_idx)
    αj = G.coefs_analytic
    αpk = G.coefs_singular

    p = pole_idx
    zp = poles(G)[p]
    Kp = length(αpk[p])
    J = degree(G)

    Gval = zero(Float64)
    Gval += Kp*abs(αpk[p][end]) * r^(-Kp-1) * min(1/sqrt(2), cos(Kp*θ))
    if Kp>1 Gval -= sum([k*abs(αpk[p][k]) * r^(-k-1) for k = 1:Kp-1]) end
    for pp in setdiff(1:length(poles(G)), p)
        Kpp = length(αpk[pp])
        zpp = poles(G)[pp]
        Gval -= sum([k*abs(αpk[pp][k]) * (r - abs(zp-zpp))^(-k-1) for k = 1:Kpp])
    end
    Gval -= sum([j*abs(αj[j+1]) * (r+abs(zp))^(j-1) for j = 1:J])
    return Gval
end