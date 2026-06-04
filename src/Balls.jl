
"""
    NonOscillatoryBall

Struct representing a non-oscillatory ball in the complex plane. 
"""

struct NonOscillatoryBall 
    c :: ComplexF64 # centre, e.g. a stationary point
    r :: Float64    # radius
end
centre_and_radius(B::NonOscillatoryBall) = (B.c, B.r)

ballradius(G::AbstractPhase, ξ, Cω; Nrays) = _compute_ballradius(G, ξ, Cω; Nrays)

function _compute_ballradius(G::AbstractPhase, ξ, Cω; Nrays)
    r = zeros(Nrays)
    for n = 1:Nrays
        θ = 2*n/Nrays
        r[n] = findradius(G,ξ,Cω,θ)
    end
    return minimum(r)
end

function findradius(G::AbstractPhase, ξ, Cω, θ)
    ### THIS IS A SLOW PART OF THE CODE!!!
    g(z) = evalphase(z,G)
    ray(r) = ξ + r*cispi(θ)
    un(r) = abs(g(ray(r)) - g(ξ))^2 - Cω^2
    guess = find_zeros_range(G,ξ) # heuristic choice - can we do better?
    # return find_zero(un, guess, Bisection())
    rs = find_zeros(un, guess[1], guess[2], no_pts = 41)
    if isempty(rs) return Inf else return minimum(rs) end
end

"""
    Construct the non-oscillatory region Ω

G      : phase function (subtype of AbstractPhaseFunction)
Cball  : constant defining the size of the non-oscillatory region
ω      : frequency parameter
δball  : overlap tolerance between non-oscillatory balls (default 1.0)

Returns a Vector of NonOscillatoryBalls representing the non-oscillatory region.
"""


function NonOscillatoryRegion(G::AbstractPhase, ω; 
                              Cball, δball, Nrays)
    # δball = 1e-3 / 2 / max(degree(G)-2,1)
    Ω = Vector{NonOscillatoryBall}()
    # Assign a NonOscillatoryBall to each stationary point
    for ξ in G.ξ
        r = ballradius(G, ξ, Cball/ω; Nrays)
        push!(Ω, NonOscillatoryBall(ξ, r)) 
    end

    Ω = unique!(Ω) # remove duplicate balls (if any)
    # Remove balls with significant overlap
    n = length(Ω)
    idx = Int[]
    for i = 1:n
        ci,ri = centre_and_radius(Ω[i])
        for j = i+1:n
            cj,rj = centre_and_radius(Ω[j])
            dist = abs(ci - cj)/max(ri, rj)
            if dist < δball # overlap detected, remove the smaller ball
                Ω[i].r ≤ Ω[j].r ? push!(idx, i) : push!(idx, j)
            end
        end
    end
    deleteat!(Ω, sort(unique!(idx)))
    return Ω
end

function get_Pstat(Ω::Vector{NonOscillatoryBall})
    # Get (relevant) stationary points in non-oscillatory region
    Pstat = ComplexF64[]
    for Ball in Ω
        c, _ = centre_and_radius(Ball)
        push!(Pstat, c)
    end
    return Pstat
end

dist(hη, Pstat :: Vector) = minimum(abs.(hη .- Pstat)) # dist(hn, P_statpoint)

"""
    Determine the set of exit points of the non-oscillatory region
"""

function exitpoints(G::AbstractPhase, Ω::Vector{NonOscillatoryBall})
    Pexit = ComplexF64[]
    g(z) = evalphase(z,G)
    for Ball in Ω
        c,r = centre_and_radius(Ball) 
        trig   = θ -> imag(g(c + r*cis(θ)))
        dtrig  = θ -> ForwardDiff.derivative(trig,θ)  # first derivative of Im(g)
        ddtrig = θ -> ForwardDiff.derivative(dtrig,θ) # second derivative of Im(g)
        θ = find_zeros(dtrig, 0,  2π) # find the roots of dtrig

        maxima = Float64[] # second-derivative test
        [ddtrig(θ) < 0.0 ? push!(maxima, θ) : nothing for θ in θ] 

        exits = [c+ r*cis(θ) for θ in maxima]
        [!isinΩ(setdiff(Ω, [Ball]), z) ? push!(Pexit, z) : nothing for z in exits] # add exit point
    end
    return Pexit
end


"""
    Some key functions must be adapted to the structure of the phase
    Others can be improved by exploiting the structure
"""

###
# Linear phase specific
###

# find_zeros_range(::LinearPhase) = nothing

exitpoints(::LinearPhase, ::Vector{NonOscillatoryBall}) = ComplexF64[]

###
# Polynomial specific functions
###

find_zeros_range(::PolynomialPhase,::Number) = (0.0, 10.0)

function findradius(G::PolynomialPhase, ξ, Cω, θ)  
    # specialised method for polynomials
    g(z) = evalphase(z,G)
    coef = coeffs(G.p)
    gpoly = Polynomial(0.0)
    for j = 0:degree(G) # compute coefficients of g(ξ + r e^(iθ))
        gpoly += Polynomial(coef[j+1] * [binomial(j,k) * ξ^(j-k) * cispi(k*θ) for k=0:j])
    end
    # construct G(r) = |g(ξ+re^{iθ})-g(ξ)|^2 - Cω^2
    Gξ = (gpoly-g(ξ)) * conj(gpoly-g(ξ)) - Cω^2
    allroots = roots(Gξ)
    posroots = allroots[real.(allroots) .> 0.0]
    # @show posroots
    rvals = real.(posroots[ abs.(imag.(posroots)) .< 0.01*(real.(posroots)) ])
    if isempty(rvals) 
        # Do bisection method to find roots.
        @warn "Doing bisection"
        Gb(r) = abs(g(ξ+r*cispi(θ)) - g(ξ))^2 - Cω^2
        return PlanBisection(posroots, Gb)
    else
        return minimum(rvals) # keep only positive roots
    end
end

function PlanBisection(posroots,G)
    # following same ideas as in https://github.com/AndrewGibbs/NumericalSteepestDescent/blob/master/src/root_finding/planBisection.m
    # If imaginary part using companion matrix is large, try something bisection
    # use roots from companion matrix as initial guesses for bisection.
    guesses = [0; real.(posroots); 2*maximum(real.(posroots))]
    for n in eachindex(guesses)[1:end-1]
        # we know G(0) < 0, so look for sign change in G(r)
        if G(guesses[n+1]) > 0
            # @show guesses[n], guesses[n+1]
            root = find_zero(G, (guesses[n], guesses[n+1]), Bisection(), xatol=1e-8)
            return root
        end
    end
    @warn "No positive roots found for radius after Bisection. Returning Inf."
    return Inf
end

function stop() end

# function exitpoints(G::PolynomialPhase, Ω :: Vector{NonOscillatoryBall})
#     # Look for local minima at the boundary of the non-oscillatory region
#     Pexit = ComplexF64[]
#     coef = coeffs(G.p)
#     J = length(coef)-1
#     for Ball in Ω # find exit points on the boundary of each ball
#         c, r = centre_and_radius(Ball)
#         # g(c + r e^(iθ)) = ∑ aj*(c + r exp(iθ))^j = ∑ bj * exp(ijθ) for j = 0:J
#         trig = Polynomial(0.0) # store g(c + r e^(iθ)) as Polynomial
#         for j = 0:J 
#             trig += Polynomial(coef[j+1] * [binomial(j,k) * c^(j-k) * r^k for k = 0:j])
#         end
#         tc = trig.coeffs 
#         # coefficients of first derivative of Im(g) = ∑ aj cos(jθ) + bj sin(jθ)
#         dtrig_cos = collect(0:J) .* real.(tc)
#         dtrig_sin = -collect(0:J) .* imag.(tc)
#         # coefficients of second derivative of Im(g) = ∑ aj' cos(jθ) + bj' sin(jθ)
#         ddtrig_cos = -collect(0:J).^2 .* imag.(tc)
#         ddtrig_sin = -collect(0:J).^2 .* real.(tc)
        
#         # find the roots of dtrig
#         tall = roots_trig_polynomial(dtrig_cos, dtrig_sin)
#         t = real.(tall[ abs.(imag.(tall)) .< 0.01 ])
#         # second-derivataive test to keep only maxima of Im g
#         dd(t) = sum( ddtrig_cos[k+1] * cos(k*t) + ddtrig_sin[k+1] * sin(k*t) for k = 0:J )
#         maxima = Float64[]
#         [dd(t) < 0.0 ? push!(maxima, t) : nothing for t in t]
        
#         # check if exit points are already in Ω (other non-oscillatory balls)
#         exits = [c+ r*cis(t) for t in maxima]
#         [!isinΩ(setdiff(Ω, [Ball]), z) ? push!(Pexit, z) : nothing for z in exits]
#         # [push!(Pexit, z) for z in exits]
#     end
#     return Pexit
# end

function exitpoints(G::PolynomialPhase, Ω :: Vector{NonOscillatoryBall})
    # Look for η s.t. Re g(η) - Re g(ξ) = 0
    Pexit = ComplexF64[]
    g(z) = evalphase(z,G)
    coef = coeffs(G.p)
    J = length(coef)-1
    for Ball in Ω # find exit points on the boundary of each ball
        c, r = centre_and_radius(Ball)
        trig = Polynomial(-real(g(c))) # store g(c + r e^(iθ)) as Polynomial
        for j = 0:J 
            trig += Polynomial(coef[j+1] * [binomial(j,k) * c^(j-k) * r^k for k = 0:j])
        end
        tc = trig.coeffs 
        trig_cos = real.(tc)
        trig_sin = -imag.(tc)
        tall = roots_trig_polynomial(trig_cos, trig_sin)
        t = real.(tall[ abs.(imag.(tall)) .< 0.01 ])
        # Keep only descent directions
        exits = Vector{ComplexF64}(undef,0)
        for θ in t
            z = c+r*cis(θ)
            if -imag(g(z)) < -imag(g(c)) push!(exits,z) end
        end
        # check if exit points are already in Ω (other non-oscillatory balls)
        [!isinΩ(setdiff(Ω, [Ball]), z) ? push!(Pexit, z) : nothing for z in exits]
    end
    return Pexit
end

###
# Square-root specific
###

function ballradius(G::SquareRootPhase, ξ, Cω; Nrays)
    # if stationary point is too large, drop contour tracing on a half-plane
    # we do it by giving it a small radius (which won't affect the algorithm)
    if abs(ξ) == Inf return NaN
    else return _compute_ballradius(G, ξ, Cω; Nrays)
    end
end

function find_zeros_range(G::SquareRootPhase, ::Number) 
    rmin = max(1.0 /(1-abs(G.b)) - 100., 0.0)
    rmax = 1.0 /(1-abs(G.b)) + 100 # min(1.0 /(1-abs(G.b)) + 100., 1e20)
    sort([0.0, rmax])
end

function exitpoints(G::SquareRootPhase, Ω::Vector{NonOscillatoryBall})

    """ DEV here """
    # special case, phase is constant on a half-plane
    # if abs(G.b) == 1.0 
    #     return [0.0im]
    # end

    c,r = centre_and_radius(Ω[1]) 
    g(z) = evalphase(z,G)
    trig   = θ -> imag(g(c + r*cis(θ)))
    dtrig  = θ -> ForwardDiff.derivative(trig,θ)  # first derivative of Im(g)
    ddtrig = θ -> ForwardDiff.derivative(dtrig,θ) # second derivative of Im(g)
  
    θ = find_zeros(dtrig, 0,  2π) # find the roots of dtrig
    if abs(im*G.a - c) > singular_tol(G)*r # check if branch point is away from the ball,
        maxima = Float64[] # second-derivative test
        [ddtrig(θ) < 0.0 ? push!(maxima, θ) : nothing for θ in θ] 
        return [c+ r*cis(θ) for θ in maxima]
    else # move along real axis
        # @info "moving along real axis"
        # xmin = max(0,real(c-r)) # assumes integration along [0,1]
        # xmax = min(1,real(c+r))
        return [c-r, c+r] 
    end
end

###
# Rational phase-specific
###

""" is radius ∞ as ω tends to zero? """

function find_zeros_range(::RationalPhase, ::Number) 
    return (0.0, 10.0) 
    # @show maxradius =  minimum(abs.( G.p .- ξ ))
    # return (0.0, maxradius)
end

function findradius(G::RationalPhase, ξ, Cω, θ)  
    αj  = G.coefs_analytic
    αpk = G.coefs_singular
    gpoly = Polynomial(0.0) # store coefficients of g(ξ + r e^(iθ)) as ∑ cm * r^m

    Q = Polynomial(1.0) # q(z) = ∏(z-zp)^Kp
    for (p,zp) in enumerate(poles(G)) 
        Kp = length(αpk[p])

        # singular part is regularised by q(z)
        regterm = Polynomial(0.0)
        for k = 1:Kp # regularised term at index p
            regterm += αpk[p][k] * Polynomial([binomial(Kp-k,m)*(ξ-zp)^(Kp-k-m)*cispi(m*θ) for m = 0:Kp-k])
        end
        otherterm = Polynomial(1.0) # other terms at indexes p' ≠ p
        for pp in setdiff(1:length(poles(G)), p)
            Kpp = length(αpk[pp])
            zpp = poles(G)[pp]
            otherterm *= Polynomial([binomial(Kpp,m)*(ξ-zpp)^(Kpp-m)*cispi(m*θ) for m = 0:Kpp])
        end
        gpoly += regterm * otherterm

        # store q(z)
        Q *= Polynomial([binomial(Kp,m)*(ξ-zp)^(Kp-m)*cispi(m*θ) for m = 0:Kp])
    end

    polyterm = Polynomial(0.0)
    for j = 0:degree(G)
        polyterm += αj[j+1] * Polynomial([binomial(j,m)*ξ^(j-m)*cispi(m*θ) for m=0:j])
    end
    gpoly += (polyterm - evalphase(ξ,G)) * Q

    # construct G(r) = |g(z)-g(ξ)|^2*|q(z)|^2 - Cω^2*|q(z)|^2, with z = ξ+re^{iθ}
    # gξ = evalphase(ξ,G)
    # G = (gpoly-gξ*Q) * conj(gpoly-gξ*Q) - Cω^2* Q*conj(Q)
    G = gpoly * conj(gpoly) - Cω^2* Q*conj(Q)

    rvals = roots(G)
    rvals = real.(rvals[ abs.(imag.(rvals)) .< 0.1])

    return minimum(rvals[rvals .> 0.0]) # keep only positive roots
end

function exitpoints(G::RationalPhase, Ω :: Vector{NonOscillatoryBall})
    # Look for η s.t. Re g(η) - Re g(ξ) = 0
    Pexit = ComplexF64[]
    g(z) = evalphase(z,G)
    αj  = G.coefs_analytic
    αpk = G.coefs_singular

    for Ball in Ω # find exit points on the boundary of each ball
        c, r = centre_and_radius(Ball)
        Q = Polynomial(1.0) # q(z) = ∏(z-zp)^Kp
        for (p,zp) in enumerate(poles(G))
            Kp = length(αpk[p])
            Q *= Polynomial([binomial(Kp,m)*(c-zp)^(Kp-m)*r^m for m = 0:Kp])
        end
        gpoly = 0.0
        for (p,zp) in enumerate(poles(G))
            Kp = length(αpk[p])
            regterm = Polynomial(0.0)
            for k = 1:Kp # regularised term at index p
                regterm += αpk[p][k] * Polynomial([binomial(Kp-k,m)*(c-zp)^(Kp-k-m)*r^m for m = 0:Kp-k])
            end
            for pp in setdiff(1:length(poles(G)), p)
                Kpp = length(αpk[pp])
                zpp = poles(G)[pp]
                regterm *= Polynomial([binomial(Kpp,m)*(c-zpp)^(Kpp-m)*r^m for m = 0:Kpp])
            end
            gpoly += regterm
        end
        polyterm = Polynomial(0.0)
        for j = 0:degree(G)
            polyterm += αj[j+1] * Polynomial([binomial(j,m)*c^(j-m)*r^m for m=0:j])
        end
        gpoly += polyterm * Q

        # Now, create the equation Re g(z) = Re g(ξ) for z(θ) = ξ + r*cis(θ) 
        gpoly -= g(c)*Q 
        trig_cos = real.(gpoly.coeffs)
        trig_sin = -imag.(gpoly.coeffs)
        # Solve the equation for θ (and remove imaginary solutions)
        tall = roots_trig_polynomial(trig_cos, trig_sin)
        t = real.(tall[ abs.(imag.(tall)) .< 0.01 ])

        # Keep only descent directions
        exits = Vector{ComplexF64}(undef,0)
        for θ in t
            z = c+r*cis(θ)
            if -imag(g(z)) < -imag(g(c)) push!(exits,z) end
        end
        # check if exit points are already in Ω (other non-oscillatory balls)
        [!isinΩ(setdiff(Ω, [Ball]), z) ? push!(Pexit, z) : nothing for z in exits]
        # [push!(Pexit, z) for z in exits]
    end
    return Pexit
end



