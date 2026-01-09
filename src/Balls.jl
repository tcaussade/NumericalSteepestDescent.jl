
"""
    Struct representing a non-oscillatory ball in the complex plane
"""

struct NonOscillatoryBall 
    c :: ComplexF64 # centre, e.g. a stationary point
    r :: Float64    # radius
end
centre_and_radius(B::NonOscillatoryBall) = (B.c, B.r)

function ballradius(G::AbstractPhaseFunction, ξ, Cω; Nrays)
    r = zeros(Nrays)
    g(z) = evalphase(G,z)
    for n = 1:Nrays
        ray(r) = ξ + r*cispi(2*n/Nrays)
        un(r) = abs(g(ray(r)) - g(ξ))^2 - Cω^2
        r_end = find_zeros_range(G) # heuristic choice - can we do better?
        r[n] = minimum(find_zeros(un, 0,  r_end)) 
    end
    return minimum(r)
end

find_zeros_range(G::PolynomialPhaseFunction) = 10
find_zeros_range(G::SquareRootPhaseFunction) = 1.0 /(1-abs(G.b)) + 10

"""
    Construct the non-oscillatory region Ω

    G      : phase function (subtype of AbstractPhaseFunction)
    Cball  : constant defining the size of the non-oscillatory region
    ω      : frequency parameter
    δball  : overlap tolerance between non-oscillatory balls (default 1.0)

    Returns a Vector of NonOscillatoryBalls representing the non-oscillatory region.
"""

function NonOscillatoryRegion(G::AbstractPhaseFunction, ω; 
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

function exitpoints(G::PolynomialPhaseFunction, Ω :: Vector{NonOscillatoryBall})
    Pexit = ComplexF64[]
    coef = coeffs(G.p)
    J = length(coef)-1
    for Ball in Ω # find exit points on the boundary of each ball
        c, r = centre_and_radius(Ball)
        # g(c + r e^(iθ)) = ∑ aj*(c + r exp(iθ))^j = ∑ bj * exp(ijθ) for j = 0:J
        trig = Polynomial(0.0) # store g(c + r e^(iθ)) as Polynomial
        for j = 0:J 
            trig += Polynomial(coef[j+1] * [binomial(j,k) * c^(j-k) * r^k for k = 0:j])
        end
        tc = trig.coeffs 
        # coefficients of first derivative of Im(g) = ∑ aj cos(jθ) + bj sin(jθ)
        dtrig_cos = collect(0:J) .* real.(tc)
        dtrig_sin = -collect(0:J) .* imag.(tc)
        # coefficients of second derivative of Im(g) = ∑ aj' cos(jθ) + bj' sin(jθ)
        ddtrig_cos = -collect(0:J).^2 .* imag.(tc)
        ddtrig_sin = -collect(0:J).^2 .* real.(tc)
        
        # find the roots of dtrig
        tall = roots_trig_polynomial(dtrig_cos, dtrig_sin)
        t = Vector{Float64}()
        for ti in tall # keep only real roots 
            if abs(imag(ti)) < 1e-1 # anything with larger imag part is dicarded
                push!(t, real(ti))
            end
        end
        # second-derivataive test to keep only maxima of Im g
        dd(t) = sum( ddtrig_cos[k+1] * cos(k*t) + ddtrig_sin[k+1] * sin(k*t) for k = 0:J )
        maxima = Float64[]
        [dd(t) < 0.0 ? push!(maxima, t) : nothing for t in t]
        
        # check if exit points are already in Ω (other non-oscillatory balls)
        exits = [c+ r*cis(t) for t in maxima]
        [!isinΩ(setdiff(Ω, [Ball]), z) ? push!(Pexit, z) : nothing for z in exits]
        # [push!(Pexit, z) for z in exits]
    end
    return Pexit
end

function exitpoints(G::SquareRootPhaseFunction, Ω::Vector{NonOscillatoryBall})
    c,r = centre_and_radius(Ω[1]) 
    g(z) = evalphase(G,z)
    trig   = θ -> imag(g(c + r*cis(θ)))
    dtrig  = θ -> ForwardDiff.derivative(trig,θ)  # first derivative of Im(g)
    ddtrig = θ -> ForwardDiff.derivative(dtrig,θ) # second derivative of Im(g)
  
    θ = find_zeros(dtrig, 0,  2π) # find the roots of dtrig
    if G.a > r # !isempty(θ) 
        maxima = Float64[] # second-derivative test
        [ddtrig(θ) < 0.0 ? push!(maxima, θ) : nothing for θ in θ] 
        return [c+ r*cis(θ) for θ in maxima]
    else # if branch point is inside ball, move along real axis
        return [c-r, c+r]
    end

end
