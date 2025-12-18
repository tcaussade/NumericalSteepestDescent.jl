
"""
    Struct representing a non-oscillatory ball in the complex plane
"""

struct NonOscillatoryBall 
    c :: ComplexF64 # centre, e.g. a stationary point
    r :: Float64    # radius
end
centre_and_radius(B::NonOscillatoryBall) = (B.c, B.r)

function ballradius(G::PolynomialPhaseFunction, ξ, Cω; Nrays)
    r = zeros(Nrays)
    for n = 1:Nrays
        ray(r) = ξ + r*cispi(2*n/Nrays)
        un(r) = abs(evalphase(G, ray(r)) - evalphase(G, ξ))^2 - Cω^2
        r[n] = minimum(find_zeros(un, 0,  10)) # shoudl we give some other value?
    end
    return minimum(r)
end

"""
    Construct the non-oscillatory region Ω

    G      : phase function (subtype of AbstractPhaseFunction)
    Cball  : constant defining the size of the non-oscillatory region
    ω      : frequency parameter
    δball  : overlap tolerance between non-oscillatory balls (default 1.0)

    Returns a Vector of NonOscillatoryBalls representing the non-oscillatory region.
"""

function NonOscillatoryRegion(G::AbstractPhaseFunction, ω; Cball, δball, Nrays)
    # δball = 1e-3 / 2 / max(degree(G)-2,1)
    Ω = Vector{NonOscillatoryBall}()
    # Assign a NonOscillatoryBall to each stationary point
    for ξ in G.ξ
        r = ballradius(G, ξ, Cball/ω; Nrays)
        push!(Ω, NonOscillatoryBall(ξ, r)) 
    end

    # Remove balls with significant overlap
    n = length(G.ξ)
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
    deleteat!(Ω, idx)
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

function exitpoints(G::AbstractPhaseFunction, Ω :: Vector{NonOscillatoryBall})
    Pexit = ComplexF64[]
    coef = coeffs(G.p)
    J = length(coef)-1
    for Ball in Ω # find exit points on the boundary of each ball
        c, r = centre_and_radius(Ball)
        trig = Polynomial(0.0) # trig = ∑ aj*(c + r exp(iθ))^j = ∑ bj * exp(ijθ) for j = 0:J
        for j = 0:J 
            trig += Polynomial(coef[j+1] * [binomial(j,k) * c^(j-k) * r^k for k = 0:j])
        end
        tc = trig.coeffs # convert to -imag(trig) = ∑ cj cos(jθ) + dj sin(jθ)
        dtrig_cos = -imag.(im * tc) .* collect(0:J)
        dtrig_sin = -imag.(tc) .* collect(0:J)
        
        # find the roots of the imaginary part of the derivative of trig
        tall = roots_trig_polynomial(dtrig_cos, dtrig_sin)
        t = Vector{Float64}()
        for ti in tall # keep only real roots 
            if abs(imag(ti)) < 1e-12 
                push!(t, real(ti))
            end
        end
        # second-derivataive test to keep only minima
        ddtrig_cos = -imag( -tc ) .* collect(0:J).^2
        ddtrig_sin = -imag( -im * tc ) .* collect(0:J).^2
        dd(t) = sum( ddtrig_cos[k+1] * cos(k*t) + ddtrig_sin[k+1] * sin(k*t) for k = 0:J )
           
        minima = Float64[]
        [dd(t) > 0.0 ? push!(minima, t) : nothing for t in t]
        
        # check if exit points are already in Ω (other non-oscillatory balls)
        exits = [c+ r*cis(t) for t in minima]
        [!isinΩ(setdiff(Ω, [Ball]), z) ? push!(Pexit, z) : nothing for z in exits]
    end
    return Pexit
end

