
"""
    Struct representing a non-oscillatory ball in the complex plane
"""

struct NonOscillatoryBall 
    c :: ComplexF64 # centre, e.g. a stationary point
    r :: Float64    # radius
end
centre_and_radius(B::NonOscillatoryBall) = (B.c, B.r)

function ballradius(G::PolynomialPhaseFunction, ξ, Cω; Nball = 16)
    r = zeros(Nball)
    for n = 1:Nball
        ray(r) = ξ + r*cispi(2*n/Nball)
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

function NonOscillatoryRegion(G::AbstractPhaseFunction, Cball, ω; δball = 1.0)
    δball = 1e-3 / 2 / max(degree(G)-2,1)
    Ω = Vector{NonOscillatoryBall}()
    # Assign a NonOscillatoryBall to each stationary point
    for ξ in G.ξ
        r = ballradius(G, ξ, Cball/ω)
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
            #@show trig
        end
        # convert to -imag(trig) = ∑ cj cos(jθ) + dj sin(jθ)
        tc = trig.coeffs
        # trig_cos = -imag.(tc)
        # trig_sin = -imag.(im * tc)
        # get coefficients of derivative of trig polynomial
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

function roots_trig_polynomial(a,b)

    """ 
    We are following "Computing the zeros, maxima and inflection points
        of Chebyshev, Legendre and Fourier series: solving
        transcendental equations by spectral interpolation
        and polynomial rootfinding" (Boyd, 2007)
    """

    # assume a vector of coefficients of cos(nθ) and b vector of coefficients of sin(nθ)
    # length(a) = 2N + 1
    @assert length(a) == length(b) "Coefficient vectors must be of the same length"
    # @assert abs(b[1]) ≈ 0.0 "The first coefficient of b_k (sin(kθ)) must be zero"
    J = length(a) - 1

    # compute Foorier companion matrix  
    h = zeros(ComplexF64, 2J+1)
    h[1:J]      = [a[J-j+1] + im*b[J-j+1] for j = 0:J-1]
    h[J+1]      = 2*a[1]
    h[J+2:2J+1] = [a[j-J+1] - im*b[j-J+1] for j = J+1:2J]

    B = zeros(ComplexF64, 2J, 2J)
    B[1:end-1, 2:end] = I(2J-1)
    B[end, :] = - h[1:end-1] / (a[end] - im*b[end]) 

    # convert eigenvalues to roots
    z = eigvals(B)
    θ = angle.(z) - im*log.(abs.(z)) 
    return θ
end


"""
    Determine the valleys at infinity
"""

function valleys(::AbstractPhaseFunction)
    # find valleys at infinity
    Pvalley = ComplexF64[]
    @warn "not implemented yet"
    return Pvalley
end

"""
    Tracing the SD contours
"""

dist(hη, Pstat :: Vector) = minimum(abs.(hη .- Pstat)) # dist(hn, P_statpoint)

function tracecontour(G::AbstractPhaseFunction, η, Ω; δODE = 1e-1, δcoarse = 1e-2)
    Pstat = get_Pstat(Ω)
    p1 = zero(ComplexF64)
    h1 = η # initial conditions 
    n = 0  # counter of iterations
    d = dist(h1, Pstat)
    while d > 0 
        n+=1

        # predictor
        step1 = 2 * abs(G.dp(h1)^2/G.dp2(h1))
        step1 = isnan(step1) ? Inf : step1  # fix NaN instability
        p2 = p1 + δODE * min(step1, abs(G.dp(h1)) * d) # adaptative step
        h2 = h1 + (p2-p1) * im / G.dp(h1) # ode_iteration

        # corrector - ensure we are following the SD contour
        rtol = δcoarse * d
        h1 = find_zero((h->G.p(h)-G.p(η)-im*p2,G.dp),h2,Roots.Newton(); rtol)
        p1 = p2
        #Roots.newton(x -> G.p(x) - G.p(η)-im*p2, G.dp, h2; rtol) # update value
        
        # determine if we have found entrance point or a valley
        if isinΩ(Ω, h1)
            
            @info "Reached Ω from η1=$η to η2=$h1 in $n steps."
            return (h1, :entrance)

        elseif isinValley(G,h1) # define some threshold for valley
            if doublecheck_valley(G,h1)
                v = goes_to_valley(G, angle(h1))
                if v isa Nothing continue end
                # @show angle(h1)/π, v/π
                hvalley = rstar * cis(v)
                # @show v, angle(h1)
                 @info "Reached valley region at $(v/π)π from η=$η in $n steps."
                return (hvalley, :valley)
            end
        end
    end
end

function isinΩ(Ω::Vector{NonOscillatoryBall}, z)
    for Ball in Ω
        c, r = centre_and_radius(Ball)
        if abs(z - c) < r
            return true
        end
    end
    return false
end

function isinValley(G::AbstractPhaseFunction, z)
    # determine if z is in a valley region
    if abs(z) > rstar

        return true
    end
    return false
end

function goes_to_valley(G::AbstractPhaseFunction, θ) 
    # identifies the valley where θ is
    J = degree(G)
    valleys = G.v
    for v in valleys
        dist = minimum(abs.((θ-v) .- 2π*(-J:J)))
       # @show dist, θ/π, v/π
        if dist < π/(2J)
            #@show v/π
            return v
        end
    end
end

function rvalley(G::AbstractPhaseFunction)
    # define threshold distance for valley region
    α = coeffs(G.p)
    J = length(α)-1
    β = [k*abs(α[k+1]) for k = 1:J-1]
    poly = Polynomial([β; -J*abs(α[J+1])/sqrt(2)]) 
    return maximum(real.(roots(poly))) # solution is the only positive root
end

function doublecheck_valley(G::AbstractPhaseFunction, h)
    r = abs(h)
    θ = angle(h)
    return true
    g(r,θ) = 0.0

    if g(r,θ) > 0 
        return true
    end
    return false
end


function tracing_contours(G::AbstractPhaseFunction, points, Ω::Vector{NonOscillatoryBall})
    entrances     = Vector{ComplexF64}()
    valley_points = Vector{ComplexF64}()
    η_to_entrance = Dict{ComplexF64, ComplexF64}()
    η_to_valley   = Dict{ComplexF64, ComplexF64}()
    for η in points
        h_end, status = tracecontour(G, η, Ω)
        # @show h_end, status, nmax
        if status == :entrance
            push!(entrances, h_end) # (h_end, status))
            η_to_entrance[η] = h_end
        elseif status == :valley
            push!(valley_points, h_end)
            η_to_valley[η] = h_end
        end
    end
    return entrances, η_to_entrance, valley_points, η_to_valley
end


