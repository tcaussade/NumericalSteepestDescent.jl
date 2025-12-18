""" 
    Struct for complex contour 

    Type can be :finite, :infiniteSD or :finiteSD, and is needed to choose suitable quadature rules
    Orientation tells in which way the contour is traversed.

    If contour is :finite or :finiteSD, location should be two points 
    If contour is :infiniteSD, location is point and valley angle.
"""

# struct ComplexContour
#     location :: Vector{ComplexF64}
#     type :: Symbol # Finite / Infinite SD / Finite SD
#     orientation :: Symbol # 
#     function ComplexContour(location, type::Symbol, orientation::Symbol)
#         @assert type in (:finite, :infiniteSD, :finiteSD)
#         @assert orientation in (:positive, :negative)
#         new(location, type, orientation)
#     end
# end

struct ComplexContour
    contourtype :: Symbol # Finite / Infinite SD / Finite SD
    source :: ComplexF64
    destination :: ComplexF64 
    function ComplexContour(type, src, dst)
        @assert type in (:finite, :infiniteSD, :finiteSD)
        new(type, src, dst)
    end
end

contour_type(γ::ComplexContour) = γ.contourtype
_check_contour_type(γ::ComplexContour, s::Symbol ) = @assert γ.contourtype == s "Contour type mismatch: expected $s, got $(γ.contourtype)"
# contour_orientation(γ::ComplexContour) = γ.orientation

at(γ::ComplexContour) = γ.source
to(γ::ComplexContour) = γ.destination


"""
    Coarse tracing for SD contour at η
    Goal is determining if the contour goes to a valley or into the non-oscillatory region.
"""

function tracecontour_coarse(G::AbstractPhaseFunction, η, Ω; δODE, δcoarse)
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

        # determine if we have found entrance point or a valley
        if isinΩ(Ω, h1)    
            @info "Reached Ω from η1=$η to η2=$h1 in $n steps."
            # STORE AT THIS MOMENT WHATEVER WE NEED FOR LATER!!
            return (h1, :entrance)

        elseif isinValley(G,h1) # define some threshold for valley
            # if doublecheck_valley(G,h1)
            v = goes_to_valley(G, angle(h1))
            if v isa Nothing continue end
            # @show angle(h1)/π, v/π
            hvalley = G.rStar * cis(v)
            # @show v, angle(h1)
                @info "Reached valley region at $(v/π)π from η=$η in $n steps."
            return (hvalley, :valley)
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
    if abs(z) > G.rStar

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

""" Store traced contours in dictionary """

function tracing_contours(G::AbstractPhaseFunction, points, Ω::Vector{NonOscillatoryBall};
                          δODE, δcoarse)
    # entrances     = Vector{ComplexF64}()
    # valley_points = Vector{ComplexF64}()
    # η_to_entrance = Dict{ComplexF64, ComplexF64}()
    # η_to_valley   = Dict{ComplexF64, ComplexF64}()

    ve = Vector{ComplexContour}()
    vv = Vector{ComplexContour}()
    for η in points
        h_end, status = tracecontour_coarse(G, η, Ω; δODE, δcoarse)
        # @show h_end, status, nmax
        if status == :entrance
            # push!(entrances, h_end) # (h_end, status))
            γ = ComplexContour(:finiteSD, η, h_end)
            # η_to_entrance[η] = h_end
            push!(ve,γ)
        elseif status == :valley
            # push!(valley_points, h_end)
            γ = ComplexContour(:infiniteSD, η, h_end)
            # η_to_valley[η] = h_end
            push!(vv,γ)
        end
    end
    # return entrances, η_to_entrance, valley_points, η_to_valley
    # return η_to_entrance, η_to_valley
    return ve, vv
end

