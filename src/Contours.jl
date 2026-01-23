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
    g(z)   = evalphase(G,z)
    dg(z)  = evalphase_derivative(G, z)
    dg2(z) = evalphase_derivative2(G, z)
    while d > 0 
        n+=1
        # predictor
        step1 = 2 * abs(dg(h1)^2/dg2(h1))
        step1 = isnan(step1) ? Inf : step1  # fix NaN instability
        p2 = p1 + δODE * min(step1, abs(dg(h1)) * d) # adaptative step
        h2 = h1 + (p2-p1) * im / dg(h1) # ode_iteration
        # corrector - ensure we are following the SD contour
        rtol = δcoarse * d
        # println("iterating, h1 = $(abs(h1)) * cis($(angle(h1)/π))π")
        h1 = find_zero((h->g(h)-g(η)-im*p2,dg),h2,Roots.Newton(); rtol = rtol)
        p1 = p2

        # determine if we have found entrance point or a valley
        if isinΩ(Ω, h1)    
            # @info "Reached Ω from η1=$η to η2=$h1 in $n steps."
            # STORE AT THIS MOMENT WHATEVER WE NEED FOR LATER!!
            return (h1, :entrance)

        else
            # check if in valley region
            bool_valley, hvalley = isinValley(G,h1)
            if bool_valley return (hvalley, :valley) end

            # check if near a pole
            bool_pole, hpole = isnearPole(G,h1)
            if bool_pole return (hpole, :pole) end
        end
    end
end


""" Check if z is inside the non-oscillatory region"""
function isinΩ(Ω::Vector{NonOscillatoryBall}, z)
    for Ball in Ω
        c, r = centre_and_radius(Ball)
        if abs(z - c) ≤ r
            return true
        end
    end
    return false
end

""" 
    Check if z is in a valley (PolynomialPhase) 
    For SqrtPhase and LinearPhase, we know a priori this is the case 
"""

# method for Polynomial Phase function
function isinValley(G::PolynomialPhaseFunction, z)
    # determine if z is in a valley region
    if abs(z) > G.rStar
        v = goes_to_valley(G, angle(z))
        if v isa Nothing # not in valley angular region - keep tracing
            return false, nothing
        end
        hvalley = G.rStar * cis(v) 
        return true, hvalley 
        
    end
    return false, nothing
end
function goes_to_valley(G::PolynomialPhaseFunction, θ) 
    # identifies the valley where θ is
    J = degree(G)
    valleys = G.v
    for v in valleys
        dist = minimum(abs.((θ-v) .- 2π*(-J:J)))
       # @show dist, θ/π, v/π
        if dist ≤ π/(2J)
            #@show v/π
            return v
        end
    end
end

# method for Rational Phase function
function isinValley(G::RationalPhaseFunction, z)
    # check if in valley at infinity
    if abs(z) > G.rstar_valley
        v = goes_to_valley(G, angle(z))
        if !(v isa Nothing) # not in valley angular region - keep tracing
            hvalley = G.rstar_valley * cis(v) 
            return true, hvalley 
        end    
    end
    return false, nothing
end
function goes_to_valley(G::RationalPhaseFunction, θ) 
    # identifies the valley where θ is
    J = length(G.analytic)-1 # degree(G.num)
    valleys = G.vinf
    for v in valleys
        dist = minimum(abs.((θ-v) .- 2π*(-J:J)))
       # @show dist, θ/π, v/π
        if dist ≤ π/(2J)
            #@show v/π
            return v
        end
    end
end
function isnearPole(G::RationalPhaseFunction, z)
    # check if in valley at a pole
    zp = near_pole(G,z)
    if !(zp isa Nothing) 
        return true, zp
    end
    # else, not in valley at a pole - keep tracing
    return false, nothing
end
function near_pole(G::RationalPhaseFunction, z)
    # identifies nearby poles
    for (i,zp) in enumerate(G.p)
        if abs(z-zp) < G.rstar_pole[i]
            Kp = length(G.vpole[i])
            θ  = angle(z - zp)
            for v in G.vpole[i] # check angle of valleys
                # @show v/π, abs.((θ-v) .- 2π*(-2Kp:2Kp))/π
                dist = minimum(abs.((θ-v) .- 2π*(-2Kp:2Kp)))
                if dist ≤ π/(2*Kp) return zp + cis(-v) * 10*eps() end
            end
        end
    end
end

function isnearPole(G::AbstractPhaseFunction, z)
    # if the function is not RationalPhase we don't expect any poles
    return false, nothing
end

# method for sqrt phase function
function isinValley(G::SquareRootPhaseFunction, z)
    if abs(z) == Inf
        return false, nothing
    end
    if abs(G.ξ[1]) == Inf 
        ξ = 1e12 * sign(real(G.ξ[1]))
        v = π/2 * sign(real(z) - ξ)
    else
        v = π/2 * sign(real(z) - G.ξ[1])
    end
    hvalley = cis(v) # this is patch fix!! re-do
    return true, hvalley
end

function isinValley(::LinearPhaseFunction, z)
    return true, cis(π/2)
end



""" Store traced contours in dictionary 
    Specialised method are provided for some phase functions.
"""

function tracing_contours(G::AbstractPhaseFunction, points, Ω::Vector{NonOscillatoryBall};
                          δODE, δcoarse)
    # entrances     = Vector{ComplexF64}()
    # valley_points = Vector{ComplexF64}()
    # η_to_entrance = Dict{ComplexF64, ComplexF64}()
    # η_to_valley   = Dict{ComplexF64, ComplexF64}()
    ve = Vector{ComplexContour}()
    vv = Vector{ComplexContour}()
    vp = Vector{ComplexContour}()
    for η in points
        h_end, status = tracecontour_coarse(G, η, Ω; δODE, δcoarse)
        if status == :entrance
            γ = ComplexContour(:finiteSD, η, h_end)
            push!(ve,γ)
        elseif status == :valley # at infinity
            γ = ComplexContour(:infiniteSD, η, h_end)
            push!(vv,γ)
        elseif status == :pole
            γ = ComplexContour(:infiniteSD, η, h_end)
            push!(vp,γ)
        end
    end
    # return entrances, η_to_entrance, valley_points, η_to_valley
    # return η_to_entrance, η_to_valley
    return ve, vv, vp
end

function tracing_contours(G::LinearPhaseFunction, points, ::Vector{NonOscillatoryBall};
        δODE, δcoarse)
    vv = Vector{ComplexContour}()
    ve = Vector{ComplexContour}()
    for η in points
        h_end = isinValley(G,η)[2] 
        γ = ComplexContour(:infiniteSD, η, h_end)
        push!(vv,γ)
    end
    return ve, vv
end

function tracing_contours(G::SquareRootPhaseFunction, points, ::Vector{NonOscillatoryBall};
        δODE, δcoarse)
    vv = Vector{ComplexContour}()
    ve = Vector{ComplexContour}()
    for η in points
        @show bool, h_end = isinValley(G, η) 
        # bool = true if we want to trace from there, false if ξ is too large
        if bool
            @show γ = ComplexContour(:infiniteSD, η, h_end)
            push!(vv,γ)
        end
    end
    return ve, vv
end

