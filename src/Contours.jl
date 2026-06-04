""" 
    ComplexContour

Struct for complex contour 

Type can be :finite, :infiniteSD or :finiteSD, and is needed to choose suitable quadature rules
Orientation tells in which way the contour is traversed.

If contour is :finite or :finiteSD, location should be two points 
If contour is :infiniteSD, location is point and valley angle.
"""

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
    Store traced contours in dictionary 
"""

function tracing_contours(G::AbstractPhase, points, Ω::Vector{NonOscillatoryBall};
                          δODE, δcoarse)
    # entrances     = Vector{ComplexF64}()
    # valley_points = Vector{ComplexF64}()
    # η_to_entrance = Dict{ComplexF64, ComplexF64}()
    # η_to_valley   = Dict{ComplexF64, ComplexF64}()
    ve = Vector{ComplexContour}()
    vv = Vector{ComplexContour}()
    vp = Vector{ComplexContour}()
    for η in points
        # println("Tracing from η=$η")
        h_end, pointtype = tracecontour_coarse(η, G, Ω; δODE, δcoarse) 
        if pointtype == :entrance
            γ = ComplexContour(:finiteSD, η, h_end)
            push!(ve,γ)
        elseif pointtype == :infvalley
            γ = ComplexContour(:infiniteSD, η, h_end)
            push!(vv,γ)
        elseif pointtype == :pole
            γ = ComplexContour(:infiniteSD, η, h_end)
            push!(vp,γ)
        end
    end
    return ve, vv, vp
end

"""
    Coarse tracing for SD contour at η
    Goal is determining if the contour goes to a valley or into the non-oscillatory region.
"""

function tracecontour_coarse(η, G::AbstractPhase, Ω; δODE, δcoarse)
    Pstat = get_Pstat(Ω)
    p1 = zero(ComplexF64)
    h1 = η # initial conditions 
    n = 0  # counter of iterations
    d = dist(h1, Pstat)
    g(z)   = evalphase(z,G)
    dg(z)  = evalphase_derivative(z,G)
    dg2(z) = evalphase_derivative2(z,G)
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
        d = dist(h1, Pstat) # update distance to stationary point
        
        if isinΩ(Ω, h1) # check if we entered non-oscillatory region - we have found an entrance
            # @info "Reached Ω from η1=$η to η2=$h1 in $n steps."
            # STORE AT THIS MOMENT WHATEVER WE NEED FOR LATER!!
            return (h1, :entrance)

        else # check if we entered a no-return region - we are at a valley
            bool, hval, typeval = isinValley(h1,G)
            if bool return (hval, typeval) end
        end
    end
end


""" 
    Check if z is inside the non-oscillatory region
"""

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
    Below are phase specific routines to identify valleys

tracing_contours() is adapted to special cases where we can exploit phase structure
"""

###
# Linear phase
###

valleyangle(θ, ::LinearPhaseFunction) = if θ≥0 return π/2 end # check z is in the upper half-plane

function tracecontour_coarse(::Any, ::LinearPhaseFunction, ::Any; δODE, δcoarse) 
    # we don't need to trace anything...
    return (cis(π/2), :infvalley)
end

###
# Polynomial phase
###

function isinValley(z, G::PolynomialPhase)
    if abs(z) > rstar_valley(G)
        v = valleyangle(angle(z), G)
        if ~(v isa Nothing) # is in no-return region at infinity
            hval = rstar_valley(G) * cis(v)
            return true, hval, :infvalley
            # should evaluate G(r,θ) here too!
        end
    end
    return false, nothing, nothing # break if not in valley at infinity
end

# This one also applies for RationalPhase
function valleyangle(θ, G::AbstractPhase) # find the angular valley of arg(z)
    J = degree(G)
    valleys = infvalleys(G)
    for v in valleys
        dist = minimum(abs.((θ-v) .- 2π*(-J:J)))
        if dist ≤ π/(2J) return v end
    end
end


###
# Sqrt phase
###

function tracecontour_coarse(η, G::SquareRootPhase, ::Any; δODE, δcoarse)
    # we don't need to trace anything...
    ξ = real(stationary_points(G)[1])
    v = π/2 * sign(real(η) - ξ) # we are assuming η is real
    return cis(v), :infvalley
end


###
# Rational phase
###

function isinValley(z, G::RationalPhase)

    # check if z is in valley at infinity
    if abs(z) > rstar_valley(G)
        v = valleyangle(angle(z), G)
        if ~(v isa Nothing) # is in no-return region at infinity
            hval = rstar_valley(G) * cis(v) 
            # return true, hval, :infvalley
            Gval = evaluate_noreturn_Ginf(abs(z), angle(v), G)
            if Gval > 0 return true, hval, :infvalley end   
        end

    # check if z is in valley at a pole
    else
        for (i,pole) in enumerate(poles(G))
            rp = rstar_pole(G)[i]
            
            if abs(z-pole) < rp # is near pole
                vp = polevalleys(G)[i]
                Kp = length(vp)
                for valley in vp # check if in angular valley at pole
                    θ = minimum(abs.( (angle(z-pole)-valley) .- 2π*(-2Kp:2Kp) ))
                    if θ < π/(2Kp)
                        hval = pole + cis(valley) * 1e-2
                        # add perturbation to distinguish between angular valleys
                        Gval = evaluate_noreturn_Gpole(abs(z-pole), θ, G; pole_idx = i)
                        if Gval > 0 return true, hval, :pole end    
                        # return true, hval, :pole           
                    end
                end
            end
        end
    end
    return false, nothing, nothing # break if not in valley at infinity
end
