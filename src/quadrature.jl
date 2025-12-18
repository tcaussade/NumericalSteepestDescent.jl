
""" 
    Evaluate the integral a given contour type.

    We are using direct quadrature. Could be improved to adaptive quadrature
"""

function integrate(γ::ComplexContour, f, G, ω, x, w; δfine)
    if contour_type(γ) == :finite   
        return integrate_finite(γ, f, G, ω, x, w)
    elseif contour_type(γ) == :infiniteSD
        return integrate_infiniteSD(γ, f, G, ω, x, w; δfine)
    elseif contour_type(γ) == :finiteSD
        return integrate_finiteSD(γ, f, G, ω, x, w; δfine)
    end
end

""" 
    We use Gauss-Legendre for finite contours
"""

function trace_finite(a,b)
    # parametrisation of finite straight line from a to b
    u -> 0.5*((b+a) + (b-a)*u) # :: Function
end

function integrate_finite(γ::ComplexContour, f::Function, G::AbstractPhaseFunction, ω, x, w)
    # evaluate integral along finite straight line from a to b
    g(z) = evalphase(G, z)
    a,b  = at(γ), to(γ)
    h    = trace_finite(a,b)
    0.5*(b-a) * dot(w, f.(h.(x)).*cis.(ω*g.(h.(x))))
end

""" 
    We use Gauss-Laguerre for infinite SD contours
"""

function integrate_infiniteSD(γ::ComplexContour, f::Function, G::AbstractPhaseFunction, ω, x, w; δfine)
    # evaluate integral along infinite SD path at η
    g(z)  = evalphase(G, z)
    dg(z) = evalphase_derivative(G, z)
    η = at(γ)
    h = points_on_SDcontour(η, g, dg, x/ω; δfine)
    dh = im ./ dg.(h)
    cis(ω*g(η))/ω * dot(w, f.(h).*dh)    
end

function points_on_SDcontour(η, g, dg, xvec::Vector; δfine)
    # solve X in g(X) = g(η) + i x/ω
    h = zeros(ComplexF64, length(xvec))
    h[1] = Roots.newton(u -> g(u) - g(η) - im * xvec[1], dg, η) # x0 = η
    for j in 2:length(xvec)
        h[j] = Roots.newton(u -> g(u) - g(η) - im * xvec[j], dg, h[j-1]; rtol = δfine) # x0 = h[j-1] 
    end
    return h
end

# function trace_infiniteSDpath(G::AbstractPhaseFunction, η)
#     # parametrisation of infinite SD path at η
#     g(z) = evalphase(G, z)
#     u -> evalinverse(G, η, g(η) + im*u) # :: Function
# end


""" 
    We use (possibly truncated) Gauss-Legendre for finite SD contours
"""

function integrate_finiteSD(γ::ComplexContour, f::Function, G::AbstractPhaseFunction, ω, x, w)
    # evaluate integral along finite SD path at η going to Ω
    g(z)  = evalphase(G, z)
    dg(z) = evalphase_derivative(G, z)
    η = at(γ)
    h = complexified_nodes(η, g, dg, x/ω)
    dh = im ./ dg.(h)
    cis(ω*g(η))/ω * dot(w, f.(h).*dh)    
end

function eval_finiteSD(f::Function, G::AbstractPhaseFunction, γ::ComplexContour, ω, x, w)
    # evaluate integral along finite SD path from a to b
    _check_contour_type(γ, :finiteSD)
    g(z) = evalphase(G, z)
    h = γ.parametrisation
    η = h(0)
    @warn "Finite SD path evaluation not implemented yet"
    return 
end


