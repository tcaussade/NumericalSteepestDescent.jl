
""" 
    Evaluate the integral along any given complex contour
"""

function eval_finite(f::Function, G::AbstractPhaseFunction, γ::ComplexContour, ω, x, w)
    # evaluate integral along finite straight line from a to b
    # (x,w) are Gauss-Legendre points and weights
    _check_contour_type(γ, :finite)
    g(z) = evalphase(G, z)
    h = γ.parametrisation
    a,b = γ.parametrisation( -1. ), γ.parametrisation( 1. )
    0.5*(b-a) * dot(w, f.(h.(x)).*cis.(ω*g.(h.(x))))
end

function eval_infiniteSDpath(f::Function, G::AbstractPhaseFunction, γ::ComplexContour, ω, x, w)
    # evaluate integral along infinite SD path at η
    # (x,w) are Gauss-Laguerre points and weights
    _check_contour_type(γ, :infiniteSD)
    g(z) = evalphase(G, z)
    h = γ.parametrisation
    η = h(0)
    dh = u -> im / evalphase_derivative(G, h(u))
    cis(ω*g(η))/ω * dot(w, f.(h.(x/ω)).*dh.(x/ω))    
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


