
""" 
    Evaluate the integral a given contour type.

    We are using direct quadrature. Could be improved to adaptive quadrature
"""

function _does_contribute(G,ω,γ; δquad) 
    # determine is contribution from contour γ is significant
    abs(cis(ω*evalphase(G, at(γ)))) < δquad
end

function integrate(γ::ComplexContour, f, G, ω, x, w, quadtype :: Symbol; δfine, δquad, atol)
    
    if _does_contribute(G,ω,γ; δquad) return 0.0 + 0.0im end

    if quadtype == :gaussian
        # choose appropriate integration method based on contour type
        if contour_type(γ) == :finite   
            @show _is_singular(G,γ)
            if _is_singular(G,γ)
                # println("singular quad!")
                return integrate_finite_hp(γ, f, G, ω, x, w)
            else
                return integrate_finite(γ, f, G, ω, x, w)
            end
        elseif contour_type(γ) == :infiniteSD
            return integrate_infiniteSD(γ, f, G, ω, x, w; δfine)
        elseif contour_type(γ) == :finiteSD
            return integrate_finiteSD(γ, f, G, ω, x, w; δfine, δquad)
        end
    elseif quadtype == :adaptive
        if contour_type(γ) == :finite   
            return integrate_finite_gk(γ, f, G, ω; atol)
        elseif contour_type(γ) == :infiniteSD || contour_type(γ) == :finiteSD
            return integrate_SD_gk(γ, f, G, ω; δfine, δquad, atol)
        end
    else 
        @error "quadtype $quadtype is not valid.\nshould be :gaussian or :adaptive" 
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
    Singular quadrature routine (if needed). Criteria should be added to specific phase functions
"""

_is_singular(::AbstractPhaseFunction, ::ComplexContour) = false #

function _is_singular(G::SquareRootPhaseFunction, γ::ComplexContour)
    L = abs.(at(γ) - to(γ))
    tol = 0.01 # arbitrary choice
    if G.a / L < tol 
        return true
    else
        return false
    end
end

const σ = 0.17 # Grading parameter


function integrate_finite_hp(γ::ComplexContour, f::Function, G::AbstractPhaseFunction, ω, x, w)
    # evaluate integral along finite straight line from a to b
    g(z) = evalphase(G, z)
    L   = abs(to(γ)-at(γ))
    @assert abs(at(γ)) < abs(to(γ)) 

    S = zero(ComplexF64)
    n = ceil(log(2*G.a*σ/(1-σ))/log(σ)) # number of layers
    for i in range(1,n-1)
        el1,el2 = at(γ) .+ L.*([σ^i, σ^(i-1)]) # mesh sub-interval
        h = trace_finite(el1,el2)
        S += 0.5*(el2-el1) * dot(w, f.(h.(x)).*cis.(ω*g.(h.(x))))
    end
    # smallest interval
    el1, el2 = at(γ) .+ L.*[0.0, σ^(n-1)]
    h = trace_finite(el1,el2)
    S += 0.5*(el2-el1) * dot(w, f.(h.(x)).*cis.(ω*g.(h.(x))))
    return S
end


""" 
    We use Gauss-Laguerre for infinite SD contours
"""

function integrate_infiniteSD(γ::ComplexContour, f::Function, G::AbstractPhaseFunction, ω, x, w; δfine)
    # evaluate integral along infinite SD path at η
    g(z)  = evalphase(G, z)
    dg(z) = evalphase_derivative(G, z)
    @show η = at(γ)
    h = points_on_SDcontour(η, G, x/ω; δfine)
    dh = im ./ dg.(h)
    cis(ω*g(η))/ω * dot(w, f.(h).*dh)    
end

function points_on_SDcontour(η, G::AbstractPhaseFunction, xvec::Vector; δfine, η0 = η)
    # solve X in g(X) = g(η) + i x/ω
    g(z)  = evalphase(G, z)
    dg(z) = evalphase_derivative(G, z)
    h = zeros(ComplexF64, length(xvec))
    # f(u) = g(u)-g(η)-im*xvec[1]
    h[1] = Roots.newton(u -> g(u)-g(η)-im*xvec[1], dg, η0) # x0 = η
    # h[1] = Roots.newton(f,dg, η0)
    for j in 2:length(xvec)
        h[j] = Roots.newton(u -> g(u)-g(η)-im*xvec[j], dg, h[j-1]; rtol = δfine) # x0 = h[j-1] 
    end
    return h
end

function points_on_SDcontour(η, G::LinearPhaseFunction, xvec::Vector; δfine)
    g(z)  = evalphase(G, z)
    return η .+ im * xvec 
end


""" 
    We use (possibly truncated) Gauss-Legendre for finite SD contours
"""

function integrate_finiteSD(γ::ComplexContour, f::Function, G::AbstractPhaseFunction, ω, x, w; 
                            δfine, δquad)
    # evaluate integral along finite SD path at η going to Ω
    g(z)  = evalphase(G, z)
    dg(z) = evalphase_derivative(G, z)
    η = at(γ)
    umax = im * (g(η) - g(to(γ))) # pre-image of destination point
    # @assert abs(imag(umax)) < 1e-14 "umax = $umax should be real-valued \t η = $η"

    # possible truncation
    M  = 1.0 # pending: should be M = max(cis(ω * g(η))) for all η ∈ {Pstat, Pendp, Pexit}
    P = min(-log(δquad * M / abs(cis(ω*g(η)))), real(umax)*ω)

    p = trace_finite(0,P)
    h  = points_on_SDcontour(η, G, p.(x)/ω; δfine)
    dh = im ./ dg.(h)
    return cis(ω*g(η))/ω * dot(w, f.(h).*dh.*exp.(-p.(x))) * 0.5*P   
end


"""
    Adaptive quadrature routines based on Gauss-Kronrod quadrature rules
"""

# Default parameters for Gauss-Kronrod 
# Nodes and weights for G7/K15 over (-1,1) interval precomputed using QuadGK.jl
const x_gk  = [-0.9914553711208126, -0.9491079123427585, -0.8648644233597691, -0.7415311855993945, -0.5860872354676911, -0.4058451513773972, -0.2077849550078985, 0.0, 0.2077849550078985, 0.4058451513773971, 0.5860872354676911, 0.7415311855993945, 0.8648644233597691, 0.9491079123427584, 0.9914553711208125]
const w_gk  = [0.022935322010529256, 0.06309209262997842, 0.10479001032225017, 0.14065325971552592, 0.16900472663926788, 0.19035057806478559, 0.20443294007529877, 0.20948214108472793, 0.20443294007529877, 0.19035057806478559, 0.16900472663926788, 0.14065325971552592, 0.10479001032225017, 0.06309209262997842, 0.022935322010529256]
const wg_gk = [0.12948496616886981, 0.2797053914892767, 0.38183005050511887, 0.41795918367346907, 0.38183005050511887, 0.2797053914892767, 0.12948496616886981]

# Nodes and weights for G11/K23 over (-1,1) interval precomputed using QuadGK.jl
# const x_gk = [-0.9963696138895426, -0.978228658146057, -0.941677108578068, -0.8870625997680953, -0.816057456656221, -0.7301520055740494, -0.6305995201619651, -0.5190961292068118, -0.3979441409523776, -0.269543155952345, -0.13611300079936184, 0.0, 0.13611300079936184, 0.2695431559523449, 0.3979441409523776, 0.5190961292068117, 0.6305995201619652, 0.7301520055740494, 0.8160574566562211, 0.8870625997680954, 0.941677108578068, 0.9782286581460569, 0.9963696138895426]
# const w_gk = [0.009765441045960853, 0.027156554682104254, 0.04582937856442634, 0.06309742475037486, 0.07866457193222733, 0.09295309859690085, 0.1058720744813894, 0.11673950246104732, 0.12515879910031955, 0.13128068422980557, 0.1351935727998845, 0.13657779471111842, 0.1351935727998845, 0.13128068422980557, 0.12515879910031955, 0.11673950246104732, 0.1058720744813894, 0.09295309859690085, 0.07866457193222733, 0.06309742475037486, 0.04582937856442634, 0.027156554682104254, 0.009765441045960853]
# const wg_gk = [0.05566856711617362, 0.12558036946490478, 0.18629021092773423, 0.23319376459199045, 0.26280454451024665, 0.2729250867779004, 0.26280454451024665, 0.23319376459199045, 0.18629021092773423, 0.12558036946490478, 0.05566856711617362]

function integrate_finite_gk(γ::ComplexContour, f::Function, G::AbstractPhaseFunction, ω; atol)
    # evaluate integral along finite straight line from a to b
    g(z) = evalphase(G, z)
    a,b  = at(γ), to(γ)
    h    = trace_finite(a,b)
    ζ(u) = f(h(u)) * cis(ω*g(h(u)))
    # PROBLEM: specifying atol in quadgk() for some reason crashes the method...
    # @show int = quadgk_count(ζ,-1,1, atol = atol) 
    # using the QuadGK.jl here is much more efficient than handmade versions...
    return 0.5*(b-a) * quadgk(ζ,-1,1, atol = atol)[1] 
end

function integrate_SD_gk(γ::ComplexContour, f::Function, G::AbstractPhaseFunction, ω; δfine, δquad, atol)
    # evaluate integral along infinite or finite SD path at η using truncated GK
    g(z)  = evalphase(G, z)
    @show η = at(γ)

    M  = 1.0 # pending: should be M = max(cis(ω * g(η))) for all η ∈ {Pstat, Pendp, Pexit}
    umax = contour_type(γ) == :infiniteSD ? Inf : im * (g(η) - g(to(γ)))

    P = min(-log(δquad * M / abs(cis(ω*g(η)))), real(umax)*ω) # domain truncation
    # p = trace_finite(0,P)

    # Do Gauss-Kronrod routine by hand
    val = do_quadgk(0,P,f,G,ω,η; atol, δfine)
    cis(ω*g(η))/ω * val[1] # cis(ω*g(η))/ω * dot(w, f.(h).*dh.*exp.(-p.(x))) * 0.5*P  
end

# Core evaluation for adaptive quadrature
function eval_gk(a,b,f,G,ω,η; δfine)
    p = trace_finite(a,b)
    dg(z) = evalphase_derivative(G, z)
    # choose starting point, should be hη(p(-1))
    h  = points_on_SDcontour(η, G, p.(x_gk)/ω; δfine, η0 = η)
    dh = im ./ dg.(h)
    ζx =  f.(h) .* dh .* exp.(-p.(x_gk)) 
    Igk = sum(w_gk .* ζx)
    Ig  = sum(wg_gk .* ζx[2:2:end])
    abs_error = abs(Ig-Igk)
    # rel_error = abs_error / abs(Igk)
    return Igk * 0.5(b-a), abs_error #, rel_error
end

# Adapt if necessary
function do_quadgk(a,b,f,G,ω,η; atol, δfine)
    # @show (a,b)
    val, abserror = eval_gk(a,b,f,G,ω,η; δfine)
    if abserror > atol # && relerror > rtol
        mid = 0.5*(a+b)
        # should return the sum of the two intervals
        val1, e1 = do_quadgk(a,mid,f,G,ω,η; atol, δfine)
        val2, e2 = do_quadgk(mid,b,f,G,ω,η; atol, δfine)
        return val1+val2, e1+e2
    end
    return val, abserror
end


