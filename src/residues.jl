"""
    Evaluation of resiudes for phases with poles
"""

using QuadGK

function winding_number(zp,pts)
    # zp is pole position
    # assume contour is union of straight lines with endpoints given in pts
    # x,w = gausslegendre(5)
    s = zero(ComplexF64)
    for i in 1:length(pts)-1
        x0,x1 = pts[i], pts[i+1]
        # h(t) = 0.5*(x1+x0 - (x1-x0)*t)
        # s += 0.5 * (x1-x0) * sum(w ./ (h.(x) .- zp)) 
        s += quadgk(t->1.0/(t-zp), x0, x1)[1]
    end
    # s / (2π*im)
    return s / (2π*im) # should be zero or one
end

"""
    winding_contour(γ,γ0)

Let γ be a quasi-SD deformation, and γ0 the original integration contour.
"""

function winding_contour(γsd, γ0)
    # Substract original integration domain and deformed contour
    # @show γsd[1], γ0
    @assert at(γsd[1]) ≈ γ0[1]  "$(at(γsd[1])) should be equal to $(γ0[1]))"

    γclose = Vector{ComplexF64}()
    [push!(γclose,z) for z in γ0] # original integration domain
    for γ in reverse(γsd)
        # we are approximating SD contour by straight lines.
        push!(γclose, at(γ))
    end
    # To-Do: use more points from coarse contour tracing
    return γclose
end