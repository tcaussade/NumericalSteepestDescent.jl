"""
    Evaluation of resiudes for phases with poles
"""

function winding_number(zp,pts)
    # zp is pole position
    # assume contour is union of straight lines with endpoints given in pts
    x,w = gausslegendre(5)
    s = zero(ComplexF64)
    for i in 1:length(pts)-1
        x0,x1 = pts[i], pts[i+1]
        h(t) = 0.5*(x1+x0 - (x1-x0)*t)
        s += 0.5 * (x1-x0) * sum(w ./ (h.(x) .- zp)) 
    end
    return s / (2π*im) # should be zero or one
end

function residue()
   """ idea: use trapezoidal """ 
   return 0.0
end

"""
    add_residues()
If the function has poles, should be included
Else, this step is not needed and we return 0.0
"""

add_residues(G::AbstractPhaseFunction, ::Any) = 0.0 

function add_residues(G::RationalPhaseFunction, γtot)

    # Substract original integration domain and deformed contour
    xs  = at(γtot[1])
    closedcurve = [xs] # store points
    for γ in γtot 
        @show to(γ)
        @show to(γ) ≈ 0.0im
        push!(closedcurve, to(γ)) 
    end
    push!(closedcurve, xs) # close loop
    # For each poles inside closedcurve, add residue
    res = 0.0
    for zp in poles(G)
        @show ind = winding_number(zp, closedcurve)
        if abs(ind) > 0.5  # zp is inside closedcontour
            @info "adding residue"
            res += residue()
        end
    end

    return 2π*im*res
end