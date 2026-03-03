using FastGaussQuadrature
using QuadGK
# using BenchmarkTools

""" 
    TASKS:
- use trapezoidal to evaluate residues
- try quadgk or gauleg to evaluate winding number
"""

# winding number
a = 0.0im
pts = [-1, -2-im, 2, 2+im, im, -1]


s = begin
    x,w = gausslegendre(5)
    s = 0.0
    for i in 1:length(pts)-1
        x0,x1 = pts[i], pts[i+1]
        h(t) = 0.5*(x1+x0 - (x1-x0)*t)
        s += 0.5 * (x1-x0) * sum(w.* 1.0./(h.(x) .- a)) / (2π*im)
    end
    s
end
abs(s)
sgk = quadgk(z -> 1/(z-a), pts, rtol = 5e-1)[1] / (2π*im)



# residue
ω = 50
εvals = [0.1,0.5,1]
g(z) = z + 1/z
f(z) = cis(ω*g(z))

# n = 10
# xt = 2π * (range(0,1,n)[1:end-1] .+ 0.5/(n-1)) # trapezoidal nodes
xgk = [1,im,-1,-im,1]

res_trapezoidal = Vector{ComplexF64}(undef,0)
res_quadgk = Vector{ComplexF64}(undef,0)
for ε in εvals
    @show dx = 1/exp(g(ε)) # adaptive trapezoidal
    @show n = ceil(1/dx + 1) |> Int
    xt = 2π * (range(0,1,n)[1:end-1] .+ 0.5/(n-1)) # trapezoidal nodes

    he(t) = ε*cis(t)
    r1 = sum(f.(he.(xt)) .* cis.(xt)) * im*ε / (2π*im) * 2π/(n-1)
    # @show ε, r1
    push!(res_trapezoidal, r1)
end
@show res_trapezoidal
for ε in εvals
    r2 =  quadgk(f,ε * xgk)[1] / (2π*im)
    push!(res_quadgk,r2)
end
@show res_quadgk

quadgk(z-> z+1/(z)+1/(z-0.1) + 1/z^2, xgk)[1] / (2*π*im)