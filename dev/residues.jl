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

# exponential series truncation 
using SpecialFunctions
ω = 10
res = besseli(1,2*im*ω)
res = -im * besselj(1, -2ω)

si = 0.0
sj = 0.0
for m = 0:1
    # den = 1/(factorial(m)*factorial(m+1))
    # num = big((im*ω)^(2m+1))
    si += (im*ω)^(2m+1) / (factorial(m)*factorial(m+1))
    sj += (-1)^m * (-ω)^(2m+1) / (factorial(m)*factorial(m+1))
end
@show si, sj

"""
    try recurrence relation (Miller algorithm)
    Bessel function
"""

function eval_bessel(n,x)
    N = 100 # if n>N-5, perhaps one should use asymptotic formulas
    @assert n < N-5
    σ  = 0.0 # store the sum: σ = (2 ∑J_{2m}) / J0
    h = zeros(N) # stores h[k] = J_k/J_(k-1)
    for k = N-1:-1:1
        h[k] = x/(2k - x*h[k+1])
        σ = h[k] * (σ + 2*iseven(k))
    end 
    j0 = 1/(1+σ) # J0 = 1/(1+σ)
    if n == 0 return j0 end
    jk = zeros(n)
    jk[1] = h[1] * j0
    if n == 1 return jk[1] end

    for k = 2:n
        jk = h[k] * h[k+1]
    end
    return jk
end
x = 20
@time val = eval_bessel(1,x) # this is almost as good as SpecialFunctions package!!
@time exact = besselj(1,x)


asym(n,z) = sqrt(2/(π*z)) * cos(z-0.5*n*π - 0.25π)
@show ap = asym(1,x)
exact - ap

""" What if we didnt know that it is a Bessel function """
# Take g(z) = z+1/z and f(z) = 1.0
ω = π
g(z) = z + 1/z
F(z) = cis(ω*(g(z) - 1/z)) # in general, we can just substract the singularity

S = 0.0im
for m = 1:20
    Cm = (im*ω)^(m-1) * cis(0) / factorial(m-1)
    @show S += (im*ω)^m / factorial(m)
end
@show S
@show besseli(1,2*im*ω)