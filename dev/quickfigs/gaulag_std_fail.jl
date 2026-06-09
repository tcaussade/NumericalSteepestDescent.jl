using FastGaussQuadrature
using SpecialFunctions
using CairoMakie

function test_gaulag(δ,x,w)
    f(x) = 1/sqrt(Complex(x-δ))
    i1 = sum(w .* f.(x))
    i2 = erfc(sqrt(Complex(-δ)))*exp(-δ)*sqrt(π)
    return abs(i1-i2) / abs(i2)
end

δvals = range(0.0, 50.0, length = 10001)
x,w = gausslaguerre(20)

e1 = []
e2 = []
e3 = []
for δ in δvals
    # δ1 = δ
    # δ2 = δ + 0.1im
    # δ3 = δ *im
    δ1 = δ + 0.01im
    δ2 = δ + 0.1im
    δ3 = δ + im
    push!(e1, test_gaulag(δ1,x,w))
    push!(e2, test_gaulag(δ2,x,w))
    push!(e3, test_gaulag(δ3,x,w))  
end
fig = Figure()
ax = Axis(fig[1, 1],
          xlabel = "δ", ylabel = "Relative error", yscale = log10,
          title = "GauLag convergence for f(x) = 1/sqrt(x-δ)",
          xticks = 0:5:50)
# scatterlines!(δvals, e1, label = "real asingularity", color = :blue)
# scatterlines!(δvals, e2, label = "close to real singularity", color = :red)
# scatterlines!(δvals, e3, label = "pure imag singularity", color = :green)
scatterlines!(δvals, e1, label = "δ + 0.01i", color = :blue)
scatterlines!(δvals, e2, label = "δ + 0.1i", color = :red)
scatterlines!(δvals, e3, label = "δ + i", color = :green)
axislegend(ax)
fig

# Another experiment would be to have the singularity at \delta + i 
# for a few values of \delta (e.g. 1,5,10,...) and plot convergence of GauLag for increasing N.


Nvals = 1:100
imval = 0.
e1,δ1 = [], 10 + imval*im
e2,δ2 = [], 20 + imval*im
e3,δ3 = [], 30 + imval*im
for n in Nvals
    @assert typeof.([δ1, δ2, δ3]) == [ComplexF64, ComplexF64, ComplexF64]
    xn,wn = gausslaguerre(n)
    push!(e1, test_gaulag(δ1, xn, wn))
    push!(e2, test_gaulag(δ2, xn, wn))
    push!(e3, test_gaulag(δ3, xn, wn))
end
fig = Figure()
ax = Axis(fig[1, 1],
          xlabel = "√N", ylabel = "Relative error", yscale = log10, xscale = sqrt,
          title = "GauLag convergence for f(x) = 1/sqrt(x-δ)")
scatterlines!(Nvals, e1, label = "δ = $δ1", color = :blue)
scatterlines!(Nvals, e2, label = "δ = $δ2", color = :red)
scatterlines!(Nvals, e3, label = "δ = $δ3", color = :green)
axislegend(ax, position = :lb)
limits!(ax, 1,100, 1e-16, 10)
fig

