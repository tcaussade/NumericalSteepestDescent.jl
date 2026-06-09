using CairoMakie
using FastGaussQuadrature

J = 3
ω = 10
f(z) = 1.0   
g(z) = z^J
dg(z) = J*z^(J-1)
hη(η) = u -> (g(η) + im*u)^(1/J)

# exit point
η = cis(π/(2J)) * (2π/ω)^(1/J)
# η = (2π/ω)^(1/J)

# gauleg parameters
M = 1.0
δquad = 1e-16


nmax = 50
gaulag_vals = zeros(ComplexF64, nmax)     
gauleg_vals = zeros(ComplexF64, nmax)
for n = 1:nmax
    # do gauss laguerre quadrature
    ilag = begin
        xlag,wlag = gausslaguerre(n)
        hη_vals = [hη(η)(p/ω) for p in xlag]
        dhη_vals = [im / dg(hp) for hp in hη_vals]
        sum( wlag .* f.(hη_vals).* dhη_vals ) * cis(ω*g(η))/ω
    end
    # println("n = $n, ilag = $ilag")
    gaulag_vals[n] = ilag

    # do gauss legendre quadrature
    ileg = begin
        
        xleg,wleg = gausslegendre(n)
        P = -log(δquad * M / abs(cis(ω*g(η))))
        xp = 0.5*((P+0) .+ (P-0)*xleg)
        hvals = [hη(η)(p/ω) for p in xp]
        dhvals = [im / dg(hp) for hp in hvals]
        cis(ω*g(η))/ω * sum(wleg .* f.(hvals).*dhvals.*exp.(-xp)) * 0.5*P   
    end
    # println("n = $n, ileg = $ileg")
    gauleg_vals[n] = ileg

end

refval = gaulag_vals[end]
gaulag_errors = abs.(gaulag_vals .- refval) #/abs(refval)
gauleg_errors = abs.(gauleg_vals .- refval) #/abs(refval)

fig = Figure()
ax = Axis(fig[1,1], xscale = sqrt, yscale = log10, 
    title = "Gauss-Laguerre vs Gauss-Legendre convergence", 
    xlabel = "Number of quadrature points", ylabel = "Error")
scatterlines!(ax, 1:nmax, gaulag_errors, label = "Gauss-Laguerre")
scatterlines!(ax, 1:nmax, gauleg_errors, label = "TruncatedGauss-Legendre")
axislegend(ax, position = :rt)
limits!(ax, 1, 35, 1e-18, 1)
fig

       