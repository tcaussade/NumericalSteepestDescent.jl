using NumericalSteepestDescent
using CairoMakie

framerate = 50

xvals = range(-2,2, step = 1/framerate)
G(x) = PolynomialPhase([0,-3x,0,1])
ω = 5
# γ0, infcontour = [π/2+π/4, π/4], [true,true]
γ0, infcontour = [5π/6, π/6], [true,true]
γ0, infcontour = [-1, 1], [false, false]

frame_iteration(x) = NumericalSteepestDescent.quasiSDdeformation!(fig, ax, γ0, G(x),ω; infcontour,
                                                    umax = 150,
                                                    color_lim = 100)

fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)
frame_iteration(first(xvals))
frame_iteration(0.05)
limits!(-4,4,-4,4)
fig


# CREATE ANIMATION 
record(fig, "quickfigs/cubic_coalescence.gif", xvals;
        framerate = 15) do x
    empty!(ax)
    ax.title = "a = $(round(x, digits = 1))"
    frame_iteration(x)
    limits!(-4,4,-4,4)
end

# Add parabola
a = -2
g = G(a).p
η = 1.0 # sqrt(abs(a))*im + 0.5
e1,e2 = [1+im*sqrt(3), 1-im*sqrt(3)]*0.5
θ = 0
mysqrt(z) = cis(θ/2)*sqrt(z*cis(-θ))
v(u) = begin
    ( 0.5 * (u + mysqrt(u^2-4a^3)) )^(1/3)
    # 0.5^(1/3) * (-u - im*sqrt(-u^2+4a^3))^(1/3)
end
ginv(s) = v(s) + a/v(s)
# ginv(s) = e1*v(s) + e2*a/v(s)
hη(u) = ginv(g(η) + im*u)
ρ = 1
Pρ(t) = ρ^2*(t^2-1) + im*2ρ^2*t
t = -12:0.001:12
parabola = Pρ.(t)    
mapped_parabola = hη.(Pρ.(t))
empty!(ax)
frame_iteration(a)
lines!(reim.(mapped_parabola), color = :grey, linewidth = 2)
limits!(-4,4,-4,4)
fig


## image and preimage for cubic
t = -20:1:20
Z = [x+im*y for x in t for y in t]
θ = 0

ginv(s) = e1*v(s) + e2*a/v(s)
fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)
lines!(reim.(ginv.(Z)))
limits!(-4,4,-4,4)
fig


####
# Gp = PolynomialPhase([3,5,6,2,9,5,1,4,1,3])
Gp = PolynomialPhase([0,0,1])
ω = 100
γ0 = [-1,1]

fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)
NumericalSteepestDescent.quasiSDdeformation!(fig, ax, γ0, Gp,ω; 
                                                    umax = 15,
                                                    color_lim = 10)

# frame_iteration(1)
δ = 2.5
limits!(-δ,δ,-δ,δ)
fig

#### Fractional monomial
J = 3/2
g(z) = z^J
hη