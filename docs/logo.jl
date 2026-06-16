using CairoMakie
logocolors = CairoMakie.Colors.JULIA_LOGO_COLORS

ω = 17.3
f(z) =  1 * cis(ω * z^2)
t = range(-1,1, length = 2000)
z1 = t * cis(π/128)
z2 = t * cis(π/32)
z3 = t * cis(π/4)

fig = Figure(backgroundcolor = :transparent)
ax = Axis(fig[1, 1], title = "", aspect = DataAspect())
lines!(t,real.(f.(z1)), color = logocolors.red, linewidth = 6)
lines!(t,real.(f.(z2)), color = logocolors.green, linewidth = 6)
lines!(t,real.(f.(z3)), color = logocolors.purple, linewidth = 8)
hidedecorations!(ax)
limits!(-1,1,-1,1.1)
fig

save("docs/src/assets/logo.svg", fig)