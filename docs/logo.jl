using CairoMakie
logocolors = CairoMakie.Colors.JULIA_LOGO_COLORS

ω = 100
f(z) =  1/2 * cis(ω * z^2)
t = range(-1,1, length = 3000)
z1 = t * cis(π/64)
z2 = t * cis(π/32)
z3 = t * cis(π/4)

fig = Figure(backgroundcolor = :transparent, size = (500,500))
ax = Axis(fig[1, 1], title = "", backgroundcolor = :transparent)
lines!(t,real.(f.(z1)), color = logocolors.purple, linewidth = 8)
lines!(t,real.(f.(z2)), color = logocolors.red, linewidth = 8)
lines!(t,real.(f.(z3)), color = logocolors.green, linewidth = 8)
hidedecorations!(ax)

save("docs/src/assets/logo.svg", fig)