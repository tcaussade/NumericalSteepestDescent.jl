using PathFinder
# using Polynomials

f(z) = 1.0
ω    = 1.0

X = range(-10,10, length = 100)
Y = range(-10,2, length = 100)
Z = zeros(length(X), length(Y))

a,b = (9π/8, π/8)
a,b = π/1, 0.0
for (i,x) in enumerate(X)
    for (j,y) in enumerate(Y)
        Phase = PolynomialPhaseFunction([0, x, y, 0, 1])
        Ψ, _ = integrate(a,b,f,Phase,ω; infcontour = [true,true])
        Z[i,j] = abs.(Ψ)
    end 
end

fig = Figure(size = (500,400),)
ax = Axis(fig[1, 1], title = "Pearcey catastrophe integral", 
            xlabel = "Re", ylabel = "Im")
            #xticks = -xmin:1:xmax, yticks = -xmin:1:)
levelset = contourf!(ax,X,Y,Z; levels = 20, #levels = range(0.0, 2.5, 200), 
                        colormap = :hot, extendlow = :auto, extendhigh = :auto)
Colorbar(fig[1,2], levelset)
hidedecorations!(ax, grid = false)
fig
