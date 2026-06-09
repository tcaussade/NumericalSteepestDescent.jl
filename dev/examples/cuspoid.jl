using NumericalSteepestDescent
using CairoMakie

f(z) = 1.0
ω    = 1.0

X = range(-7,7, length = 120)
Y = range(-7,2, length = 120)
Z = zeros(length(X), length(Y))

a,b = π/1, 0.0
for (i,x) in enumerate(X)
    for (j,y) in enumerate(Y)
        Phase = PolynomialPhase([0, x, y, 0, 1])
        @time Ψ = nsd([a,b],f,Phase,ω; infcontour = [true,true])
        Z[i,j] = abs.(Ψ)
    end 
end


fig = Figure() #Figure(size = (500,400),)
ax = Axis(fig[1, 1], title = "Pearcey catastrophe integral", 
            xlabel = "x", ylabel = "y")
levelset = contourf!(ax,X,Y,Z; levels = 200, 
                        colormap = :jet, # colormap = :hot
                        extendlow = :auto, extendhigh = :auto)
Colorbar(fig[1,2], levelset)
# hidedecorations!(ax, grid = false)
fig