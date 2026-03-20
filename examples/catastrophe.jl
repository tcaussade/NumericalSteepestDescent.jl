using NumericalSteepestDescent

"""
    Trying to evaluate 36.2.6 in DLMF to replicate Fig 36.3.6
"""

f(z) = 1.0
ω    = 1.0

X = range(-12,12, length = 50)
Y = range(-12,12, length = 50)
Z = zeros(length(X), length(Y))

z = 0.0
amp(x,y) = 2*sqrt(π/3)*cis(4/27*z^3+1/3*x*z-π/4)

a,b = (-7π/12, π/12) 
for (i,x) in enumerate(X)
    for (j,y) in enumerate(Y)
        println("Evaluating at ($x,$y,$z)")
        Phase = PathFinder.RationalPhaseFunction([0,0,(z^2+x),0,2z,0,1],[0.],[[0.0, y^2/12]])
        # try 
        @time Ψ = PathFinder.integrate([a,b],f,Phase,ω; infcontour = [true,true])
        Z[i,j] = abs.(amp(x,y) * Ψ)
        # catch
        #     println("(x,y,z) = ($x,$y,$z) did not throw a value")
        # end
    end 
end

using WGLMakie
fig = Figure() #Figure(size = (500,400),)
ax = Axis(fig[1, 1], title = "Elliptic umbilic catastrophe integral", 
            xlabel = "x", ylabel = "y")
levelset = contourf!(ax,X,Y,Z; levels = range(0, 4, 200), 
                        colormap = :jet, # colormap = :hot
                        extendlow = :auto, extendhigh = :auto)
Colorbar(fig[1,2], levelset)
# hidedecorations!(ax, grid = false)
fig