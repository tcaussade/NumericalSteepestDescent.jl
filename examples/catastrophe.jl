using PathFinder

f(z) = 1.0
ω    = 1.0

X = range(-2,2, length = 20)
Y = range(-2,2, length = 20)
Z = zeros(length(X), length(Y))

z = 0.0
amp(x,y) = 2*sqrt(π/3)*cis(4/27*z^3+1/3*x*z-π/4)

a,b = (-7π/12, π/12) .- π/24
for (i,x) in enumerate(X)
    for (j,y) in enumerate(Y)
        Phase = PathFinder.RationalPhaseFunction(im*[0,0,(z^2+x),0,2z,0,1],[0.],[[0.0, im*y^2/12]])
        try 
            Ψ, _ = PathFinder.integrate(a,b,f,Phase,ω; infcontour = [true,true])
            Z[i,j] = abs.(amp(x,y) * Ψ)
        catch
            println("(x,y,z) = ($x,$y,$z) did not throw a value")
        end
    end 
end

using WGLMakie
fig = Figure() #Figure(size = (500,400),)
ax = Axis(fig[1, 1], title = "Pearcey catastrophe integral", 
            xlabel = "Re", ylabel = "Im")
levelset = contourf!(ax,X,Y,Z; range(0, 5, 200), 
                        colormap = :jet, # colormap = :hot
                        extendlow = :auto, extendhigh = :auto)
Colorbar(fig[1,2], levelset)
# hidedecorations!(ax, grid = false)
fig