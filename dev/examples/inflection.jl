using NumericalSteepestDescent
# using Polynomials

f(z) = 1.0
ω    = 1.0
k = 40

a,b = (9π/10, π/2)
G(x,y) = PolynomialPhase([0, 0, -y, 0, -x/2, 2/5])


# Compute amplitude
X = range(-10,10, length = 200)
Y = range(-10,10, length = 200)
Z = zeros(length(X), length(Y))
for (i,x) in enumerate(X)
    for (j,y) in enumerate(Y)
        Phase = G(x,y)
        Ψ, _ = nsd([a,b],f,Phase,ω; infcontour = [true,true])
        Z[i,j] = abs.(Ψ)
        Z0[i,j] = real.(Ψ .* cis(k*x)) # Compute approximate solution to Helmholtz equation
    end 
end

X0 = range(-5,5, length = 200)
Y0 = range(-1,1, length = 200)
Z0= zeros(length(X), length(Y))
for (i,x) in enumerate(X)
    for (j,y) in enumerate(Y)
        Phase = G(x,y)
        Ψ, _ = nsd([a,b],f,Phase,ω; infcontour = [true,true])
        Z0[i,j] = real.(Ψ .* cis(k*x)) # Compute approximate solution to Helmholtz equation
    end 
end

# Plots
using CairoMakie
fig = Figure() #  Figure(size = (800,1000)) # Figure(size = (500,600),)

gtop = fig[1,1] = GridLayout()
gbot = fig[2,1] = GridLayout()

ax = Axis(gtop[1, 1], title = "Inflection point problem", 
            xlabel = "Re", ylabel = "Im", aspect = DataAspect())
levelset = contourf!(ax,X,Y,abs.(Z); levels = range(0,2.5, 200), 
                        colormap = :viridis, 
                        extendlow = :auto, extendhigh = :auto)
Colorbar(gtop[1,2], levelset)
limits!(-10,10,-10,10)
colgap!(gtop, 10)
# Label(gtop[1,1], tellwidth = false)

ax0 = Axis(gbot[1, 1], title = "Helmholtz solution", 
            xlabel = "Re", ylabel = "Im", aspect = DataAspect())
helmholtz = contourf!(ax0,X0,Y0,Z0; levels = range(-2, 2, 100), 
                        colormap = :jet, # colormap = :hot
                        extendlow = :auto, extendhigh = :auto)
limits!(-5,5,-1,1)
Colorbar(gbot[1,2], helmholtz)
colgap!(gbot, 10)
# Label(gbot[1,1], tellwidth = false)

# resize_to_layout!(fig)
fig
