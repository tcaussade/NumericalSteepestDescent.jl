using CairoMakie
using ColorSchemes

# colorschemes[:viridis][1.0]
cmap = :viridis
a = Complex(-1)
G(z) = z^3 - 3*a*z
e1,e2 = [1+im*sqrt(3), 1-im*sqrt(3)]*0.5
ξ1,ξ2 = [-sqrt(Complex(a)),sqrt(Complex(a))]

function mlstar(η)
    @show a
    @show vm,vp = G(ξ1), G(ξ2)
    θm = angle(vm - G(η)) 
    θp = angle(vp - G(η))
    return θm, θp
end
mlstar(2im)./π
mlstar(0)

mysqrt(z;θ) = cis(θ/2)*sqrt(Complex(z)*cis(-θ))
function ginv(s;θm,θp)
    v(u) = ( 0.5 * (u + mysqrt(u-2a^(3/2); θ=-θm) *  mysqrt(u+2a^(3/2); θ=-θp)) )^(1/3)
    g1(s) = v(s) + a/v(s)
    g2(s) = -e1*v(s) - e2*a/v(s)
    g3(s) = -e2*v(s) - e1*a/v(s)
    return g1.(s), g2.(s), g3.(s)
end
t = -60:0.35:60
Z = [x+im*y for x in t for y in t]


### ADD MAPPED PARABOLA
function choose_branch(η)
    for i = 1:3
        h(u) = ginv(G(η) + im*u; θm, θp)[i]
        if abs(h(0) - η) < 1e-14
            @show h(0), η, h(0)-η
            return i
        end
    end
end


Cubic = PolynomialPhase([0,-3a,0,1])
Ω = NonOscillatoryRegion(Cubic, ω; Cball=2π, δball=1e-3,  Nrays=16)
exits = NumericalSteepestDescent.exitpoints(Cubic,Ω)

η = exits[1]
η = exits[2]

θm,θp = mlstar(η)
ginv1,ginv2,ginv3 = ginv(Z; θm, θp) # ginv1,ginv2,ginv3 = ginv(Z; θm = π, θp = -π)

idx = choose_branch(η)
hη(u) = ginv(G(η) + im*u; θm, θp)[idx]
@assert abs(hη(0) - η) < 1e-14

ρ = 0.8
fig = Figure()
ax  = Axis(fig[1, 1], title = "a=$a, η=$η", aspect = DataAspect(),
            xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)
create_fig!(ρ)
fig

function create_fig!(ρ)   
    tp = range(-15,15, length = 150)
    Pρ = hη.( ρ^2*(tp.^2.0.-1).+im*2*ρ^2*tp )

    ### ADD APOLLONIUS CIRCLE
    θAC = range(0, 2π, length = 150)
    AC = cis.(θAC) * 2*abs(ξ1) .+ 3*ξ1
        
    ms= 8
    scatter!(reim.(ginv1), color = colorschemes[cmap][0.9], markersize = ms)
    scatter!(reim.(ginv2), color = colorschemes[cmap][0.5], markersize = ms)
    scatter!(reim.(ginv3), color = colorschemes[cmap][0.1], markersize = ms)
    scatter!(reim.([ξ1,ξ2]), color = "black", markersize = 10)
    scatter!(reim.(η), color = :red, markersize = 10)
    limits!(-3,3,-3,3)
    lines!(reim.(Pρ), color = :red, linewidth = 3)
    lines!(reim.(AC), color = :blue, linewidth = 3)
end

### Overlap two parabolas in same plot, note that ML star "division"
# when \eta is exit point is the same.

a=+1
Cubic = PolynomialPhase([0,-3a,0,1])
Ω = NonOscillatoryRegion(Cubic, ω; Cball=2π, δball=1e-3,  Nrays=16)
exits = NumericalSteepestDescent.exitpoints(Cubic,Ω)

η1 = exits[3]
η2 = exits[4]

# η1 = exits[1]
# η2 = exits[2]

create_fig(ρ,η1,η2)

function create_fig(ρ,η1,η2)
    θm,θp = mlstar(η1) # or η2
    @show mlstar(η1).-mlstar(η2)
    ginv1,ginv2,ginv3 = ginv(Z; θm, θp) # ginv1,ginv2,ginv3 = ginv(Z; θm = π, θp = -π)

    idx1 = choose_branch(η1)
    hη1(u) = ginv(G(η1) + im*u; θm, θp)[idx1]
    @show abs(hη1(0) - η1)
    @assert abs(hη1(0) - η1) < 1e-14
    idx2 = choose_branch(η2)
    hη2(u) = ginv(G(η2) + im*u; θm, θp)[idx2]
    @assert abs(hη2(0) - η2) < 1e-14

    tp = range(-15,15, length = 150)
    Pρ1 = hη1.( ρ^2*(tp.^2.0.-1).+im*2*ρ^2*tp )
    Pρ2 = hη2.( ρ^2*(tp.^2.0.-1).+im*2*ρ^2*tp )
    ### ADD APOLLONIUS CIRCLE
    θAC = range(0, 2π, length = 150)
    AC = cis.(θAC) * 2*abs(ξ1) .+ 3*ξ1
        
    fig = Figure()
    ax  = Axis(fig[1, 1], title = "a=$a, η=$η", aspect = DataAspect(),
                xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)
    create_fig!(ρ)
    ms= 8
    scatter!(reim.(ginv1), color = colorschemes[cmap][0.9], markersize = ms)
    scatter!(reim.(ginv2), color = colorschemes[cmap][0.5], markersize = ms)
    scatter!(reim.(ginv3), color = colorschemes[cmap][0.1], markersize = ms)
    scatter!(reim.([ξ1,ξ2]), color = "black", markersize = 10)
    scatter!(reim.([η1,η2]), color = :red, markersize = 10)
    limits!(-3,3,-3,3)
    lines!(reim.(Pρ1), color = :red, linewidth = 3)
    lines!(reim.(Pρ2), color = :red, linewidth = 3)
    lines!(reim.(AC), color = :blue, linewidth = 3)

    return fig
end
|