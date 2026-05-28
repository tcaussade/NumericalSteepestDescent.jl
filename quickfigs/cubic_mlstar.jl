using CairoMakie
using ColorSchemes

# colorschemes[:viridis][1.0]
cmap = :viridis


a = Complex(-1)
G(z) = z^3 - 3*a*z
η = 1.0 # sqrt(abs(a))*im + 0.5
e1,e2 = [1+im*sqrt(3), 1-im*sqrt(3)]*0.5

function mlstar(η)
    @show v1,v2 = G(sqrt(a)), G(-sqrt(a))
    # θ1 = -angle(G(η)) + angle(v1)
    # θ2 = -angle(G(η)) + angle(v2)
    θ1 = angle(-G(η)+v1)
    θ2 = angle(-G(η)+v2)
    return θ1, θ2
end
mlstar(1)./π
mlstar(0)

mysqrt(z;θ) = cis(θ/2)*sqrt(Complex(z)*cis(-θ))
function ginv(s;θ1,θ2)

    # v(u) = real(a) > 0 ?
    #     0.5^(1/3) * (-u - im*mysqrt(u-2a^(3/2); θ=θ1) *  mysqrt(u+2a^(3/2); θ=θ2))^(1/3) :
    v(u) = ( 0.5 * (u + mysqrt(u-2a^(3/2); θ=θ1) *  mysqrt(u+2a^(3/2); θ=θ2)) )^(1/3)

    # g1(s) = (real(a)>0 ? -1 : 1) * ( +v(s) + a/v(s) )
    g1(s) = -v(s) - a/v(s)
    g2(s) = e1*v(s) + e2*a/v(s)
    g3(s) = e2*v(s) + e1*a/v(s)
    return g1.(s), g2.(s), g3.(s)
end
t = -80:0.5:80
Z = [x+im*y for x in t for y in t]

η = 0.673287 + 0.809952im
η = 1 + 1.6im
@show g(η)
# θ1,θ2 = mlstar(1)
θ1,θ2 = mlstar(η)
@show θ1/π, θ2/π
ginv1,ginv2,ginv3 = ginv(Z; θ1, θ2)
# ginv1,ginv2,ginv3 = ginv(Z; θ1 = π, θ2 = -π)


fig = Figure()
ax  = Axis(fig[1, 1], title = "", aspect = DataAspect(),
              xlabel = "Re", ylabel = "Im", xticks = -8:2:8, yticks = -8:2:8)

ms= 10
scatter!(reim.(ginv1), color = colorschemes[cmap][0.9], markersize = ms)
scatter!(reim.(ginv2), color = colorschemes[cmap][0.5], markersize = ms)
scatter!(reim.(ginv3), color = colorschemes[cmap][0.1], markersize = ms)
scatter!(reim.([-sqrt(Complex(a)),sqrt(Complex(a))]), color = "black", markersize = 15)
scatter!(reim.(η), marker = :star5, color = :red, markersize = 20)
limits!(-3,3,-3,3)
fig