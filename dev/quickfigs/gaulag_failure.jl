""" This code has been AI generated """

using FastGaussQuadrature
using CairoMakie
using QuadGK

J = 2
f(z) = 1.0
g(z) = z^J
dg(z) = J*z^(J-1)
ω_vals = [5.0, 10.0, 20.0, 80.0]

x, w = gausslaguerre(20)
xref, wref = gausslaguerre(50)

function contour_min_radius(η)
    hη(u) = (g(η) + im * u)^(1 / J)
    u_end = max(10.0, 10 * abs(g(η)))
    u_vals = range(0.0, u_end, length = 501)
    return minimum(abs.(hη.(u_vals)))
end

function sd_integral(η, ω)
    hη(u) = (g(η) + im * u)^(1 / J)
    integrand(p) = im * f(hη(p / ω)) / dg(hη(p / ω)) * cis(ω * g(η)) / ω
    ival = sum(w .* integrand.(x))
    iref = sum(wref .* integrand.(xref))
    relerr = abs(ival - iref) / abs(iref)
    minrad = contour_min_radius(η)
    r_nonosc = (2π / ω)^(1 / J)
    enter_nonosc = minrad <= r_nonosc
    return ival, iref, relerr, minrad, enter_nonosc
end

δ_vals = range(0.0, 1.0, length = 41)
results_dict = Dict()

for ω in ω_vals
    results = [begin
            η = cis(-π/4) + δ
            ival, iref, relerr, minrad, enter_nonosc = sd_integral(η, ω)
            r_nonosc = (2π / ω)^(1 / J)
            (δ = δ, η = η, ival = ival, iref = iref, relerr = relerr,
                minrad = minrad, enter_nonosc = enter_nonosc,
                use_gausslaguerre = !enter_nonosc, r_nonosc = r_nonosc)
        end for δ in δ_vals]
    results_dict[ω] = results

    println("η sweep results for J=$J, ω=$ω")
    for r in results
        println("δ=$(round(r.δ, digits=3)), η=$(round(r.η, digits=3)), ival=$(round(r.ival, digits=5)), relerr=$(round(r.relerr, digits=4))")
    end
    println()
end

function sd_contour(η; n = 200)
    hη(u) = (g(η) + im * u)^(1 / J)
    u_end = max(10.0, 10 * abs(g(η)))
    u_vals = range(0.0, u_end, length = n)
    return hη.(u_vals)
end

fig = Figure(size = (600, 500), layout = (2, 2), colwidths = [0.4, 0.6], rowheights = [0.6, 0.4])
ax1 = Axis(fig[1, 1], xlabel = "Re(z)", ylabel = "Im(z)", title = "SD contour(s) in the complex plane", aspect = DataAspect())
ax2 = Axis(fig[2, 1], xlabel = "δ", ylabel = "Relative error", title = "Relative change vs δ", yscale = log10)

# Color palette for different ω values
colors = [:blue, :red, :green, :orange]
plot_handles = Vector{Any}()
plot_labels = String[]

# Plot contours for selected δ values for each ω
sample_idxs = [1, 11, 21, 31, 41]
for (ω_idx, ω) in enumerate(ω_vals)
    results = results_dict[ω]
    r_nonosc = (2π / ω)^(1 / J)
    
    # Plot selected contours for this ω
    for i in sample_idxs
        r = results[i]
        pts = sd_contour(r.η)
        color = colors[ω_idx]
        if i == sample_idxs[1]
            p = lines!(ax1, real.(pts), imag.(pts), color = color, label = "ω=$ω")
            push!(plot_handles, p)
            push!(plot_labels, "ω=$ω")
        else
            lines!(ax1, real.(pts), imag.(pts), color = color, alpha = 0.5)
        end
    end
    
    # Plot non-oscillation radius circle for this ω
    θ = range(0, 2π, length = 200)
    circle_pts = r_nonosc * cis.(θ)
    lines!(ax1, real.(circle_pts), imag.(circle_pts), color = colors[ω_idx], linestyle = :dash, alpha = 0.5)
end

Legend(fig[1, 2], plot_handles, plot_labels)

# Plot relative error curves for each ω
for (ω_idx, ω) in enumerate(ω_vals)
    results = results_dict[ω]
    lines!(ax2, δ_vals, getfield.(results, :relerr), color = colors[ω_idx], label = "ω=$ω")
end

axislegend(ax2; position = :rt)

fig

