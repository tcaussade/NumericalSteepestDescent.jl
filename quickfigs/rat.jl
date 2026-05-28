using CairoMakie

g(z) = 1 / (z^2 + 1)

# grid in the complex plane
xmin, xmax, ymin, ymax = -2.0, 2.0, -2.0, 2.0
nx, ny = 300, 300
X = range(xmin, xmax, length = nx)
Y = range(ymin, ymax, length = ny)

# compute -Im(g(z)) on the grid (matrix with dims (length(X), length(Y)))
Z = [ -imag(g(x + im*y)) for x in X, y in Y ]

# --- Steepest-descent vector field and integrator ---
gprime(z) = -2z / (z^2 + 1)^2
gdoubleprime(z) = 2 * (3z^2 - 1) / (z^2 + 1)^3
fvec(z) = 1im * conj(gprime(z))

function rk4_step(z, dt)
	k1 = fvec(z)
	k2 = fvec(z + 0.5 * dt * k1)
	k3 = fvec(z + 0.5 * dt * k2)
	k4 = fvec(z + dt * k3)
	return z + dt * (k1 + 2k2 + 2k3 + k4) / 6
end

function integrate_path(z0; dt = 0.01, nsteps = 4000, direction = 1)
	z = z0
	path = ComplexF64[z]
	margin = 0.5
	for i in 1:nsteps
		step = direction * dt
		znew = rk4_step(z, step)
		if !(isfinite(real(znew)) && isfinite(imag(znew)))
			break
		end
		# stop if outside plotting window (+ margin) or near poles at ±im
		if real(znew) < xmin - margin || real(znew) > xmax + margin || imag(znew) < ymin - margin || imag(znew) > ymax + margin
			break
		end
		if abs(znew - im) < 1e-4 || abs(znew + im) < 1e-4
			break
		end
		push!(path, znew)
		z = znew
	end
	return path
end

# build figure
fig = Figure(size = (900, 600))
ax = Axis(fig[1, 1], xlabel = "Re(z)", ylabel = "Im(z)", title = "-Im(g(z))")
cf = contourf!(ax, X, Y, Z; levels = 20)
limits!(ax, xmin, xmax, ymin, ymax)
Colorbar(fig[1, 2], cf)

# find stationary point (minimize |g'| on the grid and refine with Newton)
function find_stationary_point(; tol = 1e-12)
	minval = Inf
	minidx = (1, 1)
	for (i, x) in enumerate(X)
		for (j, y) in enumerate(Y)
			z = x + im*y
			val = abs(gprime(z))
			if val < minval
				minval = val
				minidx = (i, j)
			end
		end
	end
	zg = X[minidx[1]] + im*Y[minidx[2]]
	# Newton refinement
	z = zg
	for k in 1:50
		gp = gprime(z)
		gpp = gdoubleprime(z)
		if abs(gpp) < 1e-14
			break
		end
		dz = gp / gpp
		znew = z - dz
		if abs(znew - z) < tol
			z = znew
			break
		end
		z = znew
	end
	return z
end

z0 = find_stationary_point()
println("stationary point z0 = ", z0)

# tangent vector for level sets Im(g)=const: dz/ds = conj(g')
flevel(z) = conj(gprime(z))

function rk4_level_step(z, ds)
	k1 = flevel(z)
	k2 = flevel(z + 0.5 * ds * k1)
	k3 = flevel(z + 0.5 * ds * k2)
	k4 = flevel(z + ds * k3)
	return z + ds * (k1 + 2k2 + 2k3 + k4) / 6
end

function integrate_level_path(z0; ds = 0.005, nsteps = 4000, direction = 1)
	z = z0
	path = ComplexF64[z]
	margin = 0.6
	for i in 1:nsteps
		step = direction * ds
		znew = rk4_level_step(z, step)
		if !(isfinite(real(znew)) && isfinite(imag(znew)))
			break
		end
		if real(znew) < xmin - margin || real(znew) > xmax + margin || imag(znew) < ymin - margin || imag(znew) > ymax + margin
			break
		end
		if abs(znew - im) < 1e-3 || abs(znew + im) < 1e-3
			break
		end
		push!(path, znew)
		z = znew
	end
	return path
end

# trace the level-set SD contours from the stationary point using g'' directions
g2 = gdoubleprime(z0)
phi = angle(g2)
thetas = [ -phi/2 + k * (π/2) for k in 0:3 ]
eps = 1e-4
sd_paths = Vector{Vector{ComplexF64}}()
for (k, θ) in enumerate(thetas)
	zinit = z0 + eps * cis(θ)
	pback = integrate_level_path(zinit; ds = 0.0025, nsteps = 8000, direction = -1)
	pforw = integrate_level_path(zinit; ds = 0.0025, nsteps = 8000, direction = 1)
	combined = vcat(reverse(pback), pforw[2:end])
	push!(sd_paths, combined)
	lines!(ax, real.(combined), imag.(combined); color = :blue, linewidth = 2)
end

# highlight endpoints and mark which SD path reaches them
endpoint_targets = [-1.0 + 0im, 1.0 + 0im]
for (tidx, zt) in enumerate(endpoint_targets)
	scatter!(ax, [real(zt)], [imag(zt)]; color = :black, markersize = 6)
	for p in sd_paths
		if any(abs.(p .- zt) .< 5e-3)
			lines!(ax, real.(p), imag.(p); color = :red, linewidth = 3)
		end
	end
end

# save and return the figure
save("quickfigs/rat_phase_sd_from_stationary.png", fig)
fig