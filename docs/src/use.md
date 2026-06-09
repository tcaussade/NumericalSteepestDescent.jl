
# Basic use
The basic interface is the following
```@docs
nsd
```

**Example:** determine the quasi-SD contour deformation for an integral representation of the Airy function [DLMF, eq.(9.5.4)]
```@repl
Phase = PolynomialPhase(-im*[0,+2,0,1/3]) 
a, b = -π/3, π/3 # (∞exp(-im*π/3), ∞exp(im*π/3))
f(z) = 1.0 
ω = 1.0 # Frequency parameter
ai, fig = nsd([a, b], f, Phase, ω; infcontour = [true, true],
    plot_sd = true)
display(fig[1])
```

# Advanced use

There is also a number of adjustable parameters

|  Keyword |  Default | Meaning |
|---|---|---|
| `N` | ``25`` | Number of quadrature points used on each contour |
|  `Cball` | ``2\pi`` | Control maximum number of oscillations on non-oscillatory balls (and hence the ball radius)  |  
`Nrays`|  ``16`` | Number of rays used when determining the ball radius | 
|`δball`| ``10^{-3}`` | Governs when overlapping balls should be amalgamated | 
`δODE`| ``0.1`` | Governs the local step size in the ODE solver for SD path tracing | 
`δcoarse`| ``0.01`` |  Tolerance for the increment in corrector step of SD path tracing | 
`δfine`| ``10^{-13}`` | Tolerance of Newton solver for quadrature |
`δquad`|  ``10^{-16}`` | Tolerance for truncation and dropping contours if contribution is too small | 

