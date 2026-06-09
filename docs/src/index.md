# NumericalSteepestDescent.jl

## Introduction

NumericalSteepestDescent.jl is a Julia package for the numerical evaluation of highly oscillatory integrals. The aim is to efficiently evaluate integrals of the form

```math
I = \int_{a}^b f(z)\exp(\mathrm{i}\omega g(z)) \mathrm{d}z
```

where ``g`` is phase, ``f`` is amplitude, ``\omega>0`` is a (potentially large) frequency parameter, and the endpoints ``a`` and ``b`` may be finite or infinite. Further, it is assumed that ``I`` is a convergent integral and that both ``f`` and ``g`` are slowly-varying functions.

This package is an automatic implementation of a *Regularised Numerical Steepest Descent method*, but it can be used without a deep understanding of the underlying mathematics. For a full explanation, the interested reader is referred to [1]. 

The simplest use case is integrating with a polynomial phase function. To evaluate the integral
```math
\int_{-1}^{1} e^{i\omega z^2} dz
```
try the following:
```@repl
using NumericalSteepestDescent
Phase = PolynomialPhase([0,0,1]) # Define the phase: g(z) = z²
a, b = -1.0, 1.0 # Integration endpoints
f(z) = 1.0  # # Amplitude function (can be any function)
ω = 1000.0 # Frequency parameter
nsd([a, b], f, Phase, ω)
```

For the case of arbitrary polynomial phase, there is also a MATLAB implementation available [2]. 

## References

[1] A. Gibbs, D.P. Hewett, D. Huybrechs, (2024) Numerical evaluation of oscillatory integrals via automated steepest descent contour deformation. Journal of Computational Physics, Volume 501, 112787, https://doi.org/10.1016/j.jcp.2024.112787.

[2] Gibbs, A., (2025). PathFinder: A Matlab/Octave package for oscillatory integration. Journal of Open Source Software, 10(114), 6902, https://doi.org/10.21105/joss.06902

