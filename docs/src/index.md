# NumericalSteepestDescent.jl


NumericalSteepestDescent.jl is a Julia package for the numerical evaluation of highly oscillatory integrals. The aim is to efficiently evaluate integrals of the form

```math
I = \int_{a}^b f(z)\exp(\mathrm{i}\omega g(z)) \mathrm{d}z
```

where ``g`` is phase, ``f`` is amplitude, ``\omega>0`` is a frequency parameter, and the endpoints ``a`` and ``b`` may be finite or infinite. Further, it is assumed that ``I`` is a convergent integral and that ``|f(z)|`` grows sub-exponentially as ``|z|\to\infty``.

This package is an automatic implementation of a regularised numerical steepest descent method, but it can be used without a deep understanding of the underlying mathematics. For the case of arbitrary polynomial phase, there is also a MATLAB implementation available. See [[PathFinder, '25]](https://joss.theoj.org/papers/10.21105/joss.06902). 

The simplest use case is integrating with a polynomial phase function. To evaluate the integral
```math
\int_{-1}^{1} e^{i\omega(3+5z+6z^2+2z^3)} dz, 
```
try the following:
```@repl
using NumericalSteepestDescent
Phase = PolynomialPhase([3, 5, 6, 2]) # Define the phase: g(z) = 3 + 5z + 6z² + 2z³
a, b = -1.0, 1.0 # Integration endpoints
f(z) = 1.0  # # Amplitude function (can be any function)
ω = 1000.0 # Frequency parameter
nsd([a, b], f, Phase, ω) # Compute the integral
```

For a a full explanation of the underlying mathematics, the interested reader is referred to [[Gibbs, Hewett, Huybrechs, '24]](https://www.sciencedirect.com/science/article/pii/S0021999124000366?via%3Dihub).
