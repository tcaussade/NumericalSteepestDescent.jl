# NumericalSteepestDescent.jl

## Abstract

NumericalSteepestDescent.jl is a Julia package for the numerical evaluation of highly oscillatory integrals. The aim is to efficiently evaluate integrals of the form

$$
I = \int_{a}^b f(z)\exp(\mathrm{i}\omega g(z)) \mathrm{d}z
$$

where $g$ is phase, $f$ is amplitude, $\omega>0$ is a frequency parameter, and the endpoints $a$ and $b$ may be finite or infinite. Further, it is assumed that $I$ is a convergent integral and that $|f(z)|$ grows sub-exponentially as $|z|\to\infty$.

This package is an automatic implementation of a regularised numerical steepest descent method, but it can be used without a deep understanding of the underlying mathematics. For a MATLAB implementation, see [[PathFinder, '25]](https://joss.theoj.org/papers/10.21105/joss.06902), and for a a full explanation of the underlying mathematics, the interested reader is referred to [[Gibbs, Hewett, Huybrechs, '24]](https://www.sciencedirect.com/science/article/pii/S0021999124000366?via%3Dihub).

## First example

The simplest use case is integrating with a polynomial phase function.

```julia
using NumericalSteepestDescent

# Define the phase: g(z) = 3 + 5z + 6z² + 2z³
Phase = PolynomialPhase([3, 5, 6, 2])

# Integration endpoints
a, b = -1.0, 1.0

# Amplitude function (can be any function)
f(z) = 1.0  # constant amplitude

# Frequency parameter
ω = 10.0

# Compute the integral
result = nsd(f, Phase, [a, b], ω)
```