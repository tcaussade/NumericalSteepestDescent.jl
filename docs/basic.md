# Basic Usage

## Overview

NumericalSteepestDescent.jl is a Julia package for the numerical evaluation of highly oscillatory integrals of the form:

$$I = \int_{a}^{b} f(z)\exp(\mathrm{i}\omega g(z)) \, \mathrm{d}z$$

where:
- $g(z)$ is the **phase function**
- $f(z)$ is the **amplitude function**
- $\omega > 0$ is the **frequency parameter**
- $a, b$ are the integration endpoints (finite or infinite)

The package automatically implements a **regularised numerical steepest descent method**, which can efficiently evaluate these integrals even at high frequencies.

## Installation

Add the package to your Julia environment:

```julia
using Pkg
Pkg.add("NumericalSteepestDescent")
```

Or in the Julia REPL:

```julia
]add NumericalSteepestDescent
```

## Quick Start

### Example 1: Polynomial Phase

The simplest use case is integrating with a polynomial phase function. Here's a basic example:

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

### Example 2: Visualizing the Deformation

You can visualize how the integration contour deforms at different frequencies:

```julia
using NumericalSteepestDescent
using CairoMakie

Phase = PolynomialPhase([3, 5, 6, 2, 9, 5, 1, 4, 1, 3])

fig = Figure()
ax = Axis(fig[1, 1], title = "Steepest Descent Contour", 
          aspect = DataAspect(), xlabel = "Re(z)", ylabel = "Im(z)")

# Plot the contour deformation for ω = 5
quasiSDdeformation!(fig, ax, [-1, 1], Phase, 5.0)

fig
```

## Supported Phase Functions

### Polynomial Phase

For polynomial phases:

```julia
# Coefficients: [c₀, c₁, c₂, ...]
Phase = PolynomialPhase([1.0, 2.0, 3.0])
# This represents g(z) = 1 + 2z + 3z²
```

### Rational Phase

For rational functions (limited testing):

```julia
numerator = PolynomialPhase([1.0, 2.0])
denominator = PolynomialPhase([1.0, 1.0])
Phase = RationalPhase(numerator, denominator)
```

### Square Root Phase

For phase functions of the form $g(z) = \sqrt{z^2 + a^2} + bz$:

```julia
Phase = SquareRootPhase(a=1.0, b=0.5)
```

## Main Functions

### `nsd(f, Phase, endpoints, ω)`

Compute the highly oscillatory integral using the steepest descent method.

**Arguments:**
- `f`: Amplitude function (callable)
- `Phase`: Phase function object
- `endpoints`: Integration interval `[a, b]`
- `ω`: Frequency parameter

**Returns:** The value of the integral

### `quasiSDdeformation(endpoints, Phase, ω)`

Obtain the deformed contour for visualization or manual integration.

**Returns:** A `QuasiSDcontour` object containing the deformed path

### `plot_quasiSDdeformation(endpoints, Phase, ω)`

Plot the deformed contour (requires CairoMakie or similar plotting backend).

### `quasiSDdeformation!(fig, ax, endpoints, Phase, ω)`

Plot the deformed contour onto an existing figure and axis.

**Keyword Arguments:**
- `umax`: Maximum parameter value for contour representation
- `color_lim`: Limit for coloring by phase value
- `resolution`: Number of points for drawing the contour

## Key Concepts

### Steepest Descent Method

The method works by deforming the integration contour in the complex plane to follow the **steepest descent paths** of the phase function. This makes the oscillatory integral more tractable numerically.

### Quasi-Steepest Descent

A "regularised" or "quasi" steepest descent path is used to handle practical numerical issues while maintaining the theoretical benefits.

### Non-Oscillatory Regions

The method automatically identifies regions where the integrand is non-oscillatory and handles them appropriately.

## Performance Tips

1. **Choose appropriate frequency ranges**: The method is designed for moderate to high frequencies ($\omega \gtrsim 1$).

2. **Use analytical amplitude functions** when possible to avoid unnecessary function evaluations.

3. **For very high frequencies**, consider increasing numerical precision or using adaptive quadrature.

## Examples

More advanced examples are available in the `examples/` directory:

- `generic_example.jl`: Polynomial phase animation
- `cuspoid.jl`: Cuspoid catastrophe integrals
- `catastrophe.jl`: Catastrophe theory applications
- `airyintegral.jl`: Airy function evaluation
- `coalescence.jl`: Phase stationary points coalescence

## Related Work

This package is a Julia implementation of the regularised numerical steepest descent method. For the MATLAB version, see [PathFinder](https://github.com/AndrewGibbs/PathFinder).

## References

For theoretical background and detailed methodology, see the references in the [References](bib.md) page.
