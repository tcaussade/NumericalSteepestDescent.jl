# NumericalSteepestDescent.jl

[![Documentation Status](https://img.shields.io/badge/docs-in%20progress-orange)](https://tcaussade.github.io/NumericalSteepestDescent.jl/)

NumericalSteepestDescent.jl is a Julia package for the numerical evaluation of highly oscillatory integrals. 
The aim is to efficiently evaluate integrals of the form

$$
I = \int_{a}^b f(z)\exp(\mathrm{i}\omega g(z)) \mathrm{d}z
$$

where $g$ is phase, $f$ is amplitude, $\omega>0$ is a frequency parameter, and the endpoints $a$ and $b$ may be finite or infinite. 
Further, it is assumed that $I$ is a convergent integral and that $|f(z)|$ grows sub-exponentially as $|z|\to\infty$.

This package is an automatic implementation of a regularised numerical steepest descent method, 
which can be used without a deep understanding of the underlying mathematics
for particular choices of the phase:
- $g$ is a polynomial,
- $g$ is a rational function (not fully tested),
- $g(z) = \sqrt{z^2+a^2} + bz$, where $a>0$ and $b\in[-1,1]$. 

In the case of polynomial phase, the same method is also implemented in MATLAB, see https://github.com/AndrewGibbs/PathFinder/.
