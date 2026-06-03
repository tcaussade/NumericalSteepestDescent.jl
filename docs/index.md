# NumericalSteepestDescent Documentation

NumericalSteepestDescent.jl is a Julia package for the numerical evaluation of highly oscillatory integrals. The aim is to efficiently evaluate integrals of the form

$$
I = \int_{a}^b f(z)\exp(\mathrm{i}\omega g(z)) \mathrm{d}z
$$

where $g$ is phase, $f$ is amplitude, $\omega>0$ is a frequency parameter, and the endpoints $a$ and $b$ may be finite or infinite. Further, it is assumed that $I$ is a convergent integral and that $|f(z)|$ grows sub-exponentially as $|z|\to\infty$.

This package is an automatic implementation of a regularised numerical steepest descent method, but it can be used without a deep understanding of the underlying mathematics. The same method is also implemented in MATLAB, see [[1]](https://joss.theoj.org/papers/10.21105/joss.06902), and for a a full explanation of the underlying mathematics, the interested reader is referred to [[2]](https://www.sciencedirect.com/science/article/pii/S0021999124000366?via%3Dihub).

