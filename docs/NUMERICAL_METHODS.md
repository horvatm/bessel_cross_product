# Numerical methods and scope

This project is a small, self-contained extract from the full code associated
with Martin Horvat and Tomaž Prosen's paper, *The bends on a quantum waveguide
and cross-products of Bessel functions*, J. Phys. A: Math. Theor. 40 (2007),
6349–6379. [DOI: 10.1088/1751-8113/40/24/006](https://doi.org/10.1088/1751-8113/40/24/006).
It focuses on Bessel cross products and radial modes in waveguide bends.

## Calculation schemes

The cross-product evaluator combines Taylor series, downward recurrence,
Hankel expansions, coupled real-order recurrence and large-order asymptotic
expansions. Empirical thresholds select a scheme for each combination of
order, wavenumber and radius. The manual lists these thresholds and their
limitations.

Real-order Bessel values use C++17 `std::cyl_bessel_j` and `std::cyl_neumann`.
The library requires no external numerical libraries. Numerical failures
are reported through exceptions; diagnostic dispatch variables are thread-local.
The real evanescent asymptotic phase uses `acosh`, while the imaginary-order
phase uses `asinh`. Independent reference comparisons and a regression test
check this distinction; see [validation](VALIDATION.md).

## Scope

The public API evaluates real- and imaginary-order cross products, finds the
Dirichlet radial order spectrum by bounded scans and bisection, and evaluates
radial mode functions. Examples demonstrate spectrum calculation, plotting
and reproduction of Figure 2 of the paper.

Transition matrices and full waveguide scattering calculations are outside
this package's scope. The implementation and examples are licensed under
GPL-3.0-or-later; see [LICENSE](../LICENSE).
