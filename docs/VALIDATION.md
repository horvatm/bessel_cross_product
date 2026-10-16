# Validation report

Validated on Linux with GCC 15.2, CMake 4.2.3 and C++17, using a Release build
without fast-math. Python reference comparisons use mpmath 1.4.1.

## Reproduce

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
ctest --test-dir build --output-on-failure
python -m pip install -r requirements.txt
python tests/validate_mpmath.py
```

The optional Python script compares 19 cross products to independent Bessel
calculations with at least 90 decimal digits, increasing precision for the
exponential cancellation of imaginary-order functions. It samples every active
real and imaginary dispatcher branch. All samples pass the test's relative
error threshold of `1e-6`; the measured worst relative error is `4.84e-10`.
Near a zero, the script uses `max(abs(reference),1e-12)` as the denominator.
See [reference_results.txt](reference_results.txt) for the full sampled results.

It also recomputes all nine demonstration orders independently using mpmath
root finding and requires agreement within `1e-8` in order. These comparisons
validate this demonstration, not the completeness of the general root scanner.

CTest checks boundary identities, order symmetry, the imaginary zero limit,
the three-real-plus-six-imaginary mode count at `k=20, q=0.5`, descending signed
eigenvalues, inner boundary residuals, the radial ODE by finite differences,
a no-real-mode case and invalid geometry. It includes a regression for the
real-order asymptotic phase.

## Real-order asymptotic phase check

The real evanescent action requires the phase difference
`acosh(nu/(k*r))-acosh(nu/k)`. This differs from the `asinh` expression
used for imaginary orders.

At `nu=50, k=20, r=0.5`, the high-precision Bessel cross product is
`1.6206530037722779787e12`. The real-order asymptotic calculation agrees
with relative error `2.20e-12`. Substituting `asinh` would instead give
approximately `7.89901018357e10`, a relative error of about `0.951`.
The regression test guards against confusing the two phase expressions.

## Notebook verification

The notebook was executed from start to finish, including building and running
the C++ executable. At `k=20, q=0.5`, it produces three real orders and six
imaginary orders. The radial samples use 401 points per mode.

- Largest normalized boundary residual: `1.70e-10`.
- Largest weighted overlap error relative to the identity: `2.86e-9`.
- Spectrum and radial-mode figures were rendered and visually inspected.

The notebook includes its executed outputs, the source CSVs and PNG figures.
The plotted normalization uses trapezoidal integration with weight `1/r`.

The samples do not establish a uniform accuracy bound. Empirical
switches, asymptotic truncation and high-order overflow remain limitations.
See the manual before extending the computation to more extreme parameters.

## Paper Figure 2 reproduction

`notebooks/paper_figure2.ipynb` reproduces both panels at `q=0.6`, with the
published axis limits and the dashed equation (26) approximation. The C++ example
`paper_figure2` generates 26 branches at 601 wavenumbers from zero through
`sqrt(20000)`. The zero-wavenumber endpoints are analytic. The notebook inserts
18 independently calculated zero-order crossings to resolve the change between
imaginary and real orders on the signed square-root axis.

The notebook was executed completely and its PNG was visually compared with
Figure 2 on printed page 6354 of the associated paper. All sampled branches are
strictly increasing in signed eigenvalue. Six independent mpmath root checks,
covering real and imaginary orders near `k=10,50,120`, have a largest absolute
order difference of `3.49e-11`. The general C++ CTest suite also passes.
The reproduction uses a new sampling grid with the geometry, coordinates and
approximation specified in the paper.
