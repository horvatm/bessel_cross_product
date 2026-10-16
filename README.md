# Bessel cross products and radial modes

A small, self-contained C++17 extract from the full code associated with the
paper below, focused on calculating Bessel cross products and radial modes
in waveguide bends. Licensed under GPL-3.0-or-later.

Computes

\[
Z_\nu(k,r)=J_\nu(kr)Y_\nu(k)-Y_\nu(kr)J_\nu(k),
\]

including real orders and purely imaginary orders, and the Dirichlet radial
spectrum satisfying `Z_nu(k,q)=0`, `0<q<1`, at fixed `k>0`.
The outer radius is one. This is an **order spectrum**, not the conventional
annular eigenfrequency spectrum at fixed angular order.

## Associated paper

This repository contains the Bessel cross-product and radial-mode calculations
from the full code associated with the following paper:

Martin Horvat and Tomaž Prosen (2007). **The bends on a quantum waveguide
and cross-products of Bessel functions.** *Journal of Physics A: Mathematical
and Theoretical*, **40**(24), 6349–6379.
[doi:10.1088/1751-8113/40/24/006](https://doi.org/10.1088/1751-8113/40/24/006).

The paper develops the mode-number calculations and numerical schemes used here.
The [Figure 2 notebook](notebooks/paper_figure2.ipynb) reproduces its two
mode-number panels. The implementation is described in
[numerical methods](docs/NUMERICAL_METHODS.md).

Please cite the paper when using these methods or reproducing its results.
Machine-readable citation metadata is available in [CITATION.cff](CITATION.cff).

## Build and run

Requires CMake 3.16+ and a C++17 compiler with cylindrical Bessel special functions
(tested with GCC 15.2 on Linux). The C++ library has no external numerical-library
dependencies.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
ctest --test-dir build --output-on-failure
./build/spectrum 20 0.5 6 notebooks/spectrum
```

The example writes `notebooks/spectrum.csv` and `notebooks/spectrum_radial.csv`.
Arguments are `k q imaginary_count output_prefix`. With no arguments it uses
`20 0.5 6 spectrum`. Parent directories must already exist; output files are overwritten.

## Plot in Jupyter

```sh
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
pip install jupyterlab
jupyter lab notebooks/spectrum.ipynb
```

Run all cells. The notebook builds and runs the C++ example, plots the signed
order spectrum, and plots radial modes normalized with weight `1/r`.
It also checks the boundary residuals. Sample CSVs and executed plot outputs
are included so the demonstration can be inspected without compiling first.

## Reproduce the paper's Figure 2

Open [paper_figure2.ipynb](notebooks/paper_figure2.ipynb) for both mode-number
panels at `q=0.6`, including the dashed approximation from equation (26).
The executed notebook, 26-branch CSV data and PNG/SVG exports are included.
Set `recompute=True` in the notebook to regenerate the data, or run:

```sh
./build/paper_figure2 notebooks/paper_figure2.csv
```

Read [the manual](docs/MANUAL.md) for the API, formulas and numerical limitations,
[numerical methods](docs/NUMERICAL_METHODS.md) for the calculation schemes, and
[validation](docs/VALIDATION.md) for measured accuracy and checks.

## Minimal C++ use

```cpp
#include <bessel_cross/bessel_cross.hpp>
auto modes = bessel_cross::spectrum(20.0, 0.5, 6);
double z = bessel_cross::real_cross(2.3, 20.0, 0.7);
double u = bessel_cross::radial_value(modes.front(), 20.0, 0.7);
```

Link the `bessel_cross` CMake target, or compile your program together with
`src/cross_products.cpp` and `src/modes.cpp`, adding `-std=c++17 -Iinclude`.

## License

This project is licensed under **GPL version 3 or any later version**;
see [LICENSE](LICENSE). The same license covers the source code, documentation,
notebooks and generated example data and figures in this repository.
