# User and numerical manual

## Geometry and the meaning of a mode

Lengths are scaled by the outer radius, so the annulus is `q <= r <= 1`.
For physical outer radius `b` and inner radius `a`, use `q=a/b`, dimensionless
`k=b*k_physical`, and `r=r_physical/b`.

Separation of the Helmholtz equation gives

$$
r^2 R''+rR'+(k^2r^2-\lambda)R=0,\quad R(q)=R(1)=0.
$$

Positive `lambda=nu^2` gives a real Bessel order. Negative `lambda=-eta^2`
gives imaginary order `nu=i eta`. There are finitely many positive eigenvalues
and infinitely many negative ones. The angular equation is `Theta''+lambda*Theta=0`:
positive eigenvalues are oscillatory angular modes; negative ones are evanescent.
Additional angular boundary conditions are not imposed by this library.

In logarithmic radius `t=log(r)`, the problem becomes
`R_tt + k^2 exp(2t) R = lambda R` on `[log(q),0]`.
This explains the positive-order bound `nu<k`, the asymptotic imaginary-order
spacing `pi/abs(log(q))`, and radial orthogonality with weight `dr/r`.

## Public API

Include `bessel_cross/bessel_cross.hpp`; all functions are in `bessel_cross`.

| Function | Meaning |
|---|---|
| `real_cross(nu,k,r)` | `J_nu(kr)Y_nu(k)-Y_nu(kr)J_nu(k)` |
| `imaginary_cross(eta,k,r)` | Same cross product at order `i*eta`, returned as a real number |
| `spectrum(k,q,imaginary_count)` | All detected real orders, followed by the requested imaginary orders |
| `radial_value(mode,k,r)` | Evaluate the corresponding cross product |

Orders are accepted with either sign and evaluated by symmetry. Inputs must be
finite, `k>0`, and `0<r<=1`. The imaginary-order zero limit delegates to the
real-order calculation. At `r=1`, both functions return exactly zero.
There is no imaginary-order amplitude scaling: the combination itself is real,
even though the individual Bessel functions are complex. The Taylor form is

$$
Z_{i\eta}(k,r)=-\frac{2}{\pi\eta}\operatorname{Im}
\left[r^{i\eta}V(\eta,kr/2)\overline{V(\eta,k/2)}\right],
\quad
V(\eta,x)=\sum_{m=0}^\infty\frac{(-x^2)^m}{m!(1+i\eta)_m}.
$$

`Mode` holds `order` (nonnegative magnitude), `eigenvalue` (signed square),
and `imaginary` (whether the order is imaginary). The vector is sorted by
descending eigenvalue. Values own their storage and require no manual deletion.
Functions throw `std::invalid_argument` on invalid inputs and
`std::runtime_error` on detected numerical failures. Standard-library Bessel
functions can also throw. Allocation failure can propagate.

Raw cross products are not normalized. For a mode `R`, use
`N=sqrt(integral_q^1 R(r)^2 dr/r)` and divide by `N`.
The notebook demonstrates sampled trapezoidal normalization. Values at the
inner boundary retain the computed residual; they are not clamped to zero.

## Preserved numerical schemes

The real-order dispatcher uses:

1. Standard-library Bessel products for `k<1000` and `abs(nu)<k`.
2. Large-order asymptotics for `abs(nu)>2*k`.
3. Hankel expansions for `abs(nu)<=0.1*k*r` outside the direct branch.
4. Coupled cross-product recurrence for the remaining cases, seeded by
   standard-library Bessel values near order one.

The imaginary-order dispatcher uses:

1. Hankel expansions for `eta<=0.1*k*r`.
2. Taylor series for `k<=20`.
3. Downward recurrence in double precision for `k<=50` and `eta<max(2*k,100)`.
4. Downward recurrence in long double for `k<=100` and `eta<max(2*k,200)`.
5. Large-order asymptotics otherwise.

These switching thresholds are empirical. Long double precision is platform
dependent. Compilation deliberately avoids `-ffast-math`.
The large-real-order evanescent phase uses `acosh`; the imaginary-order phase
uses `asinh`. A regression check for this distinction is documented in the
validation report.

## Spectrum search and limitations

The standalone search is new: a fixed real-order scan over `[0,k]` with at least
256 points and 32 points per estimated radial oscillation, followed by bisection.
The imaginary scan uses spacing `pi/(32*abs(log(q)))`, followed by bisection.
Root tolerance is approximately `2e-12*(1+order)`. Both scans are bounded at a
million steps; the requested imaginary count is limited to 10000.

Sign scans are practical heuristics, not a certified root-count algorithm.
Approximate dispatch transitions can create spurious sign changes. Near
`lambda=0`, a numerically almost-zero order can be classified on either side;
only an exactly zero floating-point residual is treated as a zero mode.
There is no certified detection of a threshold mode or proof that no roots
were skipped. Very thin annuli, tiny inner radii, high wavenumbers and large
orders require additional validation for the intended application.

Cancellation can destroy relative accuracy near roots. Large real orders can
overflow. Asymptotic series have truncation error and are not uniformly accurate
across turning points. Finite output is not an accuracy guarantee. There is no
arbitrary-precision fallback or error-estimate API. The validation report gives
measured samples, not an accuracy guarantee across the full input domain.

Radial values are evaluated on the entire domain without clipping mode tails.
Transition matrices, scattering calculations and fixed-order wavenumber searches
are outside the present API. To calculate an annular eigenfrequency spectrum at
fixed integer angular order, solve `real_cross(m,k,q)=0` in `k` with a separate
bracketed search; `spectrum` solves a different problem.
