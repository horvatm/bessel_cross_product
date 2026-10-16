// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <cstddef>
#include <vector>
namespace bessel_cross {
// J_v(kr)Y_v(k) - Y_v(kr)J_v(k).
double real_cross(double order, double k, double r);
// Real-valued cross product of order i*order (no additional scaling).
double imaginary_cross(double order, double k, double r);
struct Mode {
 double order;       // nonnegative real order or magnitude of imaginary order
 double eigenvalue;  // +order^2 for real; -order^2 for imaginary
 bool imaginary;
};
// Fixed k, Dirichlet boundaries r=q and r=1. Descending eigenvalues.
std::vector<Mode> spectrum(double k, double q, std::size_t imaginary_count);
double radial_value(const Mode& mode, double k, double r);
}
