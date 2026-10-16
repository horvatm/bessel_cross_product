// SPDX-License-Identifier: GPL-3.0-or-later
#include "bessel_cross/bessel_cross.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
void check(bool ok) {if(!ok) throw std::runtime_error("numerical check failed");}
int main() {
 using namespace bessel_cross;
 check(real_cross(2,20,1)==0 && imaginary_cross(3,20,1)==0);
 check(imaginary_cross(0,20,.5)==real_cross(0,20,.5));
 check(real_cross(-2,20,.5)==real_cross(2,20,.5));
 // High-precision reference: checks the acosh phase in the real-order asymptotic formula.
 check(std::abs(real_cross(50,20,.5)/1620653003772.2779787-1)<1e-10);
 auto modes=spectrum(20,.5,6);
 check(modes.size()==9);
 double previous=1e100;
 for(auto m:modes) {
  check(m.eigenvalue<previous);previous=m.eigenvalue;
  check(std::abs(radial_value(m,20,.5))<1e-9);
  // Radial ODE: r^2 R''+r R'+(k^2 r^2 - lambda)R=0.
  double r=.73,h=1e-5,y=radial_value(m,20,r);
  double yp=radial_value(m,20,r+h),ym=radial_value(m,20,r-h);
  check(std::abs(r*r*(yp-2*y+ym)/(h*h)+r*(yp-ym)/(2*h)+(400*r*r-m.eigenvalue)*y)<1e-5);
 }
 check(spectrum(1,.5,2).size()==2);
 bool rejected=false;try{spectrum(20,1,2);}catch(const std::invalid_argument&){rejected=true;}check(rejected);
 std::cout<<"Boundary, ordering, mode count, ODE and input checks passed\n";
}
