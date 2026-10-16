// SPDX-License-Identifier: GPL-3.0-or-later
#include "bessel_cross/bessel_cross.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace bessel_cross {
namespace {
template<class F> double bisect(F f,double a,double b,double fa) {
 for(int i=0;i<100;++i) {
  double m=a+(b-a)/2, fm=f(m);
  if(fm==0 || b-a<2e-12*(1+std::abs(m))) return m;
  if(std::signbit(fa)!=std::signbit(fm)) b=m;
  else {a=m;fa=fm;}
 }
 throw std::runtime_error("bisection did not converge");
}
}
std::vector<Mode> spectrum(double k,double q,std::size_t count) {
 if(!std::isfinite(k)||!std::isfinite(q)||k<=0||q<=0||q>=1||count>10000)
  throw std::invalid_argument("require k>0, 0<q<1, imaginary_count<=10000");
 const double pi=std::acos(-1.0);
 std::vector<Mode> result;
 // In log radius, the potential is k^2 exp(2t), so positive orders are < k.
 // Resolve at least 32 samples per estimated radial oscillation.
 double estimate=std::max(256.0,32*k*(1-q)/pi);
 if(estimate>1000000) throw std::invalid_argument("spectrum scan too large");
 std::size_t steps=static_cast<std::size_t>(std::ceil(estimate));
 auto real=[&](double v){return real_cross(v,k,q);};
 double a=0,fa=real(a);
 if(fa==0) result.push_back({0,0,false});
 for(std::size_t i=1;i<=steps;++i) {
  double b=k*i/steps,fb=real(b);
  if(fb==0 || (fa!=0 && std::signbit(fa)!=std::signbit(fb))) {
   double v=fb==0 ? b : bisect(real,a,b,fa);
   result.push_back({v,v*v,false});
  }
  a=b;fa=fb;
 }
 std::reverse(result.begin(),result.end());
 auto imag=[&](double v){return imaginary_cross(v,k,q);};
 double step=pi/(-std::log(q))/32;
 a=0;fa=imag(0);
 std::size_t found=0;
 for(std::size_t i=1;found<count && i<=1000000;++i) {
  double b=i*step,fb=imag(b);
  if(fb==0 || (fa!=0 && std::signbit(fa)!=std::signbit(fb))) {
   double v=fb==0 ? b : bisect(imag,a,b,fa);
   result.push_back({v,-v*v,true}); ++found;
  }
  a=b;fa=fb;
 }
 if(found!=count) throw std::runtime_error("imaginary-order scan limit exceeded");
 return result;
}
double radial_value(const Mode& mode,double k,double r) {
 return mode.imaginary ? imaginary_cross(mode.order,k,r) : real_cross(mode.order,k,r);
}
}
