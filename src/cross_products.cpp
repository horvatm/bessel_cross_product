// SPDX-License-Identifier: GPL-3.0-or-later
// Bessel cross-product schemes; see docs/NUMERICAL_METHODS.md.
#include "bessel_cross/bessel_cross.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>
namespace bessel_cross { namespace detail {
using namespace std;
using dcomplex = complex<double>;
using dlcomplex = complex<long double>;
constexpr double pi = 3.141592653589793238462643383279502884;
#define mymax(a,b) std::max<double>((a),(b))
inline dcomplex expI(double x) { return {cos(x), sin(x)}; }
inline double bess_JR(double v,double x) { return std::cyl_bessel_j(v,x); }
inline double bess_YR(double v,double x) { return std::cyl_neumann(v,x); }
#define BESS_PREC1 1e-6



dcomplex bess_V(double v, double x) {

  long k, kmax;
  long double xt = x*x, vt = v;
  dlcomplex sum = 1, b = 1;

  if (x != 0) {
    kmax = long(0.5*sqrt(4.0*xt/BESS_PREC1 + vt*vt));

    for (k = 1; k < kmax; k++){

      b *= -(xt/k)/dlcomplex(k,vt);
      sum += b;
    }
  }



  return dcomplex(real(sum),imag(sum));
}
#undef BESS_PREC1



#define func(x) (0.5*(x)*log(1 + (x)*(x)) + c*(x) - d)

long bess_rec_maxindex3(double v, double x, double log_eps){

  double c = log(2*v/x) - 1,
         d = (log_eps - 1 + 0.5*log(1+v*v) - v*atan(v))/v + M_PI/2,
         x0, x1, y0, y1, dx;

  if (d < 0) {
    cerr << "bess_U: Recursion algorithm will fail!\n";
    throw std::runtime_error("numerical scheme failed");
  }

  if (c < 0){
    x0 = sqrt(exp(-2*c)-1);  y0 = func(x0);
    x1 = x0 + 1/v;
  } else {
    x0 = d/(c + 1);  y0 = func(x0);
    x1 = x0 - 1/v;
  }


  dx  = x1 - x0;

  do {
    y1 = func(x1);

    dx *= - y1/(y1-y0);
    y0 = y1; x0 = x1;
    x1 += dx;

  } while (abs(dx) > 1e-6);



  return mymax(5, long(v*x1));
}
#undef func

template <class T> dcomplex bess_U(double v, double x) {

  long N, M, k;

  T eps;

  if (sizeof(T) == sizeof(double)) {
    eps = 1e-20*exp(-1.5*v);

    N = bess_rec_maxindex3(v, x, -log(eps));

  } else if (sizeof(T) == sizeof(long double)) {
    eps = 1e-50*exp(-1.5*v);

    N = bess_rec_maxindex3(v, x, -log(eps));

  } else {
    cerr << "bess_U: Don't know this data type of size" << sizeof(T) << '\n';
    throw std::runtime_error("numerical scheme failed");
  }


  N = N +  (N % 2);
  M = N >> 1;



  complex<T> *J = new complex<T>  [M], a,
             *p, *b, *end = &J[M-1], tmp = T(eps);

  T xt = 1.0/x, vt = T(v)*xt;

  p = b = end;
  (*p) = complex<T>((N-1)*xt,vt)*tmp;

  k = N - 2;



  do {
    b--;
    tmp = complex<T>(k*xt,vt)*(*p) - tmp;  k--;
    (*b) = complex<T>(k*xt,vt)*tmp - (*p); k--;
    p = b;

  } while (p != J);


  vt = T(v);


  a = complex<T>(1.0)/((*J)*complex<T>(0,vt));


  for (k = 1, p = J + 1; k < M; k++, p++){
    (*p) *= (a *= complex<T>( T(k - 1)/k, vt/k))*complex<T>(2*k,vt);

  }

  complex<T> sum = complex<T>(1.0);

  for (p = end; p != J; p--) sum += (*p);

  delete [] J;

  sum = complex<T>(1.0)/sum;





  return dcomplex(real(sum), imag(sum));
}

long double bess_Pi(double y, double z){

  if (y > z) {
   cerr << "bess_Pi: Is not applicable\n";
   throw std::runtime_error("numerical scheme failed");
  }

  long k = 2, l = 2,
       kmin = long(z + sqrt(z*z - 0.25*y*y));

  long double sum, sum0, a, b,
              lr = -64*z*z, mil = 4*y*y;


  sum = sum0 = b = 1;
  do {
    sum0 = sum;
    b *= a = (mil + (long double)(l-1)*(l-1))*
             (mil + (long double)(l+1)*(l+1))/(lr*k*(k-1));
    sum += b;



    k += 2; l += 4;

  } while (((k <=  kmin+1) || (abs(a) <= 1)) && (abs(b/sum) > 1e-16));

 return sum0;
}

long double bess_Qi(double y, double z){

  if (y > z) {
   cerr << "bess_Qi: Is not applicable\n";
   throw std::runtime_error("numerical scheme failed");
  }

  long k = 3, l = 4,
       kmin = long(z + sqrt(z*z - 0.25*y*y));

  long double sum, sum0, a, b,
              lr = -64*z*z, mil = 4*y*y;



  sum = b = sum0 = -(mil+1)/(8*z);

  do {
    sum0 = sum;
    b *= a = (mil + (long double)(l-1)*(l-1))*
             (mil + (long double)(l+1)*(l+1))/(lr*k*(k-1));
    sum += b;



    k += 2; l += 4;
  } while (((k <=  kmin+1) || (abs(a) <= 1)) && (abs(b/sum) > 1e-16));

  return sum0;
}

int bess_G_a_num = 8;
inline long double bess_G_a(int n, long double t){

  switch (n) {
    case 0: return 1.0;
    case 1: return (3 - t*5)/24.0;
    case 2: return (81 + t*(-462 + 385*t))/1152.0;
    case 3: return (30375 + t*(-369603 + t*(765765- t*425425)))/414720.0;
    case 4: return (4465125 + t*(-94121676 + t*(349922430 +
                     t*(-446185740 + t*185910725))))/39813120.0;
    case 5: return (1519035525. +t*(- 49286948607. + t*(284499769554. +
                    t*(- 614135872350. + t*(566098157625. - t*188699385875.)))))/
                    6688604160.0;
    case 6: return (2757049477875. + t*(- 127577298354750. +
                    t*(1050760774457901. + t*(-3369032068261860. +
                    t*(5104696716244125. + t*(-3685299006138750. +
                    t*1023694168371875.))))))/4815794995200.0;
    case 7: return (199689155040375. + t*(-12493049053044375. +
                    t*(138799253740521843. + t*(-613221795981706275. +
                    t*(1347119637570231525. + t*(-1570320948552481125. +
                    t*(931766432052080625. - t*221849150488590625.)))))))/
                     115579079884800.0;

  }
  return 0;
}



inline dcomplex bess_Gi(double t, double x){

  dlcomplex sum = 1, a = 1, tt = dlcomplex(0,-1.0/x);

  for (long m = 1; m < bess_G_a_num; m++){
    a *= tt;
    sum += a*bess_G_a(m,t);
  }

  return dcomplex(real(sum),imag(sum));
}


inline double bess_Ga(double t, double x){

  long double sum = 1, a = 1, tt = 1.0/x;

  for (long m = 1; m < bess_G_a_num; m++){
    a *= tt;
    sum += a*bess_G_a(m,t);
  }



  return sum;
}

inline double bess_Gb(double t, double x){

  long double sum = 1, a = 1, tt = -1.0/x;

  for (long m = 1; m < bess_G_a_num; m++){
    a *= tt;
    sum += a*bess_G_a(m,t);
  }


  return sum;
}



inline void bess_give_pqr_at_v(double v, double a, double b,
                               long double &p, long double &q, long double &r){

  v--;

  int n = 3;
  double *ja = new double [3],
         *ya = new double [3],
         *jb = new double [3],
         *yb = new double [3];

  for (int i=0; i<n; ++i) ja[i] = bess_JR(v+i, a);
  for (int i=0; i<n; ++i) ya[i] = bess_YR(v+i, a);
  for (int i=0; i<n; ++i) jb[i] = bess_JR(v+i, b);
  for (int i=0; i<n; ++i) yb[i] = bess_YR(v+i, b);


  p = ja[1]*yb[1] - jb[1]*ya[1];
  q = 0.5*( ja[1]*(yb[0]-yb[2]) - (jb[0]-jb[2])*ya[1]);
  r = 0.5*( (ja[0]-ja[2])*yb[1] - jb[1]*(ya[0]-ya[2]));



  delete [] ja;
  delete [] ya;
  delete [] jb;
  delete [] yb;

}



inline double bess_crossR_rec(double v, double k, double r){

  long i, n = long(v) - 1;

  long double p_1, p0, p1,
         q_1, q0, q1,
         r_1, r0, r1,
         x = 1.0/(k*r), y = 1.0/k;

  v -= n;
  bess_give_pqr_at_v(v, k*r, k, p_1, q_1, r_1);
  v += 1;
  bess_give_pqr_at_v(v, k*r, k, p0 , q0 , r0);

  for (i = 1; i < n; i++, v += 1){

    p1 = p_1 - 2*v*(x*q0 + y*r0);
    q1 = v*x*p0 - (v+1)*y*p1 - r0;
    r1 = v*y*p0 - (v+1)*x*p1 - q0;

    p_1 = p0; p0 = p1;
    q_1 = q0; q0 = q1;
    r_1 = r0; r0 = r1;
  }

  return p0;
}



long double bess_Pr(double mi, double r){

  long k = 2, l = 2, kmin = long (sqrt(mi));
  long double sum = 1, sum0 = sum, a, b = 1, lr = -r*r, mil = mi;

  do {
    sum0 = sum;
    b *= a = (mil - (long double)(l-1)*(l-1))*
             (mil - (long double)(l+1)*(l+1))/(lr*k*(k-1));
    sum += b;



    k += 2; l += 4;

  } while (((k <=  kmin+1) || (abs(a) <= 1)) && (abs(b/sum) > 1e-16));

 return sum0;
}

long double bess_Qr(double mi, double r){

  long k = 3, l = 4, kmin = long(sqrt(mi));
  long double sum = (mi-1)/r, sum0 = sum, a, b = sum, lr = -r*r, mil = mi;

  do {
    sum0 = sum;
    b *= a = (mil - (long double)(l-1)*(l-1))*
             (mil - (long double)(l+1)*(l+1))/(lr*k*(k-1));
    sum += b;



    k += 2; l += 4;
  } while (((k <=  kmin+1) || (abs(a) <= 1)) && (abs(b/sum) > 1e-16));

  return sum0;
}

inline double bess_crossR_han(double v, double k, double r){

  long double mi = 4*v*v, x = k*(1 - r),
         A = bess_Pr(mi, 8.0*k*r), B = bess_Qr(mi, 8.0*k*r),
         C = bess_Pr(mi, 8.0*k),   D = bess_Qr(mi, 8.0*k);

  return M_2_PI*((A*C + B*D)*sin(x) + (A*D - B*C)*cos(x))/(k*sqrt(r));
}


inline double bess_crossR_asim(double v, double k, double r){

  long double v2 = v*v, k2 = k*k, kr2 = k2*r*r,
         x1 = sqrt(v2 - k2), x2 = sqrt(v2 - kr2),
         t1 = 1.0/(1.0 - k2/v2), t2 = 1.0/(1.0 - kr2/v2),
         phi = x2 - x1 - v*(acosh(v/(k*r)) - acosh(v/k));

  return (-exp(phi)*bess_Ga(t2, x2)*bess_Gb(t1, x1) +
          exp(-phi)*bess_Ga(t1, x1)*bess_Gb(t2, x2))/(M_PI*sqrt(x1*x2));
}




thread_local bool bess_crossR_verbose = false;
thread_local int bess_crossR_choice = -1;

inline double bess_crossR(double v, double k, double r){

  if (r == 1) return 0;

  v = abs(v);




  if (k < 1000 && v < k){

    if (bess_crossR_verbose){
      cerr << "bess_crossR 2 k*r=" << k*r << " k=" << k << " v=" << v << '\n';
      cerr.flush();
    }

    bess_crossR_choice = 3;
    return  bess_JR(v, k*r)*bess_YR(v, k) - bess_JR(v, k)*bess_YR(v, k*r);
  }

  if (bess_crossR_verbose){
    cerr << "bess_crossR 3 k*r=" << k*r << " k=" << k << " v=" << v << '\n';
    cerr.flush();
  }

  if (v > 2*k ) {
    bess_crossR_choice = 4;
    return bess_crossR_asim(v, k, r);
  }

  if (v <= 0.1*k*r){
    bess_crossR_choice = 5;
    return bess_crossR_han(v, k, r);
  }

  bess_crossR_choice = 6;
  return bess_crossR_rec(v,k,r);
}



inline double bess_crossI_tay(double v, double k, double r){
  return -M_2_PI*imag(expI(log(r)*v)*
          bess_V(v,0.5*k*r)*conj(bess_V(v,0.5*k)))/v;
}



template <class T>
inline double bess_crossI_rec(double v, double k, double r){
  return -M_2_PI*imag(expI(log(r)*v)*
          bess_U<T>(v,0.5*k*r)*conj(bess_U<T>(v,0.5*k)))/v;
}

inline double bess_crossI_han(double v, double k, double r){
  long double x = k*(1 - r),
         A = bess_Pi(v, k*r), B = bess_Qi(v, k*r),
         C = bess_Pi(v, k),   D = bess_Qi(v, k);

  return M_2_PI*((A*C + B*D)*sin(x) + (A*D - B*C)*cos(x))/(k*sqrt(r));
}

inline double bess_crossI_asim(double v, double k, double r){
  double v2 = v*v, k2 = k*k, kr2 = k2*r*r,
         x1 = sqrt(v2 + k2), x2 = sqrt(v2 + kr2),
         t1 = 1.0/(1.0 + k2/v2), t2 = 1.0/(1.0 + kr2/v2);


  return M_2_PI/((1-exp(-2*pi*v))*sqrt(x1*x2))*
         imag( expI(x1 - x2 + v*(asinh(v/(k*r)) - asinh(v/k)) )*
               bess_Gi(t1, x1)*conj(bess_Gi(t2, x2)) );

}


thread_local int bess_crossI_choice = -1;

inline double bess_crossI(double v, double k, double r){

  if (r == 1) return 0;


  v = abs(v);


  if (v <= 0.1*k*r) {
    bess_crossI_choice = 1;
    return bess_crossI_han(v, k, r);
  }


  if (k <= 20 )  {
    bess_crossI_choice = 2;
    return bess_crossI_tay(v, k, r);
  }

  if (k <= 50 && v < mymax(2*k, 100)) {
    bess_crossI_choice = 3;
    return bess_crossI_rec<double>(v, k, r);

   }

  if (k <= 100 && v < mymax(2*k, 200)) {
    bess_crossI_choice = 4;
    return bess_crossI_rec<long double>(v, k, r);
   }

  bess_crossI_choice = 5;
  return bess_crossI_asim(v, k, r);
}


#undef mymax
} // namespace detail
namespace {
void validate(double order, double k, double r) {
 if (!std::isfinite(order) || !std::isfinite(k) || !std::isfinite(r) || k<=0 || r<=0 || r>1)
   throw std::invalid_argument("require finite order, k>0 and 0<r<=1");
}
}
double real_cross(double order,double k,double r) {
 validate(order,k,r);
 double z=detail::bess_crossR(order,k,r);
 if (!std::isfinite(z)) throw std::runtime_error("real cross product overflow or numerical failure");
 return z;
}
double imaginary_cross(double order,double k,double r) {
 validate(order,k,r);
 if (order==0) return real_cross(0,k,r);
 double z=detail::bess_crossI(order,k,r);
 if (!std::isfinite(z)) throw std::runtime_error("imaginary cross product overflow or numerical failure");
 return z;
}
} // namespace bessel_cross
