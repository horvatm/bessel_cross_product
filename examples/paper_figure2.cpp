// SPDX-License-Identifier: GPL-3.0-or-later
// Figure 2 of Horvat & Prosen, J. Phys. A 40 (2007), 6349-6379.
#include "bessel_cross/bessel_cross.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
int main(int argc,char** argv) {
 try {
  if(argc!=2) throw std::invalid_argument("usage: paper_figure2 output.csv");
  constexpr double q=.6;
  constexpr std::size_t branches=26, samples=600;
  const double pi=std::acos(-1.0), kmax=std::sqrt(20000.0);
  std::ofstream out(argv[1]);
  if(!out) throw std::runtime_error("cannot open output CSV");
  out<<std::setprecision(17)<<"n,k,order,eigenvalue,imaginary\n";
  // Exact k=0 limit in logarithmic radius; the cross-product API needs k>0.
  for(std::size_t n=1;n<=branches;++n) {
   double eta=n*pi/std::abs(std::log(q));
   out<<n<<",0,"<<eta<<','<<-eta*eta<<",1\n";
  }
  for(std::size_t j=1;j<=samples;++j) {
   double k=kmax*j/samples;
   std::size_t straight=static_cast<std::size_t>(k*(1-q)/pi);
   auto modes=bessel_cross::spectrum(k,q,branches-straight);
   if(modes.size()<branches) throw std::runtime_error("not enough mode branches");
   for(std::size_t i=0;i<branches;++i) {
    const auto& m=modes[i];
    out<<i+1<<','<<k<<','<<m.order<<','<<m.eigenvalue<<','<<m.imaginary<<'\n';
   }
  }
  if(!out) throw std::runtime_error("failed writing CSV");
  std::cout<<"Wrote "<<branches<<" mode branches at "<<samples+1<<" wavenumbers (q=0.6)\n";
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
