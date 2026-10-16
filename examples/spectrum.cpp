// SPDX-License-Identifier: GPL-3.0-or-later
#include "bessel_cross/bessel_cross.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <stdexcept>
int main(int argc,char** argv) {
 try {
  if(argc!=1 && argc!=5) throw std::invalid_argument("usage: spectrum [k q imaginary_count output_prefix]");
  double k=argc==5?std::stod(argv[1]):20, q=argc==5?std::stod(argv[2]):0.5;
  std::size_t count=argc==5?std::stoul(argv[3]):6;
  std::string prefix=argc==5?argv[4]:"spectrum";
  auto modes=bessel_cross::spectrum(k,q,count);
  std::ofstream out(prefix+".csv"),radial(prefix+"_radial.csv");
  if(!out||!radial) throw std::runtime_error("cannot open output files");
  out<<std::setprecision(17)<<"index,k,q,order,eigenvalue,imaginary\n";
  radial<<std::setprecision(17)<<"index,r,value\n";
  for(std::size_t i=0;i<modes.size();++i) {
   auto m=modes[i];
   out<<i<<','<<k<<','<<q<<','<<m.order<<','<<m.eigenvalue<<','<<m.imaginary<<'\n';
   for(int j=0;j<=400;++j) {
    double r=q+(1-q)*j/400;
    radial<<i<<','<<r<<','<<bessel_cross::radial_value(m,k,r)<<'\n';
   }
  }
  if(!out||!radial) throw std::runtime_error("failed writing output");
  std::cout<<"Wrote "<<modes.size()<<" modes to "<<prefix<<".csv and "<<prefix<<"_radial.csv\n";
 } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
