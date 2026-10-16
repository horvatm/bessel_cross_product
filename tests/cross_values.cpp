// SPDX-License-Identifier: GPL-3.0-or-later
#include "bessel_cross/bessel_cross.hpp"
#include <iomanip>
#include <iostream>
int main() {
 std::cout<<std::setprecision(17);
 double cases[][4]={
 {0,0,20,.5},{0,2.3,20,.5},{0,21,20,.5},{0,50,20,.5},
 {0,5,1000,.5},{0,600,1000,.5},{0,2100,1000,.9},
 {1,.1,20,.5},{1,2,20,.5},{1,10,20,.5},{1,50,20,.5},
 {1,8,30,.5},{1,35,50,.5},{1,80,80,.5},
 {1,150,150,.5},{1,400,100,.5},{1,10,1000,.5},
 {1,200,1000,.5},{1,500,1000,.9}};
 for(auto& c:cases) {
  bool imag=c[0];double v=c[1],k=c[2],r=c[3];
  std::cout<<imag<<','<<v<<','<<k<<','<<r<<','
   <<(imag?bessel_cross::imaginary_cross(v,k,r):bessel_cross::real_cross(v,k,r))<<'\n';
 }
}
