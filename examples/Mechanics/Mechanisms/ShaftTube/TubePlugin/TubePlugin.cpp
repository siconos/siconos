#include <math.h>
#include <stdio.h>
#include <iostream>

#include <boost/math/quaternion.hpp>

extern "C" void externalForceG(double t, double* f, unsigned int size_z, double *z){
  // std::cout << "externalForceG at time "<< t << std::endl;
  f[0]=0;
  f[1]=0;
  f[2]=0;
  
  if (t > 0.04){
    f[2] = +1000.0;
  }
}
extern "C" void externalMomentumY(double t,double *m, unsigned int size_z, double *z){
  // std::cout << "externalMomentumY"<<std::endl;
  m[0]=5000;
  m[1]=0;
  m[2]=0;
}
