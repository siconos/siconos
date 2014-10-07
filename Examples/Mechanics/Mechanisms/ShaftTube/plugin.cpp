#include <math.h>
#include <stdio.h>

#include <boost/math/quaternion.hpp>

extern "C" void externalForces(double t, double* q,double *f,double *z){
  f[0]=00;
  f[1]=00;
  f[2]=0.0;
  
}
double prevq3=0;
extern "C" void externalMomentum(double t, double* q,double *m,double *z){

    m[0]=0;
    m[1]=0;
    m[2]=0;
  
  double deltamax = 1e-4;
  if (fabs(prevq3-q[2])< deltamax){
    printf("externalMomentum doing\n");
    m[1]=50;
  }else{
    printf("externalMomentum passif\n");
  }
  prevq3=q[2];
  // }
  //prevq=
}

extern "C" void externalForceG(double t, double* q,double *f,double *z){
  // std::cout << "externalForceG"<<std::endl;
  f[0]=1000;
  f[1]=0;
  f[2]=0;
  
}
extern "C" void externalMomentumY(double t, double* q,double *m,double *z){
  // std::cout << "externalMomentumY"<<std::endl;
  m[0]=0;
  m[1]=0;
  m[2]=0;
  
}
