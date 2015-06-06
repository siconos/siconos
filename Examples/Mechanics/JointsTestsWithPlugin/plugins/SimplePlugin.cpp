#include <math.h>
#include <stdio.h>

#include <boost/math/quaternion.hpp>

extern "C" void externalForces(double t, double *f, unsigned int size_z, double *z){
  f[0]=0;
  f[1]=0;
  f[2]=0;

}

extern "C" void externalMomentum(double t,double *m, unsigned int size_z, double *z){
    m[0]=0;
    m[1]=0;
    m[2]=0;
}

extern "C" void fInt_beam1(double t, double *q, double *v, double *f, unsigned int size_z, double *z){
  f[0]=-1e4*q[0];
  f[1]=0.0;
  f[2]=0.0;-1e4*q[2];

}
extern "C" void jacobianFIntq_beam1(double t, double *q, double *v, double *jac, unsigned int size_z, double *z){
  for (int i =0; i < 3; i++)
  {
    for (int j=0; j<7; j++)
      jac[i+j*3]=0.0;
  }
  jac[0+0*3]=-1e4;
  //jac[2+2*3]=-1e4;
}

extern "C" void jacobianFIntv_beam1(double t, double *q, double *v, double *jac, unsigned int size_z, double *z){
  for (int i =0; i < 3; i++)
  {
    for (int j=0; j< 6; j++)
      jac[i+j*3]=0.0;
  }
}



extern "C" void internalMomentum_beam1(double t,double *m, unsigned int size_z, double *z){
    m[0]=0;
    m[1]=0;
    m[2]=0;
}
