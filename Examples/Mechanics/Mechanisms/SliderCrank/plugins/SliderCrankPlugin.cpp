#include <math.h>
#include <stdio.h>

#include <boost/math/quaternion.hpp>

extern "C" void externalForces(double t, double *f, unsigned int size_z, double *z){
  f[0]=0;
  f[1]=0;
  f[2]=-10;

}

double prevq3=0;
extern "C" void externalMomentum(double t,double *m, unsigned int size_z, double *z){

    m[0]=0;
    m[1]=0;
    m[2]=0;

  // double deltamax = 1e-4;
  // if (fabs(prevq3-q[2])< deltamax){
  //   printf("externalMomentum doing\n");
  //   m[1]=50;
  // }else{
  //   printf("externalMomentum passif\n");
  // }
  // prevq3=q[2];
  // }
  //prevq=
}
extern "C" void externalForcesB1(double t, double *f, unsigned int size_z,double *z){
  f[0]=0;
  f[1]=0;
  f[2]=-9.81*0.038;
  // printf("externalForcesB1 :\n");
  // printf("f[0] = %e\t f[1] = %e\t, f[2]=%e\n",f[0],f[1],f[2]);
}

extern "C" void internalForcesB1(double t, double *q, double *v, double *f, unsigned int size_z,double *z){
  f[0]=0;
  f[1]=0;
  f[2]=1e4*q[2];
  // printf("internalForcesB1 :\n");
  // printf("f[0] = %e\t f[1] = %e\t, f[2]=%e\n",f[0],f[1],f[2]);
}
extern "C" void internalForcesB1_Jacq(double t, double *q, double *v, double *jac, unsigned int size_z,double *z){
  for (int i =0; i < 3; i++)
  {
    for (int j=0; j<7; j++)
      jac[i+j*3]=0.0;
  }
  jac[2+2*3]=1e4;
  // printf("internalForcesB1_Jacq :\n");
  // printf("jac[2+2*3] = %e\n", jac[2+2*3]);
}
extern "C" void externalForcesB2(double t,double *f, unsigned  int size_z,double *z){
  f[0]=0;
  f[1]=0;
  f[2]=-9.81*0.038;

}
extern "C" void externalForcesS(double t,double *f, unsigned int size_z, double *z){
  f[0]=0;
  f[1]=0;
  f[2]=-9.81*0.076;

}
extern "C" void externalForceG(double t,double *f, unsigned  int size_z, double *z){
  // std::cout << "externalForceG"<<std::endl;
  f[0]=0;
  f[1]=0;
  f[2]=-9.81;

}
extern "C" void externalMomentumY(double t,double *m, unsigned  int size_z, double *z){
  // std::cout << "externalMomentumY"<<std::endl;
  m[0]=0;
  m[1]=10;
  m[2]=0;

}
