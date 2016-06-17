#include <math.h>
#include <stdio.h>

#include <boost/math/quaternion.hpp>

extern "C" void externalForces(double t, double *f, unsigned int size_z, double *z){
  f[0]=0;
  f[1]=0;
  f[2]=-10;
}
extern "C" void externalMoment(double t,double *m, unsigned int size_z, double *z){
    m[0]=0;
    m[1]=0;
    m[2]=0;
}
extern "C" void externalForcesB1(double t, double *f, unsigned int size_z,double *z){
  f[0]=0;
  f[1]=0;
  f[2]=-9.81*0.038;
  // printf("externalForcesB1 :\n");
  // printf("f[0] = %e\t f[1] = %e\t, f[2]=%e\n",f[0],f[1],f[2]);
}
extern "C" void internalForcesB1(double t, double *q, double *v, double *f, unsigned int size_z,double *z){
  // Simple spring in z direction.
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

extern "C" void internalMomentsB1(double t, double *q, double *v, double *m, unsigned int size_z,double *z){
  //  printf("internalMomentsB1 :\n");
  // Simple torsional spring around y axis
  // printf("q[3] = %e\n", q[3]);
  // printf("q[4] = %e\n", q[4]);
  // printf("q[5] = %e\n", q[5]);
  // printf("q[6] = %e\n", q[6]);

  double angle = 2*asin(q[5]);
  // printf("angle = %e\n", angle);
  m[0]=0.0;
  m[1]=1e3*(angle);
  m[2]=0.0;
  // printf("m[0] = %e\t m[1] = %e\t, m[2]=%e\n",m[0],m[1],m[2]);
}

extern "C" void internalMomentsB1_Jacq(double t, double *q, double *v, double *jac, unsigned int size_z,double *z){
  //printf("internalMomentsB1_Jacq :\n");
  for (int i =0; i < 3; i++)
  {
    for (int j=0; j<7; j++)
      jac[i+j*3]=0.0;
  }
  // printf("q[3] = %e\n", q[3]);
  // printf("q[4] = %e\n", q[4]);
  // printf("q[5] = %e\n", q[5]);
  // printf("q[6] = %e\n", q[6]);

  double angle = 2*asin(q[5]);

  // printf("angle = %e\n", angle);
  jac[1+5*3]=1e3 * 2.0 / sqrt(1 - q[5]*q[5]) ;
  // printf("jac[3+3*3] = %e\n", jac[3+3*3]);
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

extern "C" void prescribedvelocityB1(double time, unsigned int sizeofprescribedvelocity, double *pv)
{
  /* the plugin implements v(t) = C + A cos(omega *t) */

  double C = -150.0 ;
  double omega = M_PI / 2.0;
  double A = 10.0;

  //pv[0] =  A * cos(omega * time*100.0);
  pv[0] =  C;
  //printf("prescribed velocity = %e\n", pv[0]);
}
