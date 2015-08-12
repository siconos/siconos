#include <math.h>
#include <stdio.h>

#include <boost/math/quaternion.hpp>

extern "C" void externalForces(double t, double *f, unsigned int size_z, double *z){
  f[0]=0;
  f[1]=0;
  f[2]=0;

}

extern "C" void externalMoment(double t,double *m, unsigned int size_z, double *z){
    m[0]=0;
    m[1]=0;
    m[2]=0;
}

extern "C" void fInt_beam1(double t, double *q, double *v, double *f, unsigned int size_z, double *z){
  // printf("fInt_beam1\n");
  // printf("q[0] = %e\t q[1] = %e\t, q[2]=%e\n",q[0],q[1],q[2]);
  f[0]=1e4*q[0];
  f[1]=0.0;
  f[2]=1e4*q[2];
  // printf("f[0] = %e\t f[1] = %e\t, f[2]=%e\n",f[0],f[1],f[2]);

}
extern "C" void jacobianFIntq_beam1(double t, double *q, double *v, double *jac, unsigned int size_z, double *z){
  for (int i =0; i < 3; i++)
  {
    for (int j=0; j<7; j++)
      jac[i+j*3]=0.0;
  }
  jac[0+0*3]=1e4;
  jac[2+2*3]=1e4;
}

extern "C" void jacobianFIntv_beam1(double t, double *q, double *v, double *jac, unsigned int size_z, double *z){
  for (int i =0; i < 3; i++)
  {
    for (int j=0; j< 6; j++)
      jac[i+j*3]=0.0;
  }
}


extern "C" void mInt_beam1(double t, double *q, double *v, double *m, unsigned int size_z, double *z){
  // printf("mInt_beam1 :\n");

  // Simple torsional spring around y axis
  // printf("q[0] = %e\n", q[0]);
  // printf("q[1] = %e\n", q[1]);
  // printf("q[2] = %e\n", q[2]);
  // printf("q[3] = %e\n", q[3]);
  // printf("q[4] = %e\n", q[4]);
  // printf("q[5] = %e\n", q[5]);
  // printf("q[6] = %e\n", q[6]);

  double angle = 2*asin(q[5]);
  // printf("angle = %e\n", angle);
  m[0]=0.0;
  m[1]=1e3*(angle-1.0);
  m[2]=0.0;
  // printf("m[0] = %e\t m[1] = %e\t, m[2]=%e\n",m[0],m[1],m[2]);
}

extern "C" void jacobianMIntq_beam1(double t, double *q, double *v, double *jac, unsigned int size_z, double *z){
  // printf("jacobianMIntq_beam1:\n ");
  for (int i =0; i < 3; i++)
  {
    for (int j=0; j<7; j++)
      jac[i+j*3]=0.0;
  }
  // printf("q[0] = %e\n", q[0]);
  // printf("q[1] = %e\n", q[1]);
  // printf("q[2] = %e\n", q[2]);
  // printf("q[3] = %e\n", q[3]);
  // printf("q[4] = %e\n", q[4]);
  // printf("q[5] = %e\n", q[5]);
  // printf("q[6] = %e\n", q[6]);

  double angle = 2*asin(q[5]);

  // printf("angle = %e\n", angle);
  jac[1+5*3]=1e3 * 2.0 / sqrt(1 - q[5]*q[5]) ;
  //printf("jac[3+3*3] = %e\n", jac[3+3*3]);
  // printf("[(");
  // for (int i =0; i < 3; i++)
  // {
  //   for (int j=0; j<7; j++) printf("%e,\t",jac[i+j*3]);
  //   printf("\n");
  // }
  // printf(")]\n");
  // // Computation by finite difference
  // double epsilon = 1e-8;
  // double vector[7];
  // double m1[3],m2[3];

  // for (int j =0; j < 7; j++)
  // {
  //   for (int k =0; k < 7; k++) vector[k] =q[k];
  //   vector[j] += epsilon;
  //   mInt_beam1(t, q, v, m1,  size_z, z);
  //   mInt_beam1(t, vector, v, m2,  size_z, z);
  //   jac[0+j*3]=(m2[0]-m1[0])/epsilon ;
  //   jac[1+j*3]=(m2[1]-m1[1])/epsilon ;
  //   jac[2+j*3]=(m2[2]-m1[2])/epsilon ;
  // }

  // for (int i =0; i < 3; i++)
  // {
  //   for (int j=0; j<7; j++) printf("%e\t",jac[i+j*3]);
  //   printf("\n");
  // }

}




extern "C" void jacobianMIntv_beam1(double t, double *q, double *v, double *jac, unsigned int size_z, double *z){
  for (int i =0; i < 3; i++)
  {
    for (int j=0; j< 6; j++)
      jac[i+j*3]=0.0;
  }
}
