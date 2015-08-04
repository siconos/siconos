#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <stdio.h>
#include "donnee.h"

#undef restrict
#define restrict __restrict

// forces extérieures appliquées sur le yoyo dans la phase contrainte
SICONOS_EXPORT void force_ext(double time, unsigned int sizeOfq, double* restrict fExt, unsigned int sizeZ, double* restrict  z)
{
  fExt[0] = -m * r * g;
  fExt[1] = 0;
  // fExt[2] = accelerationmain(5,A,Cy,time);
  fExt[2] = 0;
}

// forces extérieures appliquées sur le yoyo dans la phase libre
SICONOS_EXPORT void force_extf(double time, unsigned int sizeOfq, double* restrict fExt, unsigned int sizeZ, double* restrict  z)
{
  fExt[0] = 0;
  fExt[1] = -m * g;
  //fExt[2] = accelerationmain (5,A,Cy,time);
  fExt[2] = 0;
}

// forces intérieures appliquées sur le yoyo dans la phase contrainte
SICONOS_EXPORT void F_int(double time, unsigned int sizeOfq,  double* restrict q,  double* restrict velocity, double* restrict fInt, unsigned int sizeZ, double* restrict  z)
{
  fInt[0] =  r * epsilon * (velocity[0]);
  fInt[1] =  0;
  int i = 0;
  while (temps[i] < time) i++ ;
  if (velocity[0] < 0 && q[0] < thetaset(Som, i)) fInt[2] = -w * g   ;
  else  fInt[2] = c1 * velocity[2] + c2 * q[2] ;
  //fInt[2] =0;
}


SICONOS_EXPORT void jacobianFIntq(double time, unsigned int sizeOfq,  double* restrict q,  double* restrict velocity, double* restrict jacob, unsigned int sizeZ, double* restrict  z)
{
  jacob[0] =  0;
  jacob[1] =  0;
  int i = 0;
  while (temps[i] < time) i++ ;
  if (velocity[0] < 0 && q[0] < thetaset(Som, i)) jacob[2] = 0  ;
  else  jacob[2] = c2 ;
  //jacob[2] =0;
}

SICONOS_EXPORT void jacobianVFInt(double time, unsigned int sizeOfq,  double* restrict q,  double* restrict velocity, double* restrict jacob, unsigned int sizeZ, double* restrict  z)
{
  jacob[0] = r * epsilon;
  jacob[1] =   0;
  int i = 0;
  while (temps[i] < time) i++ ;
  if (velocity[0] < 0 && q[0] < thetaset(Som, i)) jacob[2] = 0  ;
  else  jacob[2] = c1 ;
  //  jacob[2] =0;
}

SICONOS_EXPORT void h1(unsigned int sizeDS,  double* restrict  q, double time, unsigned int sizeY, double* restrict  y, unsigned int sizeOfZ, double* restrict  z)
{
  y[0] = q[1] - r * q[0] + L - q[2];
}

SICONOS_EXPORT void G10(unsigned int sizeDS,  double* restrict  q, double time, unsigned int  sizeY, double* restrict  G, unsigned int sizeOfZ, double* restrict  z)
{
  G[0] = -r;
  G[1] = 1;
  G[2] = -1;
}

SICONOS_EXPORT void G11(unsigned int sizeDS,  double* restrict  q, double time, unsigned int  sizeY, double* restrict  G, unsigned int sizeOfZ, double* restrict  z)
{
  G[0] = 0;
  //G[0]= -vitessemain(5,A,Cy,time);
}


// forces intérieures appliquées sur le yoyo dans la phase libre
SICONOS_EXPORT void F_intf(double time, unsigned int sizeOfq,  double* restrict q,  double* restrict velocity, double* restrict fInt, unsigned int sizeZ, double* restrict  z)
{
  fInt[0] =  0;
  fInt[1] =  0;
  int i = 0;
  while (temps[i] < time) i++ ;
  if (velocity[0] < 0 && q[0] < thetaset(Som, i)) fInt[2] = -w * g   ;
  else  fInt[2] = c1 * velocity[2] + c2 * q[2] ;
  //fInt[2] =0;
}


SICONOS_EXPORT void jacobianFIntqf(double time, unsigned int sizeOfq,  double* restrict q,  double* restrict velocity, double* restrict jacob, unsigned int sizeZ, double* restrict  z)
{
  jacob[0] =  0;
  jacob[1] =  0;
  int i = 0;
  while (temps[i] < time) i++ ;
  if (velocity[0] < 0 && q[0] < thetaset(Som, i)) jacob[2] = 0  ;
  else  jacob[2] = c2 ;
  //jacob[2] =0;
}

SICONOS_EXPORT void jacobianVFIntf(double time, unsigned int sizeOfq,  double* restrict q,  double* restrict velocity, double* restrict jacob, unsigned int sizeZ, double* restrict  z)
{
  jacob[0] = 0;
  jacob[1] =   0;
  int i = 0;
  while (temps[i] < time) i++ ;
  if (velocity[0] < 0 && q[0] < thetaset(Som, i)) jacob[2] = 0  ;
  else  jacob[2] = c1 ;
  //jacob[2] =0;
}
























