
#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include "Robot.h"
#include <stdio.h>
#include <math.h>

const unsigned int n0 = 21; // Real problem dimension

SICONOS_EXPORT void mass(unsigned int sizeOfq, const double *q, double *Mass, unsigned int sizeZ, double* z)
{
  // compute mass matrix
  Inertia(Mass, q);
}

SICONOS_EXPORT void FGyr(unsigned int sizeOfq, const double *q, const double *velocity, double *FGyr, unsigned int sizeZ, double* z)
{
  // compute mass matrix
  NLEffects(FGyr, q, velocity);
}

SICONOS_EXPORT void jacobianFGyrq(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // compute jacobian matrix
  JacobianQNLEffects(jacob, q, velocity);

}

SICONOS_EXPORT void jacobianVFGyr(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // compute jacobian matrix
  JacobianVNLEffects(jacob, q, velocity);

}

// SICONOS_EXPORT void FInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *fInt, unsigned int sizeZ, double* z)
// {
//   unsigned int i;
//   unsigned int n = sizeOfq;
//   for(i = 0; i<n ; i++)
//     fInt[i] = 0.0;
// }
// SICONOS_EXPORT void jacobianFIntq(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
// {
//   unsigned int i;
//   unsigned int n = sizeOfq * sizeOfq;
//   for(i = 0; i<n ; i++)
//     jacob[i] = 0.0;
// }

// SICONOS_EXPORT void jacobianVFInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
// {
//   unsigned int i;
//   unsigned int n = sizeOfq * sizeOfq;
//   for(i = 0; i<n ; i++)
//     jacob[i] = 0.0;
// }

SICONOS_EXPORT void h0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* z)
{
  unsigned int i;
  double CC[69];

  Contact(CC, q);

  for (i = 0; i < sizeOfY / 3; i++)
  {
    y[i * 3] = CC[sizeOfY / 3 + i];
    y[i * 3 + 1] = CC[i];
    y[i * 3 + 2] = CC[2 * sizeOfY / 3 + i];
  }
}

SICONOS_EXPORT void G0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* z)
{
  unsigned int i, j;
  double CJ[1449];

  ContactJacobian(CJ, q);

  for (i = 0; i < sizeOfq; i++)
  {
    for (j = 0; j < sizeOfY / 3; j++)
    {
      G[i * sizeOfY + 3 * j] = CJ[i * sizeOfY + j + sizeOfY / 3];
      G[i * sizeOfY + 3 * j + 1] = CJ[i * sizeOfY + j];
      G[i * sizeOfY + 3 * j + 2] = CJ[i * sizeOfY + j + 2 * sizeOfY / 3];
    }
  }
}


