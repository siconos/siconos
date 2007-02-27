
#include "Robot.h"
#include <stdio.h>
#include <math.h>

const unsigned int n0 = 21; // Real problem dimension

extern "C" void mass(unsigned int sizeOfq, const double *q, double *Mass, double* param)
{
  // compute mass matrix
  Inertia(Mass, q);
}

extern "C" void NNL(unsigned int sizeOfq, const double *q, const double *velocity, double *NNL, double* param)
{
  // compute mass matrix
  NLEffects(NNL, q, velocity);
}

extern "C" void jacobianQNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, double* param)
{
  // compute jacobian matrix
  JacobianQNLEffects(jacob, q, velocity);
}

extern "C" void jacobianVNNL(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacob, double* param)
{
  // compute jacobian matrix
  JacobianVNLEffects(jacob, q, velocity);
}

// extern "C" void FInt(unsigned int sizeOfq, double time, const double *q, const double *velocity, double *fInt, double* param)
// {
//   unsigned int i;
//   unsigned int n = sizeOfq;
//   for(i = 0; i<n ; i++)
//     fInt[i] = 0.0;
// }
// extern "C" void jacobianQFInt(unsigned int sizeOfq, double time, const double *q, const double *velocity, double *jacob, double* param)
// {
//   unsigned int i;
//   unsigned int n = sizeOfq * sizeOfq;
//   for(i = 0; i<n ; i++)
//     jacob[i] = 0.0;
// }

// extern "C" void jacobianVFInt(unsigned int sizeOfq, double time, const double *q, const double *velocity, double *jacob, double* param)
// {
//   unsigned int i;
//   unsigned int n = sizeOfq * sizeOfq;
//   for(i = 0; i<n ; i++)
//     jacob[i] = 0.0;
// }

extern "C" void h0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, double* param)
{
  unsigned int i;
  double CC[69];

  Contact(CC, q);

  for (i = 0; i < sizeOfY; i++)
    y[i] = CC[23 + i];
}

extern "C" void G0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, double* param)
{
  unsigned int i, j;
  double CJ[1449];

  ContactJacobian(CJ, q);

  for (i = 0; i < sizeOfq; i++)
  {
    for (j = 0; j < sizeOfY; j++)
      G[i * sizeOfY + j] = CJ[i * 3 * sizeOfY + j + 23];
  }
}


