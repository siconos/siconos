#include<iostream>

// ===== Lagrangian DS  =====

// Plugins for Fext, Fint, QNLInertia (vectors), Mass, JacobianQQNLInertia, JacobianVelocityQNLInertia,
// JacobianQFint and JacobianVelocityFint (matrices)

extern "C" void computeFInt(unsigned int *sizeOfq, const double *time, double *q, double *velocity, double *fInt)
{
  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
    fInt[i] = *time * q[i];
}

extern "C" void computeFExt(unsigned int *sizeOfq, const double *time, double *fExt)
{
  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
    fExt[i] = 2 * *time;
}

extern "C" void computeQNLInertia(unsigned int *sizeOfq, double *q, double *velocity, double *Q)
{
  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
    Q[i] = q[i];

}


extern "C" void computeMass(unsigned int *sizeOfq, const double *time, double *q, double *mass)
{
  /* input parameter : sizeOfq (size of the vector q); time ; q (pointer to q vector);
   * output parameter : mass (pointer to mass matrix)
   */

  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = 0; j < n; j++)
      mass[i + j * n] = q[i];
  }

}


extern "C" void computeJacobianQFInt(unsigned int *sizeOfq, const double *time, double *q, double *velocity, double *jacob)
{

  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = 0; j < n; j++)
      jacob[i + j * n] = *time * velocity[i];
  }

}

extern "C" void computeJacobianVelocityFInt(unsigned int *sizeOfq, const double *time, double *q, double *velocity, double *jacob)
{
  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = 0; j < n; j++)
      jacob[i + j * n] = *time * velocity[i];
  }

}

extern "C" void computeJacobianQQNLInertia(unsigned int *sizeOfq, double *q, double *velocity, double *jacob)
{
  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = 0; j < n; j++)
      jacob[i + j * n] = 2 * velocity[i];
  }
}

extern "C" void computeJacobianVelocityQNLInertia(unsigned int *sizeOfq, double *q, double *velocity, double *jacob)
{
  unsigned int n = *sizeOfq;
  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = 0; j < n; j++)
      jacob[i + j * n] = velocity[i];
  }

}

