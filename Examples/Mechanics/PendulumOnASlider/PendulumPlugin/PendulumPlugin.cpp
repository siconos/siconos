
#ifdef _WIN32
#define SICONOS_EXPORT extern "C" __declspec(dllexport)
#else
#define SICONOS_EXPORT extern "C"
#endif
#include <math.h>
#include <iostream>

#include "RuntimeException.hpp"

using namespace std;

// Inertial parameters
double m1 = 1.0;
double m2 = 1.0;
double l = 1.0;
double a = 0.025;
double d = 0.25;

// force elements
double gravity = 9.81;

SICONOS_EXPORT void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  // columnwise definition of mass matrix
  mass[0] = m1 + m2;
  mass[1] = 0.0;
  mass[2] = cos(q[2]);

  mass[3] = 0.0;
  mass[4] = 1.0;
  mass[5] = 0.0;

  mass[6] = m2*l*cos(q[2]);
  mass[7] = 0.0;
  mass[8] = l;
}

SICONOS_EXPORT void FGyr(unsigned int sizeOfq, const double *q, const double *velocity, double *FGyr, unsigned int sizeZ, double* z)
{
  // nonlinear inertia terms (negative in h according to LagrangianDS)
  FGyr[0] = -m2*l*velocity[2]*velocity[2]*sin(q[2]);
  FGyr[1] = 0.0;
  FGyr[2] = 0.0;
}

SICONOS_EXPORT void jacobianFGyrq(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // Jacobian of nonlinear inertia terms with respect to q (columnwise)
  jacob[0] = 0.0;
  jacob[1] = 0.0;
  jacob[2] = -m2*l*velocity[2]*velocity[2]*cos(q[2]);

  jacob[3] = 0.0;
  jacob[4] = 0.0;
  jacob[5] = 0.0;

  jacob[6] = 0.0;
  jacob[7] = 0.0;
  jacob[8] = 0.0;
}

SICONOS_EXPORT void jacobianFGyrqDot(unsigned int sizeOfq, const double *q, const  double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // Jacobian of nonlinear inertia terms with respect to velocity (columnwise)
  jacob[0] =  0.;
  jacob[1] =  0.;
  jacob[2] =  -2.0*m2*velocity[2]*sin(q[2]);

  jacob[3] = 0.;
  jacob[4] = 0.;
  jacob[5] = 0.;

  jacob[6] = 0.;
  jacob[7] = 0.;
  jacob[8] = 0.;
}

SICONOS_EXPORT void FInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *fInt, unsigned int sizeZ, double* z)
{
  // internal forces (negative in h according to LagrangianDS)
  fInt[0] = 0.0;
  fInt[1] = gravity*q[1];
  fInt[2] = gravity*sin(q[2]);
}

SICONOS_EXPORT void jacobianFIntq(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // Jacobian of internal forces with respect to q (columnwise)
  jacob[0] = 0.0;
  jacob[1] = gravity;
  jacob[2] = 0.;

  jacob[3] = 0.;
  jacob[4] = 0.;
  jacob[5] = 0.;

  jacob[6] = 0.;
  jacob[7] = 0.;
  jacob[8] = gravity*cos(q[2]);
}

SICONOS_EXPORT void jacobianFIntqDot(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeZ, double* z)
{
  // Jacobian of internal forces with respect to velocity (columnwise)
  jacob[0] = 0.;
  jacob[1] = 0.;
  jacob[2] = 0.;

  jacob[3] = 0.;
  jacob[4] = 0.;
  jacob[5] = 0.;

  jacob[6] = 0.;
  jacob[7] = 0.;
  jacob[8] = 0.;
}

SICONOS_EXPORT void g1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* g, unsigned int sizeZ, double* z)
{
  g[0] = q[1];
  g[1] = 0.0;
}

SICONOS_EXPORT void W1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* W, unsigned int sizeZ, double* z)
{
  // Jacobian of g1 (columnwise)
  if (sizeOfY == 1)
  {
    W[0] = 0.0;
    W[1] = 1.0;
    W[2] = 0.0;
  }
  else if (sizeOfY == 2)
  {
    W[0] = 0.0;
    W[1] = 1.0;

    W[2] = 1.0;
    W[3] = 0.0;

    W[4] = 0.0;
    W[5] = 0.0;
  }
  else
    RuntimeException::selfThrow("W1 - not implemented!");
}

SICONOS_EXPORT void g2(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* g, unsigned int sizeZ, double* z)
{
  g[0] = d - q[0] - 0.5*a; // normal
}

SICONOS_EXPORT void W2(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* W, unsigned int sizeZ, double* z)
{
    W[0] = -1.0;
    W[1] = 0.0;
    W[2] = 0.0;
}

SICONOS_EXPORT void g3(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* g, unsigned int sizeZ, double* z)
{
  g[0] = q[0] - 0.5*a; // normal
}

SICONOS_EXPORT void W3(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* W, unsigned int sizeZ, double* z)
{
    W[0] = 1.0;
    W[1] = 0.0;
    W[2] = 0.0;
}
