/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include <stdio.h>

// ==== Dynamical System ====

extern "C" void computeF(double time, unsigned int sizeOfX, const double* x, double* f, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < sizeOfX; i++)
    f[i] = time * x[i];
}

extern "C" void computeJacobianXF(double time, unsigned int sizeOfX, const double* x, double* jacob, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < sizeOfX * sizeOfX; i++)
    jacob[i] = i + 1;
}

// ===== Lagrangian DS  =====

extern "C" void computeMass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    mass[i] = 0;
  mass[0] = 1;
  mass[4] = 2;
  mass[8] = 3;
}

extern "C" void computeFInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *fInt, unsigned int sizeZ, double * z)
{
  for (unsigned int i = 0; i < sizeOfq; ++i)
    fInt[i] = i * q[i];
}

extern "C" void computeFExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeOfZ, double *z)
{
  for (unsigned int i = 0; i < sizeOfq; ++i)
    fExt[i] = i * time;
}

extern "C" void computeNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *NNL, unsigned int sizeOfZ, double *z)
{
  for (unsigned int i = 0; i < sizeOfq; ++i)
    NNL[i] = i * q[i];
}

extern "C" void computeJacobianQFInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}

extern "C" void computeJacobianVelocityFInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}

extern "C" void computeJacobianQNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}

extern "C" void computeJacobianVelocityNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}


//==================  FirstOrderLinearDS ==================

extern "C" void computeB(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double *z)
{
  for (unsigned int i = 0; i < sizeOfB; i++)
    b[i] = time * i ;

}
extern "C" void computeA(double time, unsigned int  sizeOfA, double* A, unsigned int sizeOfZ, double *z)
{
  for (unsigned int j = 0; j < sizeOfA; j++)
  {
    for (unsigned int i = 0; i < sizeOfA; i++)
      A[i + j * sizeOfA] = 4 * (i + 1);
  }
}

//==================  LagrangianScleronomousR ==================
extern "C" void hSclero(unsigned int sizeDS, const double* q, unsigned int sizeY, double* y, unsigned int sizeZ, double* z)
{
  printf("Call of the function 'hSclero' of the test plugin.\n");
}

extern "C" void G0Sclero(unsigned int sizeDS, const double* q, unsigned int sizeY, double* G0, unsigned int sizeZ, double* z)
{
  printf("Call of the function 'G0' of the test plugin.\n");
}

//==================  LagrangianRheonomousR ==================

extern "C" void hRheo(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*)
{
  printf("Call of the function 'hRheo' of the test plugin.\n");
}

extern "C" void G0Rheo(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*)
{
  printf("Call of the function 'G0Rheo' of the test plugin.\n");
}

extern "C" void hDot(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*)
{
  printf("Call of the function 'hDot' of the test plugin.\n");
}

//==================  LagrangianCompliantR ==================

extern "C" void hCompl(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Call of the function 'hCompl' of the test plugin.\n");
}

extern "C" void G0Compl(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Call of the function 'G0Compl' of the test plugin.\n");
}

extern "C" void G1Compl(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Call of the function 'G1Compl' of the test plugin.\n");
}

// ========== FirstOrderType1R ==========

extern "C" void hT1(unsigned int, const double*, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'hT1' of the test plugin.\n");
}

extern "C" void gT1(unsigned int, const double*, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'gT1' of the test plugin.\n");
}

extern "C" void Jh0T1(unsigned int, const double*, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'Jh0T1' of the test plugin.\n");
}

extern "C" void Jg0T1(unsigned int, const double*, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'Jg0T1' of the test plugin.\n");
}

//==================  FirstOrderLinearR ==================

extern "C" void C(double, unsigned int, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'C' of the test plugin.\n");
}

extern "C" void D(double, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'D' of the test plugin.\n");
}

extern "C" void F(double, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'F' of the test plugin.\n");
}

extern "C" void e(double, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'e' of the test plugin.\n");
}

extern "C" void B(double, unsigned int, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'B' of the test plugin.\n");
}
