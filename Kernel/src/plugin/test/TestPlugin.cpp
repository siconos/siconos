/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include <cstdio>

#if defined(_MSC_VER)
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

// ==== Dynamical System ====

extern "C" DLLEXPORT void computef(double time, unsigned int sizeOfX, const double* x, double* f, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < sizeOfX; i++)
    f[i] = time * x[i];
}

extern "C" DLLEXPORT void computeJacobianfx(double time, unsigned int sizeOfX, const double* x, double* jacob, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < sizeOfX * sizeOfX; i++)
    jacob[i] = i + 1;
}

// ===== Lagrangian DS  =====

extern "C" DLLEXPORT void computeMass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    mass[i] = 0;
  mass[0] = 1;
  mass[4] = 2;
  mass[8] = 3;
}

extern "C" DLLEXPORT void computeFInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *fInt, unsigned int sizeZ, double * z)
{
  for (unsigned int i = 0; i < sizeOfq; ++i)
    fInt[i] = i * q[i];
}

extern "C" DLLEXPORT void computeFExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeOfZ, double *z)
{
  for (unsigned int i = 0; i < sizeOfq; ++i)
    fExt[i] = i * time;
}

extern "C" DLLEXPORT void computeNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *NNL, unsigned int sizeOfZ, double *z)
{
  for (unsigned int i = 0; i < sizeOfq; ++i)
    NNL[i] = i * q[i];
}

extern "C" DLLEXPORT void computeJacobianFIntq(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}

extern "C" DLLEXPORT void computeJacobianFintVelocity(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}

extern "C" DLLEXPORT void computeJacobianNNLq(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}

extern "C" DLLEXPORT void computeJacobianNNLVelocity(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}


//==================  FirstOrderLinearDS ==================

extern "C" DLLEXPORT void computeb(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double *z)
{
  for (unsigned int i = 0; i < sizeOfB; i++)
    b[i] = time * i ;

}
extern "C" DLLEXPORT void computeA(double time, unsigned int rowA, unsigned int colA, double* A, unsigned int sizeOfZ, double *z)
{
  for (unsigned int j = 0; j < rowA; j++)
  {
    for (unsigned int i = 0; i < rowA; i++)
      A[i + j * rowA] = 4 * (i + 1);
  }
}

//==================  LagrangianScleronomousR ==================
extern "C" DLLEXPORT void hSclero(unsigned int sizeDS, const double* q, unsigned int sizeY, double* y, unsigned int sizeZ, double* z)
{
  printf("Call of the function 'hSclero' of the test plugin.\n");
}

extern "C" DLLEXPORT void G0Sclero(unsigned int sizeDS, const double* q, unsigned int sizeY, double* G0, unsigned int sizeZ, double* z)
{
  printf("Call of the function 'G0' of the test plugin.\n");
}

//==================  LagrangianRheonomousR ==================

extern "C" DLLEXPORT void hRheo(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*)
{
  printf("Call of the function 'hRheo' of the test plugin.\n");
}

extern "C" DLLEXPORT void G0Rheo(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*)
{
  printf("Call of the function 'G0Rheo' of the test plugin.\n");
}

extern "C" DLLEXPORT void hDot(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*)
{
  printf("Call of the function 'hDot' of the test plugin.\n");
}

//==================  LagrangianCompliantR ==================

extern "C" DLLEXPORT void hCompl(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Call of the function 'hCompl' of the test plugin.\n");
}

extern "C" DLLEXPORT void G0Compl(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Call of the function 'G0Compl' of the test plugin.\n");
}

extern "C" DLLEXPORT void G1Compl(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Call of the function 'G1Compl' of the test plugin.\n");
}

// ========== FirstOrderType1R ==========

extern "C" DLLEXPORT void hT1(unsigned int, const double*, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'hT1' of the test plugin.\n");
}

extern "C" DLLEXPORT void gT1(unsigned int, const double*, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'gT1' of the test plugin.\n");
}

extern "C" DLLEXPORT void Jh0T1(unsigned int, const double*, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'Jh0T1' of the test plugin.\n");
}

extern "C" DLLEXPORT void Jg0T1(unsigned int, const double*, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'Jg0T1' of the test plugin.\n");
}

//==================  FirstOrderLinearR ==================

extern "C" DLLEXPORT void C(double, unsigned int, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'C' of the test plugin.\n");
}

extern "C" DLLEXPORT void D(double, unsigned int, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'D' of the test plugin.\n");
}

extern "C" DLLEXPORT void F(double, unsigned int, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'F' of the test plugin.\n");
}

extern "C" DLLEXPORT void e(double, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'e' of the test plugin.\n");
}

extern "C" DLLEXPORT void B(double, unsigned int, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'B' of the test plugin.\n");
}
