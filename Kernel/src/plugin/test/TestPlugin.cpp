/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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

// function to compute u
extern "C" void computeU(const unsigned int& sizeOfU, const unsigned int& sizeOfX, const double* time, const double* xPtr, double* UPtr, double* param)
{
  /* input parameter :
   *  sizeOfU (size of the vector u)
   *  sizeOfX (size of the vector x)
   *  time
   *  pointer to x
   *  pointer to u (in-out parameter)
   */

  for (unsigned int i = 0; i < sizeOfU; i++)
  {
    if (i < sizeOfX)
      UPtr[i] = *time * xPtr[i];
    else
      UPtr[i] = 0;
  }
}

// function to compute T
extern "C" void computeT(const unsigned int& sizeOfU, const unsigned int& sizeOfX, const double* xPtr, double* TPtr, double* param)
{
  /* input parameter :
   *  sizeOfU (size of the vector u)
   *  sizeOfX (size of the vector X)
   *  pointer to x
   *  pointer to T (in-out parameter) - T is a (sizeOfX x sizeOfU) matrix
   */
  for (unsigned int j = 0; j < sizeOfU; j++)
  {
    for (unsigned int i = 0; i < sizeOfX; i++)
      TPtr[i + j * sizeOfX] =  xPtr[i];
  }
}

// vectorField
extern "C" void vectorField(const unsigned int& sizeOfX, const double *time, const double *x, double *xdot, double* param)
{
  /* input parameter : sizeOfX (size of the vector X); time ; x (pointer to X vector);
   * output parameter : xdot (pointer to Xdot vector)
   */
  for (unsigned int i = 0; i < sizeOfX; i++)
    xdot[i] = *time * x[i];
}



extern "C" void computeJacobianX(const unsigned int &sizeOfX, const double *time, const double *x, double *jacob, double* param)
{
  /* input parameter : sizeOfX (size of the vector X); time; x (pointer to x vector);
   * output parameter : jacob (pointer to JacobianX matrix)
   */

  printf("Call of the function 'computeJacobianX' of the basic plugin.\nYou have to implement this function.\n");
}


// ===== Lagrangian DS  =====

// Plugins for Fext, Fint, NNL (vectors), Mass, JacobianQNNL, JacobianVelocityNNL,
// JacobianQFint and JacobianVelocityFint (matrices)

extern "C" void computeFInt(const unsigned int&sizeOfq, const double *time, const double *q, const double *velocity, double *fInt, double * param)
{
  for (unsigned int i = 0; i < sizeOfq; ++i)
    fInt[i] = i * q[i];
}


extern "C" void computeFExt(const unsigned int&sizeOfq, const double *time, double *fExt, double *param)
{
  for (unsigned int i = 0; i < sizeOfq; ++i)
    fExt[i] = i**time;
}

extern "C" void computeNNL(const unsigned int&sizeOfq, const double *q, const double *velocity, double *NNL, double *param)
{
  for (unsigned int i = 0; i < sizeOfq; ++i)
    NNL[i] = i * q[i];
}


extern "C" void computeMass(const unsigned int&sizeOfq, const double *q, double *mass, double* param)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    mass[i] = 0;
  mass[0] = 1;
  mass[4] = 2;
  mass[8] = 3;
}


extern "C" void computeJacobianQFInt(const unsigned int&sizeOfq, const double *time, const double *q, const double *velocity, double *jacob, double* param)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}

extern "C" void computeJacobianVelocityFInt(const unsigned int&sizeOfq, const double *time, const double *q, const double *velocity, double *jacob, double* param)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}

extern "C" void computeJacobianQNNL(const unsigned int&sizeOfq, const double *q, const double *velocity, double *jacob, double* param)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}

extern "C" void computeJacobianVelocityNNL(const unsigned int&sizeOfq, const double *q, const  double *velocity, double *jacob, double* param)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}


// ===== Linear DS  ====

// Plugins for A, B (matrices), u and f. See LinearDS.h

extern "C" void computeB(const unsigned int &sizeOfB, const double *time, double* b, double *param)
{
  for (unsigned int i = 0; i < sizeOfB; i++)
    b[i] = *time * i ;

}
extern "C" void computeA(const unsigned int &sizeOfX, const double *time, const double *x, double *jacob, double* param)
{
  /* input parameter : sizeOfX (size of the vector X); time; x (pointer to x vector);
   * output parameter : jacob (pointer to JacobianX matrix)
   */
  for (unsigned int j = 0; j < sizeOfX; j++)
  {
    for (unsigned int i = 0; i < sizeOfX; i++)
      jacob[i + j * sizeOfX] = 4 * (i + 1);
  }
}

// ===== RELATIONS ====

extern "C" void y(const unsigned int* sizeOfX, const double* x, const double* time, const unsigned int* sizeOfY, const double* lambda,
                  const unsigned int* sizeOfU, const double* u, double* y, double* param)
{
  /* input parameter : sizeOfX (size of the vector X); x (pointer to x vector); time; lambda (pointer to lambda vector)
   * output parameter : y (pointer to vector y )
   */
  printf("Warning: call of the function 'computeOutput' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

extern "C" void R(const unsigned int* sizeY, const double* lambda, const double* time, double* r, double* param)
{
  /* input parameter : sizeOfX (size of the vector X); x (pointer to x vector); time; lambda (pointer to lambda vector)
   * output parameter : r (pointer to vector r )
   */
  printf("Warning: call of the function 'computeInput' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

// === Lagrangian Relations ===

extern "C" void h0(const unsigned int* sizeDS, const double* q, const unsigned int* sizeY, double* y, double* param)
{
  printf("Call of the function 'h0' of the default plugin.\nYou have to implement this function.\n");
}

extern "C" void G0(const unsigned int* sizeDS, const double* q, const unsigned int* sizeY, double* G0, double* param)
{
  printf("Call of the function 'G0' of the default plugin.\nYou have to implement this function.\n");
}

extern "C"  void h1(const unsigned int* sizeDS, const double* q, const double* time, const unsigned int* sizeY, double* y, double* param)
{
  printf("Call of the function 'h1' of the default plugin.\nYou have to implement this function.\n");
}


extern "C"    void G10(const unsigned int* sizeDS, const double* q, const double* time, const unsigned int* sizeY, double* G0, double* param)
{
  printf("Call of the function 'G10' of the default plugin.\nYou have to implement this function.\n");
}


extern "C"    void G11(const unsigned int* sizeDS, const double* q, const double* time, const unsigned int* sizeY, double* G1, double* param)
{
  printf("Call of the function 'G11' of the default plugin.\nYou have to implement this function.\n");
}

extern "C"  void h2(const unsigned int* sizeDS, const double* q, const double* lambda, const unsigned int* sizeY, double* y, double* param)
{
  printf("Call of the function 'h2' of the default plugin.\nYou have to implement this function.\n");
}


extern "C"    void G20(const unsigned int* sizeDS, const double* q, const double* lambda, const unsigned int* sizeY, double* y, double* param)
{
  printf("Call of the function 'G20' of the default plugin.\nYou have to implement this function.\n");
}


extern "C"    void G21(const unsigned int* sizeDS, const double* q, const double* lambda, const unsigned int* sizeY, double* y, double* param)
{
  printf("Call of the function 'G21' of the default plugin.\nYou have to implement this function.\n");
}
