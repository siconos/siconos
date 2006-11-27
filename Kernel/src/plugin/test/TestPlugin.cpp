/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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

/** DynamicalSystem plug-in to compute f(x,t) - id="f".
 *  @param sizeOfX : the size of the vector x
 *  @param time : current time
 *  @param x : the pointer to the first element of the vector x
 *  @param[in,out] f : the pointer to the first element of the vector f(x,t)
 *  @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void computeF(unsigned int sizeOfX, double time, const double* x, double* f, double* param)
{
  for (unsigned int i = 0; i < sizeOfX; i++)
    f[i] = time * x[i];
}

/** DynamicalSystem plug-in to compute the gradient of f(x,t) with respect to the state: \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$ - id = "jacobianXF"
 * @param sizeOfX : size of vector x
 * @param time : current time
 * @param x : pointer to the first element of x
 * @param[in,out] jacob : pointer to the first element of jacobianXF matrix
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void computeJacobianXF(unsigned int sizeOfX, double time, const double *x, double *jacob, double* param)
{
  printf("Call of the function 'computeJacobianX' of the basic plugin.\nYou have to implement this function.\n");
}

/** DynamicalSystem plug-in to compute u(x,t) - id = "u"
 * @param sizeOfU : size of vector u
 * @param sizeOfX : size of vector x
 * @param time : current time
 * @param x : pointer to the first element of x
 * @param[in,out] u : pointer to the first element of u vector
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void computeU(unsigned int sizeOfU, unsigned int sizeOfX, double time, const double* x, double* u, double* param)
{
  for (unsigned int i = 0; i < sizeOfU; i++)
  {
    if (i < sizeOfX)
      u[i] = time * x[i];
    else
      u[i] = 0;
  }
}

/** DynamicalSystem plug-in to compute T(x) - id = "T"
 * @param sizeOfU : size of vector u
 * @param sizeOfX : size of vector X
 * @param x : pointer to the first element of X
 * @param[in,out] T : pointer to the first element of T matrix
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void computeT(unsigned int sizeOfU, unsigned int  sizeOfX, const double* x, double* T, double* param)
{
  for (unsigned int j = 0; j < sizeOfU; j++)
  {
    for (unsigned int i = 0; i < sizeOfX; i++)
      T[i + j * sizeOfX] =  x[i];
  }
}

// ===== Lagrangian DS  =====

// Plugins for Fext, Fint, NNL (vectors), Mass, JacobianQNNL, JacobianVelocityNNL,
// JacobianQFint and JacobianVelocityFint (matrices)

/** LagrangianDS plug-in to compute mass(q,t) - id = "mass"
 * @param sizeOfq : size of vector q
 * @param q : pointer to the first element of q
 * @param[in,out] mass : pointer to the first element of mass
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void computeMass(unsigned int sizeOfq, const double *q, double *mass, double* param)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    mass[i] = 0;
  mass[0] = 1;
  mass[4] = 2;
  mass[8] = 3;
}

/** LagrangianDS plug-in to compute internal forces \f$F_{int}(t,q,\dot q)\f$ - id = "fInt"
 * @param sizeOfq : size of vector q
 * @param time : current time
 * @param q : pointer to the first element of q
 * @param velocity : pointer to the first element of velocity
 * @param[in,out] fInt : pointer to the first element of fInt
 * @param[in,out] param  : a vector of user-defined parameters
 */
extern "C" void computeFInt(unsigned int sizeOfq, double time, const double *q, const double *velocity, double *fInt, double * param)
{
  for (unsigned int i = 0; i < sizeOfq; ++i)
    fInt[i] = i * q[i];
}

/** LagrangianDS plug-in to compute external forces \f$F_{Ext}(t)\f$, id = "fExt"
 * @param sizeOfq : size of vector q
 * @param time : current time
 * @param[in,out] fExt : pointer to the first element of fExt
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void computeFExt(unsigned int sizeOfq, double time, double *fExt, double *param)
{
  for (unsigned int i = 0; i < sizeOfq; ++i)
    fExt[i] = i * time;
}

/** LagrangianDS plug-in to compute \f$NNL(\dot q, q)\f$, id = "NNL"
 * @param sizeOfq : size of vector q
 * @param q : pointer to the first element of q
 * @param velocity : pointer to the first element of velocity
 * @param[in,out] NNL : pointer to the first element of NNL
 * @param[in,out] param  : a vector of user-defined parameters
 */
extern "C" void computeNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *NNL, double *param)
{
  for (unsigned int i = 0; i < sizeOfq; ++i)
    NNL[i] = i * q[i];
}

/** LagrangianDS plug-in to compute \f$\nabla_qF_{Int}(\dot q, q, t)\f$, id = "jacobianQFInt"
 * @param sizeOfq : size of vector q
 * @param time : current time
 * @param q : pointer to the first element of q
 * @param velocity : pointer to the first element of velocity
 * @param[in,out] jacob : pointer to the first element of the jacobian
 * @param[in,out] param  : a vector of user-defined parameters
 */
extern "C" void computeJacobianQFInt(unsigned int sizeOfq, double time, const double *q, const double *velocity, double *jacob, double* param)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}

/** LagrangianDS plug-in to compute \f$\nabla_{\dot q}F_{Int}(\dot q, q, t)\f$, id = "jacobianVelocityFInt"
 * @param sizeOfq : size of vector q
 * @param time : current time
 * @param q : pointer to the first element of q
 * @param velocity : pointer to the first element of velocity
 * @param[in,out] jacob : pointer to the first element of the jacobian
 * @param[in,out] param  : a vector of user-defined parameters
 */
extern "C" void computeJacobianVelocityFInt(unsigned int sizeOfq, double time, const double *q, const double *velocity, double *jacob, double* param)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}

/** LagrangianDS plug-in to compute \f$\nabla_qNNL(\dot q, q)\f$, id = "jacobianQNNL"
 * @param sizeOfq : size of vector q
 * @param q : pointer to the first element of q
 * @param velocity : pointer to the first element of velocity
 * @param[in,out] jacob : pointer to the first element of the jacobian
 * @param[in,out] param  : a vector of user-defined parameters
 */
extern "C" void computeJacobianQNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, double* param)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}

/** LagrangianDS plug-in to compute \f$\nabla_{\dot q}NNL(\dot q, q)\f$, id = "jacobianVelocityNNL"
 * @param sizeOfq : size of vector q
 * @param q : pointer to the first element of q
 * @param velocity : pointer to the first element of velocity
 * @param[in,out] jacob : pointer to the first element of the jacobian
 * @param[in,out] param  : a vector of user-defined parameters
 */
extern "C" void computeJacobianVelocityNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, double* param)
{
  for (unsigned int i = 0; i < (sizeOfq * sizeOfq); ++i)
    jacob[i] = i * q[0];
}


// ===== Linear DS  ====

// Plugins for A, B (matrices), u and f. See LinearDS.h

/** LinearDS plug-in to compute b(t), id = "b"
 * @param sizeOfB : size of vector b
 * @param time : current time
 * @param[in,out] b : pointer to the first element of b
 * @param[in,out] param  : a vector of user-defined parameters
 */
extern "C" void computeB(unsigned int sizeOfB, double time, double* b, double *param)
{
  for (unsigned int i = 0; i < sizeOfB; i++)
    b[i] = time * i ;

}

/** LinearDS plug-in to compute A(t), id = "A"
 * @param sizeOfA : size of square-matrix A
 * @param time : current time
 * @param[in,out] A : pointer to the first element of A
 * @param[in,out] param  : a vector of user-defined parameters
 */
extern "C" void computeA(unsigned int  sizeOfA, double time, double* A, double *param)
{
  for (unsigned int j = 0; j < sizeOfA; j++)
  {
    for (unsigned int i = 0; i < sizeOfA; i++)
      A[i + j * sizeOfA] = 4 * (i + 1);
  }
}

// ===== RELATIONS ====

/** Relation plug-in to compute y(x,t) - id="output".
 *  @param sizeOfX : the size of the vector x
 *  @param x : the pointer to the first element of the vector x
 *  @param time : current time
 *  @param sizeOfY : the size of the vector y and lambda (ie of the interaction)
 *  @param lambda : the pointer to the first element of the vector lambda
 *  @param sizeOfU : the size of the vector u
 *  @param u : the pointer to the first element of the vector u
 *  @param[in,out]  y : the pointer to the first element of the vector y
 *  @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void y(unsigned int sizeOfX, const double* x, double time, unsigned int sizeOfY, const double* lambda,
                  unsigned int sizeOfU, const double* u, double* y, double* param)
{
  printf("Warning: call of the function 'computeOutput' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

/** Relation plug-in to compute r(lambda,t) - id="input".
 *  @param sizeY : the size of the vector y and lambda (ie of the interaction)
 *  @param lambda : the pointer to the first element of the vector lambda
 *  @param time : current time
 *  @param[in,out] r : the pointer to the first element of the vector y
 *  @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void R(unsigned int sizeY, const double* lambda, double time, double* r, double* param)
{
  printf("Warning: call of the function 'computeInput' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

// === Lagrangian Relations ===

/** LagrangianR plug-in to compute h0(q) (scleronomic case) - id="h"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] y : pointer to the first element of y
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void h0(unsigned int sizeDS, const double* q, unsigned int sizeY, double* y, double* param)
{
  printf("Call of the function 'h0' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianR plug-in to compute G0(q), gradient of h0 according to q (scleronomic case) - id="G0"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] G0 : pointer to the first element of G0 (sizeY X sizeDS matrix)
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void G0(unsigned int sizeDS, const double* q, unsigned int sizeY, double* G0, double* param)
{
  printf("Call of the function 'G0' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianR plug-in to compute h1(q,t) (non-scleronomic case) - id="h"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param time : current time
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] y : pointer to the first element of y
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C"  void h1(unsigned int sizeDS, const double* q, double time, unsigned int sizeY, double* y, double* param)
{
  printf("Call of the function 'h1' of the default plugin.\nYou have to implement this function.\n");
}


/** LagrangianR plug-in to compute G10(q,t), gradient of h1 accoring to q (non-scleronomic case) - id="G10"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param time : current time
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] G10 : pointer to the first element of G10
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void G10(unsigned int sizeDS, const double* q, double time, unsigned int  sizeY, double* G10, double* param)
{
  printf("Call of the function 'G10' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianR plug-in to compute G11(q,t), gradient of h1 according to time (non-scleronomic case) - id="G11"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param time : current time
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] G11 : pointer to the first element of G11
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void G11(unsigned int sizeDS, const double* q, double time, unsigned int  sizeY, double* G11, double* param)
{
  printf("Call of the function 'G11' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianR plug-in to compute h2(q,lambda) (scleronomic+lambda case) - id="h"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param lambda : current time
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] y : pointer to the first element of y
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void h2(unsigned int  sizeDS, const double* q, const double* lambda, unsigned int  sizeY, double* y, double* param)
{
  printf("Call of the function 'h2' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianR plug-in to compute G20(q,lambda), gradient of h2 according to q (scleronomic+lambda case) - id="G20"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param lambda : current time
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] G20 : pointer to the first element of G20
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void G20(unsigned int  sizeDS, const double* q, const double* lambda, unsigned int  sizeY, double* y, double* param)
{
  printf("Call of the function 'G20' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianR plug-in to compute G21(q,lambda), gradient of h2 according to lambda (scleronomic+lambda case) - id="G21"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param lambda : current time
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] G21 : pointer to the first element of G21
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void G21(unsigned int  sizeDS, const double* q, const double* lambda, unsigned int  sizeY, double* y, double* param)
{
  printf("Call of the function 'G21' of the default plugin.\nYou have to implement this function.\n");
}
