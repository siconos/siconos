/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
extern "C" void computeF(unsigned int sizeOfX, double time, const double* x, double* f, double* z)
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
extern "C" void computeJacobianXF(unsigned int sizeOfX, double time, const double *x, double *jacob, double* z)
{
  for (unsigned int i = 0; i < sizeOfX * sizeOfX; i++)
    jacob[i] = i + 1;
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
extern "C" void computeMass(unsigned int sizeOfq, const double *q, double *mass, double* z)
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
extern "C" void computeJacobianQFInt(unsigned int sizeOfq, double time, const double *q, const double *velocity, double *jacob, double* z)
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
extern "C" void computeJacobianVelocityFInt(unsigned int sizeOfq, double time, const double *q, const double *velocity, double *jacob, double* z)
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
extern "C" void computeJacobianQNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, double* z)
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
extern "C" void computeJacobianVelocityNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, double* z)
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

// === Lagrangian Relations ===

// Scleronomous

/** LagrangianR plug-in to compute h(q,z)
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] y : pointer to the first element of y
 * @param sizeZ : size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void hSclero(unsigned int sizeDS, const double* q, unsigned int sizeY, double* y, unsigned int sizeZ, double* z)
{
  printf("Call of the function 'hSclero' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianR plug-in to compute G0(q,z), gradient of h0 according to q
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] G0 : pointer to the first element of G0 (sizeY X sizeDS matrix)
 * @param sizeZ : size of vector z
 * @param[in,out] param : a vector of user-defined parameters
 */
extern "C" void G0Sclero(unsigned int sizeDS, const double* q, unsigned int sizeY, double* G0, unsigned int sizeZ, double* z)
{
  printf("Call of the function 'G0' of the default plugin.\nYou have to implement this function.\n");
}

// Rheonomous

/** LagrangianRheonomousR plug-in to compute h(q,t,z).
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param time : current time
 * @param sizeY : size of vector y (ie of the intercation)
 * @param[in,out] y : pointer to the first element of y
 * @param sizeZ : size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void hRheo(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*)
{
  printf("Call of the function 'hRheo' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianRheonomousR plug-in to compute G0(q,t,z) gradient of h according to q
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param time : current time
 * @param sizeY : size of vector y (ie of the intercation)
 * @param[in,out] pointer to the first element of G0
 * @param sizeZ : size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void G0Rheo(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*)
{
  printf("Call of the function 'G0Rheo' of the default plugin.\nYou have to implement this function.\n");
}


/** LagrangianRheonomousR plug-in to compute hDot(q,t,z)
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param time : current time
 * @param sizeY : size of vector y (ie of the intercation)
 * @param[in,out] pointer to the first element of hDot
 * @param sizeZ : size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void hDot(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*)
{
  printf("Call of the function 'hDot' of the default plugin.\nYou have to implement this function.\n");
}

// Compliant

/** LagrangianR plug-in to compute h(q,lambda,z)
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param sizeY : size of vector y (ie of lambda and of the interaction)
 * @param lambda : pointer to lambda of the interaction
 * @param[in,out] y : pointer to the first element of y
 * @param sizeZ : size of vector z.
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void hCompl(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Call of the function 'hCompl' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianR plug-in to compute G0(q,lambda,z), gradient of hCompl according to q.
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param sizeY : size of vector y (ie of lambda and of the interaction)
 * @param lambda : pointer to lambda of the interaction
 * @param[in,out] G0 : pointer to the first element of G0
 * @param sizeZ : size of vector z.
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void G0Compl(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Call of the function 'G0Compl' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianR plug-in to compute G1(q,lambda), gradient of hCompl according to lambda.
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param sizeY : size of vector y (ie of lambda and of the interaction)
 * @param lambda : pointer to lambda of the interaction
 * @param[in,out] G1 : pointer to the first element of G1
 * @param sizeZ : size of vector z.
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void G1Compl(unsigned int, const double*, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Call of the function 'G1Compl' of the default plugin.\nYou have to implement this function.\n");
}

// ========== FirstOrderRelations ==========

extern "C" void y(unsigned int, const double*, double, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'y' of the test plugin.\n");
}

extern "C" void R(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'R' of the test plugin.\n");
}

extern "C" void Jh0(unsigned int, const double*, double, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'Jh0' of the test plugin.\n");
}
extern "C" void Jh1(unsigned int, const double*, double, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'Jh1' of the test plugin.\n");
}
extern "C" void Jg0(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'Jg0' of the test plugin.\n");
}

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
