/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
/*! \file DefaultPlugin.cpp
  A list of all the available plug-in and their default behavior.

  Each of them corresponds to a pointer to private member function.

  To use a plug-in:

  - define your own function in myPlugin.cpp, with the correct signature (corresponding to the function as described in the list below)
  - call setXXXX("myPlugin.so", functionName)
  - then each time you call computeXXX(..), your function will be used.

  Example: set internal forces for a LagrangianDS
  in MyPlugin.cpp
  \code
  extern "C" void myFint(double time, unsigned int sizeQ, const double *q,
  const double *velocity, double *fInt, unsigned int sizeZ, double * z)
  {
  // work on fInt ...
  }
  \endcode

  And in main.cpp
  \code
  DynamicalSystem * myDS = new LagrangianDS(...);

  myDS->setComputeFIntFunction("myPlugin.so", "myFint");

  myDS->computeFInt(time); // call myFint(...) with current values of q ...
  \endcode

*/
#include <stdio.h>

//================== Dynamical Systems ==================

/** DynamicalSystem plug-in to compute \f$ g(t,\dot x,x,z) \f$
 *  @param  : current time
 *  @param  : the size of the vector x
 *  @param  : the pointer to the first element of the vector x[0]=\f$ x \f$
 *  @param  : the pointer to the first element of the vector x[1]=\f$ \dot x \f$
 *  @param  : the pointer to the first element of the vector g(t, ...)
 *  @param  : the size of the vector z
 *  @param  : a vector of parameters, z
 */
extern "C" void computeF(double time, unsigned int sizeOfX, const double* x, double* f, unsigned int sizeOfZ, double* z)
{
  printf("Warning: call of the function 'f' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

/** DynamicalSystem plug-in to compute jacobianG; computeJacobianGPtr[i] for jacobianG[i].
 *  @param  : current time
 *  @param  : the size of the vector x
 *  @param  : the pointer to the first element of the vector x[0]=\f$ x \f$
 *  @param  : the pointer to the first element of the vector x[1]=\f$ \dot x \f$
 *  @param  : the pointer to the first element of the vector g(t, ...)
 *  @param  : the size of the vector z
 *  @param  : a vector of parameters, z
 */
extern "C" void computeJacobianXF(double time, unsigned int sizeOfX, const double *x, double *jacob, unsigned int sizeOfZ, double* z)
{
  printf("Warning: call of the function 'jacobianXF' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

//================== LagrangianDS ==================

/** LagrangianDS plug-in to compute mass(q,t) - id = "mass"
 * @param sizeOfq : size of vector q
 * @param q : pointer to the first element of q
 * @param[in,out] mass : pointer to the first element of mass
 * @param  size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void computeMass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeOfZ, double* z)
{
  printf("Warning: call of the function 'computeMass' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

/** LagrangianDS plug-in to compute internal forces \f$F_{int}(t,q,\dot q)\f$ - id = "fInt"
 * @param time : current time
 * @param sizeOfq : size of vector q
 * @param q : pointer to the first element of q
 * @param velocity : pointer to the first element of velocity
 * @param[in,out] fInt : pointer to the first element of fInt
 * @param  size of vector z
 * @param[in,out] param  : a vector of user-defined parameters
 */
extern "C" void computeFInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *fInt, unsigned int sizeOfZ, double * z)
{
  printf("Warning: call of the function 'computeFInt' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

/** LagrangianDS plug-in to compute external forces \f$F_{Ext}(t)\f$, id = "fExt"
 * @param time : current time
 * @param sizeOfq : size of vector q
 * @param[in,out] fExt : pointer to the first element of fExt
 * @param  size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void computeFExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeOfZ, double *z)
{
  printf("Warning: call of the function 'computeFExt' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

/** LagrangianDS plug-in to compute \f$NNL(\dot q, q)\f$, id = "NNL"
 * @param sizeOfq : size of vector q
 * @param q : pointer to the first element of q
 * @param velocity : pointer to the first element of velocity
 * @param[in,out] NNL : pointer to the first element of NNL
 * @param  size of vector z
 * @param[in,out] z  : a vector of user-defined parameters
 */
extern "C" void computeNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *NNL, unsigned int sizeOfZ, double *z)
{
  printf("Warning: call of the function 'computeNNL' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

/** LagrangianDS plug-in to compute \f$\nabla_qF_{Int}(\dot q, q, t)\f$, id = "jacobianQFInt"
 * @param time : current time
 * @param sizeOfq : size of vector q
 * @param q : pointer to the first element of q
 * @param velocity : pointer to the first element of velocity
 * @param[in,out] jacob : pointer to the first element of the jacobian
 * @param  size of vector z
 * @param[in,out] z  : a vector of user-defined parameters
 */
extern "C" void computeJacobianQFInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  printf("Warning: call of the function 'computeJacobianQFInt' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

/** LagrangianDS plug-in to compute \f$\nabla_{\dot q}F_{Int}(\dot q, q, t)\f$, id = "jacobianVelocityFInt"
 * @param time : current time
 * @param sizeOfq : size of vector q
 * @param q : pointer to the first element of q
 * @param velocity : pointer to the first element of velocity
 * @param[in,out] jacob : pointer to the first element of the jacobian
 * @param  size of vector z
 * @param[in,out] z  : a vector of user-defined parameters
 */
extern "C" void computeJacobianVelocityFInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  printf("Warning: call of the function 'computeJacobianVelocityFInt' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

/** LagrangianDS plug-in to compute \f$\nabla_qNNL(\dot q, q)\f$, id = "jacobianQNNL"
 * @param sizeOfq : size of vector q
 * @param q : pointer to the first element of q
 * @param velocity : pointer to the first element of velocity
 * @param[in,out] jacob : pointer to the first element of the jacobian
 * @param  size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void computeJacobianQNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  printf("Warning: call of the function 'computeJacobianQNNL' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

/** LagrangianDS plug-in to compute \f$\nabla_{\dot q}NNL(\dot q, q)\f$, id = "jacobianVelocityNNL"
 * @param sizeOfq : size of vector q
 * @param q : pointer to the first element of q
 * @param velocity : pointer to the first element of velocity
 * @param[in,out] jacob : pointer to the first element of the jacobian
 * @param  size of vector z
 * @param[in,out] z  : a vector of user-defined parameters
 */
extern "C" void computeJacobianVelocityNNL(unsigned int sizeOfq, const double *q, const double *velocity, double *jacob, unsigned int sizeOfZ, double* z)
{
  printf("Warning: call of the function 'computeJacobianVelocityNNL' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

//==================  FirstOrderLinearDS ==================

/** FirstOrderLinearDS plug-in to compute A(t), id = "A"
 * @param time : current time
 * @param sizeOfA : size of square-matrix A
 * @param[in,out] A : pointer to the first element of A
 * @param  size of vector z
 * @param[in,out] z  : a vector of user-defined parameters
 */
extern "C" void computeA(double time, unsigned int  sizeOfA, double* A, unsigned int sizeOfZ, double *z)
{
  printf("Warning: call of the function 'computeB' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

/** FirstOrderLinearDS plug-in to compute b(t), id = "b"
 * @param time : current time
 * @param sizeOfB : size of vector b
 * @param[in,out] b : pointer to the first element of b
 * @param  size of vector z
 * @param[in,out] z  : a vector of user-defined parameters
 */
extern "C" void computeB(double time, unsigned int sizeOfB, double* b, unsigned int sizeOfZ, double *z)
{
  printf("Warning: call of the function 'computeB' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

//==================  FirstOrderR ==================

/** FirstOrderR plug-in to compute y=h(x,t) - id="h".
 *  @param the size of the vector x.
 *  @param x : the pointer to the first element of the vector x.
 *  @param time : current time.
 *  @param the size of the vectors y and lambda.
 *  @param lambda : the pointer to the first element of the vector lambda.
 *  @param[in,out]  y : the pointer to the first element of the vector y.
 *  @param the size of the vectors z.
 *  @param[in,out] z : a vector of user-defined parameters.
 */
extern "C" void computeOutput(unsigned int sizeOfX, const double* x, double time, unsigned int sizeOfY, const double* lambda,
                              double* y, unsigned int sizeOfZ, double* z)
{
  printf("Warning: call of the function 'computeOutput' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

/** Plug-in to compute \f$ \nabla_x h(x,t,lambda,z)\f$
 *  @param the size of the vector x.
 *  @param x : the pointer to the first element of the vector x.
 *  @param time : current time.
 *  @param the size of the vectors y and lambda.
 *  @param lambda : the pointer to the first element of the vector lambda.
 *  @param[in,out]  jacob : the pointer to the first element of the jacobian.
 *  @param the size of the vectors z.
 *  @param[in,out] z : a vector of user-defined parameters.
 */
extern "C" void Jh0(unsigned int, const double*, double, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'Jh0' of the default plugin.\n");
}

/** Plug-in to compute \f$ \nabla_\lambda h(x,t,lambda,z)\f$
 *  @param the size of the vector x.
 *  @param x : the pointer to the first element of the vector x.
 *  @param time : current time.
 *  @param the size of the vectors y and lambda.
 *  @param lambda : the pointer to the first element of the vector lambda.
 *  @param[in,out]  jacob : the pointer to the first element of the jacobian
 *  @param the size of the vectors z.
 *  @param[in,out] z : a vector of user-defined parameters.
 */
extern "C" void Jh1(unsigned int, const double*, double, unsigned int, const double*, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'Jh1' of the default plugin.\n");
}

/** FirstOrderR plug-in to compute r = g(lambda,t) - id="g".
 *  @param sizeY : the size of the vector y and lambda.
 *  @param lambda : the pointer to the first element of the vector lambda.
 *  @param time : current time.
 *  @param[in,out] r : the pointer to the first element of the vector r.
 *  @param the size of the vectors z
 *  @param[in,out] z : a vector of user-defined parameters.
 */
extern "C" void computeInput(unsigned int sizeY, const double* lambda, double time, double* r, unsigned int sizeZ, double* z)
{
  printf("Warning: call of the function 'computeInput' of the default plugin, which is not implemented. Add it in yourPlugin.cpp.\n");
}

/** Plug-in to compute \f$ \nabla_X g(lambda,t,z)\f$.
 *  @param sizeY : the size of the vector y and lambda.
 *  @param lambda : the pointer to the first element of the vector lambda.
 *  @param time : current time.
 *  @param[in,out] jacob : the pointer to the first element of the jacobian.
 *  @param the size of the vectors z
 *  @param[in,out] z : a vector of user-defined parameters.
 */
extern "C" void Jg0(unsigned int, const double*, double, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'Jg0' of the default plugin.\n");
}

//==================  FirstOrderLinearR ==================

/** Plug-in list to compute C(t,z)
 * @param time: current time
 * @param rowOfC, number of rows of in-out matrix C
 * @param colOfC, number of columns of in-out matrix C
 * @param[in,out] C : pointer to the first element of C
 * @param[in] sizeOfZ: size of vector z
 * @param[in,out] z: a vector of user-defined parameters
 */
extern "C" void C(double, unsigned int, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'C' of the default plugin.\n");
}

/** Plug-in to compute D(t,z)
 * @param time: current time
 * @param rowOfD, number of rows of in-out square matrix D
 * @param[in,out] D : pointer to the first element of D
 * @param[in] sizeOfZ: size of vector z
 * @param[in,out] z: a vector of user-defined parameters
 */
extern "C" void D(double, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'D' of the default plugin.\n");
}

/** Plug-in to compute F(t,z)
 * @param time: current time
 * @param rowOfF, number of rows of in-out matrix F
 * @param[in,out] F : pointer to the first element of F
 * @param[in] sizeOfZ: size of vector z
 * @param[in,out] z: a vector of user-defined parameters
 */
extern "C" void F(double, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'F' of the default plugin.\n");
}

/** Plug-in to compute e(t,z)
 * @param time: current time
 * @param sizeOfE, size of in-out vector e
 * @param[in,out] e : pointer to the first element of e
 * @param[in] sizeOfZ: size of vector z
 * @param[in,out] z: a vector of user-defined parameters
 */
extern "C" void e(double, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'e' of the default plugin.\n");
}

/** Plug-in to compute B(t,z)
 * @param time: current time
 * @param rowOfB, number of rows of in-out matrix B
 * @param colOfB, number of columns of in-out matrix B
 * @param[in,out] B : pointer to the first element of B
 * @param[in] sizeOfZ: size of vector z
 * @param[in,out] z: a vector of user-defined parameters
 */
extern "C" void B(double, unsigned int, unsigned int, double*, unsigned int, double*)
{
  printf("Warning: call of the function 'B' of the default plugin.\n");
}

//==================  LagrangianScleronomousR ==================

/** LagrangianScleronomousR plug-in to compute h0(q)  - id="h"
 * @param sizeQ: size of q = sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] y : pointer to the first element of y
 * @param sizeZ : size of vector z
 * @param[in,out] z: pointer to z vector(s) from DS.
 */
extern "C" void h0(unsigned int sizeDS, const double* q, unsigned int sizeY, double* y, unsigned int sizeZ, double* z)
{
  printf("Call of the function 'h0' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianScleronomousR plug-in to compute G0(q), gradient of h0 according to q - id="G0"
 * @param sizeQ: size of q = sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param sizeY : size of vector y (ie of the intercation)
 * @param[in,out] G0 : pointer to the first element of G0 (sizeY X sizeDS matrix)
 * @param sizeZ : size of vector z
 * @param[in,out] z: pointer to z vector(s) from DS.
 */
extern "C" void G0(unsigned int sizeDS, const double* q, unsigned int sizeY, double* G0, unsigned int sizeZ, double* z)
{
  printf("Call of the function 'G0' of the default plugin.\nYou have to implement this function.\n");
}

//==================  LagrangianRheonomousR ==================

/** LagrangianRheonomousR plug-in to compute h1(q,t) - id="h"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param time : current time
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] y : pointer to the first element of y
 * @param  size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C"  void h1(unsigned int sizeDS, const double* q, double time, unsigned int sizeY, double* y, unsigned int sizeOfZ, double* z)
{
  printf("Call of the function 'h1' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianRheonomousR plug-in to compute G10(q,t), gradient of h1 accoring to q - id="G10"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param time : current time
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] G10 : pointer to the first element of G10
 * @param  size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void G10(unsigned int sizeDS, const double* q, double time, unsigned int  sizeY, double* G10, unsigned int sizeOfZ, double* z)
{
  printf("Call of the function 'G10' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianRheonomousR plug-in to compute G11(q,t), gradient of h1 according to time - id="G11"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param time : current time
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] G11 : pointer to the first element of G11
 * @param  size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void G11(unsigned int sizeDS, const double* q, double time, unsigned int  sizeY, double* G11, unsigned int sizeOfZ, double* z)
{
  printf("Call of the function 'G11' of the default plugin.\nYou have to implement this function.\n");
}

//==================  LagrangianCompliantR ==================

/** LagrangianCompliantR plug-in to compute h2(q,lambda) - id="h"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param lambda : current time
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] y : pointer to the first element of y
 * @param  size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void h2(unsigned int  sizeDS, const double* q, unsigned int sizeY, const double* lambda, double* y, unsigned int sizeOfZ, double* z)
{
  printf("Call of the function 'h2' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianCompliantR plug-in to compute G20(q,lambda), gradient of h2 according to q - id="G20"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param lambda : current time
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] G20 : pointer to the first element of G20
 * @param  size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void G20(unsigned int  sizeDS, const double* q, unsigned int sizeY, const double* lambda, double* y, unsigned int sizeOfZ, double* z)
{
  printf("Call of the function 'G20' of the default plugin.\nYou have to implement this function.\n");
}

/** LagrangianCompliantR plug-in to compute G21(q,lambda), gradient of h2 according to lambda - id="G21"
 * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
 * @param q : pointer to the first element of q
 * @param lambda : current time
 * @param sizeY : size of vector y (ie of the interaction)
 * @param[in,out] G21 : pointer to the first element of G21
 * @param  size of vector z
 * @param[in,out] z : a vector of user-defined parameters
 */
extern "C" void G21(unsigned int  sizeDS, const double* q, unsigned int  sizeY, const double* lambda, double* y, unsigned int sizeOfZ, double* z)
{
  printf("Call of the function 'G21' of the default plugin.\nYou have to implement this function.\n");
}
