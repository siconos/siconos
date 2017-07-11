/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
/*! \file
  Primal Fricton-Contact Non-Smooth Problem
*/
#ifndef GlobalFrictionContact_H
#define GlobalFrictionContact_H

#include "LinearOSNS.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "GlobalFrictionContactProblem.h"
#include "Friction_cst.h"

/** Pointer to function of the type used for drivers for GlobalFrictionContact problems in Numerics */
typedef int (*GFC3D_Driver)(GlobalFrictionContactProblem*, double*, double*, double*, SolverOptions*);
TYPEDEF_SPTR(GlobalFrictionContactProblem)

/** Formalization and Resolution of a Friction-Contact Problem
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Dec 15, 2005
 *
 * This class is devoted to the formalization and the resolution of
 * primal friction contact problems defined by :
 * \f{eqnarray*}
 *  M velocity =  q + H reaction \\
 *  globalVelocities = H^T velocity + tildeGlobalVelocities\\
 * \f}
 * and \f$globalVelocities, reaction\f$ belongs to the Coulomb friction law with unilateral contact.
 *
 * With:
 *    - \f$velocity \in R^{n} \f$  and \f$reaction \in R^{n} \f$ the unknowns,
 *    - \f$M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
 *    - \f$globalVelocities \in R^{m} \f$  and \f$reaction \in R^{m} \f$ the unknowns,
 *    - \f$tildeGlobalVelocities \in R^{m} \f$ is the modified local velocity (\f$ e U_{N,k}\f$)
 *    - \f$M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
 *    - \f$H \in R^{n \times m } \f$
 *
 * The dimension of the problem (2D or 3D) is given by the variable contactProblemDim and the right
 * Numerics driver will be called according to this value.
 *
 * \b Construction:
 *   - Constructor from data (inputs = Simulations*, id, SP::NonSmoothSolver) - The solver is optional.
 * Main functions:
 *
 * \b Main functions:
 *  - formalization of the problem: computes M,q using the set of "active" Interactions from the simulation and \n
 *  the interactionBlock-matrices saved in the field interactionBlocks.\n
 *  Functions: initialize(), computeInteractionBlock(), preCompute()
 *  - solving of the GlobalFrictionContact problem: function compute(), used to call solvers from Numerics through \n
 * the frictionContact2D_driver() or frictionContact3D_driver() interface of Numerics.
 *  - post-treatment of data: set values of y/lambda variables of the active Interaction (ie Interactions) using \n
 *  ouput results from the solver (velocity,reaction); function postCompute().
 *
 */
class GlobalFrictionContact : public LinearOSNS
{
private:
  /** default constructor */
  GlobalFrictionContact() {};
  
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(GlobalFrictionContact);


  /** Type (dimension) of the contact problem (2D or 3D) */
  int _contactProblemDim;

  /** size of the local problem to solve */
  size_t _sizeGlobalOutput;

  /** contains the vector globalVelocities of a GlobalFrictionContact system */
  SP::SiconosVector _globalVelocities;

  /** contains the impact contributions */
  SP::SiconosVector _b;

  /** contains the matrix H of a GlobalFrictionContact system */
  SP::OSNSMatrix _H;

  /** friction coefficients */
  SP::MuStorage _mu;

  /** Pointer to the function used to call the Numerics driver to solve the problem */
  GFC3D_Driver _gfc_driver;


public:

  /** constructor from data
   *  \param dimPb dimension (2D or 3D) of the friction-contact problem
   *  \param numericsSolverId solver to be used (see the documentation of siconos/numerics)
   */
  GlobalFrictionContact(int dimPb, int numericsSolverId = SICONOS_GLOBAL_FRICTION_3D_NSGS);

  /** destructor
   */
  virtual ~GlobalFrictionContact();

  // GETTERS/SETTERS

  /** get the type of GlobalFrictionContact problem (2D or 3D)
   *  \return an int (2 or 3)
   */
  inline int getGlobalFrictionContactDim() const
  {
    return _contactProblemDim;
  }

  /** get dimension of the problem
   *  \return an unsigned ing
   */
  inline unsigned int getGlobalSizeOutput() const
  {
    return _sizeGlobalOutput;
  }

  /** get globalVelocities
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector globalVelocities() const
  {
    return _globalVelocities;
  }

  /** set globalVelocities to pointer newPtr
   *  \param newPtr the new vector
   */
  inline void setGlobalVelocities(SP::SiconosVector newPtr)
  {
    _globalVelocities = newPtr;
  }

  // --- H ---

  /** get H
   *  \return pointer on a OSNSMatrix
   */
  inline SP::OSNSMatrix H() const
  {
    return _H;
  }

  /** set the value of H
   *  \param H the new matrix
   */
  void setH(SP::OSNSMatrix H) { _H = H;}

  /** get a pointer to mu, the list of the friction coefficients
   *  \return pointer on a std::vector<double>
   */
  inline SP::MuStorage mu() const
  {
    return _mu;
  }

  /** get the value of the component number i of mu, the vector of the friction coefficients
   *  \return the friction coefficient for the ith contact
   */
  inline double getMu(unsigned int i) const
  {
    return (*_mu)[i];
  }

  // --- Others functions ---
  /** initialize the _M and _H matrix */
  virtual void initOSNSMatrix();

  /** Memory allocation or resizing for z,w,q,b, globalVelocities */
  void initVectorsMemory();

  /** initialize the GlobalFrictionContact problem(compute topology ...)
    \param the simulation, owner of this OSNSPB
    */
  virtual void initialize(SP::Simulation sim);

  /** Construction of the problem
   *  \param time current time
   */
  virtual bool preCompute(double time);

  /** Compute the unknown reaction and velocity and update the Interaction (y and lambda )
   *  \param double current time
   *  \return int information about the solver convergence (0: ok, >0 problem, see Numerics documentation)
   */
  virtual int compute(double time);

  /** post-treatment of output from Numerics solver: \n
   *  set values of the unknowns of Interactions using (velocity,reaction)
   */
  virtual void postCompute();

  /** print the data to the screen */
  void display() const;
};

#endif // GlobalFrictionContact_H
