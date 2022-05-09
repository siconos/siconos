/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#ifndef GlobalRollingFrictionContact_H
#define GlobalRollingFrictionContact_H

#include "GlobalFrictionContact.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "GlobalRollingFrictionContactProblem.h"
#include "Friction_cst.h"

/** Pointer to function of the type used for drivers for GlobalRollingFrictionContact problems in Numerics */
typedef int (*GRFC3D_Driver)(GlobalRollingFrictionContactProblem*, double*, double*, double*, SolverOptions*);
TYPEDEF_SPTR(GlobalRollingFrictionContactProblem)

/** Formalization and Resolution of a Friction-Contact Problem
 *
 * This class is devoted to the formalization and the resolution of
 * primal friction contact problems defined by :
 *

 \rst

 .. math::

   M velocity =  q + H reaction \\
   globalVelocities = H^T velocity + tildeGlobalVelocities

 \endrst

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
 *  - solving of the GlobalRollingFrictionContact problem: function compute(), used to call solvers from Numerics through \n
 * the frictionContact2D_driver() or frictionContact3D_driver() interface of Numerics.
 *  - post-treatment of data: set values of y/lambda variables of the active Interaction (ie Interactions) using \n
 *  ouput results from the solver (velocity,reaction); function postCompute().
 *
 */
class GlobalRollingFrictionContact : public GlobalFrictionContact
{
private:
  /** default constructor */
  GlobalRollingFrictionContact() = default;

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(GlobalRollingFrictionContact);

  /** Type (dimension) of the contact problem (2D or 3D) */
  int _contactProblemDim = 5;

  /** rolling friction coefficients */
  SP::MuStorage _mu_r;

  /** Pointer to the function used to call the Numerics driver to solve the problem */
  GRFC3D_Driver _g_rolling_driver;

  GlobalRollingFrictionContactProblem _numerics_problem;
public:

  /** constructor (solver id and dimension)
      \param dimPb dimension (2D or 3D) of the friction-contact problem
      \param numericsSolverId id of the solver to be used, optional,
      default : SICONOS_GLOBAL_FRICTION_3D_NSGS
      \rst
      see :ref:`problems_and_solvers` for details.
      \endrst
  */
  GlobalRollingFrictionContact(int dimPb, int numericsSolverId = SICONOS_GLOBAL_ROLLING_FRICTION_3D_NSGS_WR);

  /**  constructor from a pre-defined solver options set.
       \param options, the options set,
       \rst
       see :ref:`problems_and_solvers` for details.
       \endrst
  */
  GlobalRollingFrictionContact(int dimPb, SP::SolverOptions options);

  /** destructor
   */
  virtual ~GlobalRollingFrictionContact(){};

  // GETTERS/SETTERS

  /** get the type of GlobalRollingFrictionContact problem (2D or 3D)
   *  \return an int (2 or 3)
   */
  inline int getGlobalRollingFrictionContactDim() const
  {
    return _contactProblemDim;
  }

  /** get dimension of the problem
   *  \return an unsigned ing
   */
  inline size_t getGlobalSizeOutput() const
  {
    return _sizeGlobalOutput;
  }

  /** get a pointer to mu, the list of the friction coefficients
   *  \return pointer on a std::vector<double>
   */
  inline SP::MuStorage mur() const
  {
    return _mu_r;
  }

  /** get the value of the component number i of mu, the vector of the friction coefficients
   *  \return the friction coefficient for the ith contact
   */
  inline double getMur(unsigned int i) const
  {
    return (*_mu_r)[i];
  }

  // --- Others functions ---

  /** initialize the GlobalRollingFrictionContact problem(compute topology ...)
   * \param sim the simulation, owner of this OSNSPB
   */
   void initialize(SP::Simulation sim);

  /**
   * \return the friction contact problem from Numerics
   */
  SP::GlobalRollingFrictionContactProblem globalRollingFrictionContactProblem();

  /**
   * \return the friction contact problem from Numerics (raw ptr, do not free)
   */
  GlobalRollingFrictionContactProblem *globalRollingFrictionContactProblemPtr();

  /** solve a friction contact problem
   * \param problem the friction contact problem
   * \return info solver information result
   */
  int solve(SP::GlobalRollingFrictionContactProblem problem = SP::GlobalRollingFrictionContactProblem());

  /** Construction of the problem
   *  \param time current time
   */
  virtual bool preCompute(double time);

  /** Compute the unknown reaction and velocity and update the Interaction (y and lambda )
   *  \param time current time
   */
  virtual int compute(double time);

  /* Check the compatibility fol the nslaw with the targeted OSNSP */
  bool checkCompatibleNSLaw(NonSmoothLaw& nslaw);

  void updateMu();
  void updateMur();
  void display() const;
};

#endif // GlobalRollingFrictionContact_H
