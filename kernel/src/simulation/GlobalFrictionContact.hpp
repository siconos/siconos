/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "Friction_cst.h"  // contains only enum. Ok.
#include "LinearOSNS.hpp"

struct GlobalFrictionContactProblem;
struct SolverOptions;

namespace siconos::nonsmooth_formulations {

/**
   Formalization and Resolution of a Friction-Contact Problem

   This class is devoted to the formalization and the resolution of
   primal friction contact problems defined by :

   \f[

   M velocity =  q + H reaction \\
   globalVelocities = H^T velocity + tildeGlobalVelocities
   \f]

   and \f$ globalVelocities, reaction \f$ belongs to the Coulomb friction law with unilateral
   contact.

   With:
   - \f$ M \in R^{n \times n } \f$  and \f$ q \in R^{n} \f$
   - \f$ globalVelocities \in R^{m} and \f$ reaction \in R^{m} \f$ the unknowns,
   - \f$ tildeGlobalVelocities \in R^{m} \f$ is the modified local velocity (\f$ e U_{N,k} \f$)
   - \f$ H \in R^{n \times m } \f$

   n is given by the DS graph (sum of all dynamical systems sizes)
   m is given by the Interaction graph (sum of all interactions sizes)

   The dimension of the problem (2D or 3D) is given by the variable contactProblemDim and the
   right Numerics driver will be called according to this value.

   \b Construction:
   - Constructor from data (inputs = Simulations*, id, shared_ptr<NonSmoothSolver>) - The
   solver is optional. Main functions:

   \b Main functions:
   - formalization of the problem: computes M,q using the set of "active" Interactions from the
   simulation and the interactionBlock-matrices saved in the field interactionBlocks.
   Functions: initialize(), computeInteractionBlock(), preCompute()
   - solving of the GlobalFrictionContact problem: function compute(), used to call solvers
   from Numerics through the frictionContact2D_driver() or frictionContact3D_driver() interface
   of Numerics.
   - post-treatment of data: set values of y/lambda variables of the active Interaction (ie
   Interactions) using ouput results from the solver (velocity,reaction); function
   postCompute().

   For details regarding the available options, see Nonsmooth problems formulations and
   available solvers in users' guide.


*/
class GlobalFrictionContact : public LinearOSNS {
 protected:
  // /** default constructor */
  // GlobalFrictionContact() = default;

  /** Pointer to function of the type used for drivers for GlobalFrictionContact problems in
   * Numerics */
  typedef int (*GFC3D_Driver)(GlobalFrictionContactProblem*, double*, double*, double*,
                              SolverOptions*);

 protected:
  ACCEPT_SERIALIZATION(GlobalFrictionContact);

  /** Type (dimension) of the contact problem (2D or 3D) */
  int _contactProblemDim = 3;

  /** size of the local problem to solve */
  siconos::algebra::Index _sizeGlobalOutput = 0;

  /** contains the vector globalVelocities of a GlobalFrictionContact system */
  std::shared_ptr<siconos::algebra::SiconosVector> _globalVelocities{nullptr};

  /** contains the impact contributions */
  std::shared_ptr<siconos::algebra::SiconosVector> _b{nullptr};

  /** friction coefficients */
  std::shared_ptr<std::vector<double>> _mu{nullptr};

  /** Pointer to the function used to call the Numerics driver to solve the problem */
  GFC3D_Driver _gfc_driver;

  // GlobalFrictionContactProblem _numerics_problem;
 public:
  /** constructor (solver id and dimension)
   *
   *  \param dimPb dimension (2D or 3D) of the friction-contact problem
   *  \param numericsSolverId id of the solver to be used, optional,
   *  default : SICONOS_GLOBAL_FRICTION_3D_NSGS
   */
  GlobalFrictionContact(int dimPb, int numericsSolverId = SICONOS_GLOBAL_FRICTION_3D_NSGS);

  /** constructor from a pre-defined solver options set
   *
   *  \param options the options set
   */
  GlobalFrictionContact(int dimPb, std::shared_ptr<SolverOptions> options);

  /** destructor
   */
  virtual ~GlobalFrictionContact() noexcept = default;

  /** get the type of GlobalFrictionContact problem (2D or 3D)
   *
   *  \return an int (2 or 3)
   */
  inline int getGlobalFrictionContactDim() const { return _contactProblemDim; }

  /** get dimension of the problem
   *
   *  \return an unsigned ing
   */
  inline size_t getGlobalSizeOutput() const { return _sizeGlobalOutput; }

  /** get globalVelocities
   *
   *  \return pointer on a siconos::algebra::SiconosVector
   */
  inline std::shared_ptr<siconos::algebra::SiconosVector> globalVelocities() const {
    return _globalVelocities;
  }

  /** set globalVelocities to pointer newPtr
   *
   *  \param newPtr the new vector
   */
  inline void setGlobalVelocities(std::shared_ptr<siconos::algebra::SiconosVector> newPtr) {
    _globalVelocities = newPtr;
  }

  /** get a pointer to mu, the list of the friction coefficients
   *
   *  \return pointer on a std::vector<double>
   */
  inline std::shared_ptr<std::vector<double>> mu() const { return _mu; }

  /** get the value of the component number i of mu, the vector of the friction coefficients
   *
   *  \return the friction coefficient for the ith contact
   */
  inline double getMu(unsigned int i) const { return (*_mu)[i]; }

  // --- Others functions ---

  /** Memory allocation or resizing for z,w,q,b, globalVelocities */
  void initVectorsMemory();

  /** initialize the GlobalFrictionContact problem(compute topology ...)
   *
   *  \param sim the simulation, owner of this OSNSPB
   */
  virtual void initialize(std::shared_ptr<siconos::simulation::Simulation> sim) override;

  /** \return the friction contact problem from Numerics
   */
  std::shared_ptr<GlobalFrictionContactProblem> globalFrictionContactProblem();

  // /** \return the friction contact problem from Numerics (raw ptr, do not f
  // ree)
  //  */
  // GlobalFrictionContactProblem* globalFrictionContactProblemPtr();

  /** solve a friction contact problem
   *
   *  \return info solver information result
   */
  int solve();
  //  int solve(std::shared_ptr<GlobalFrictionContactProblem> problem = nullptr);

  /** compute vector q */
  virtual void compute_q() override;

  /** Construction of the problem
   *
   *  \param time current time
   */
  virtual bool preCompute(double time) override;

  /** Compute the unknown reaction and velocity and update the Interaction (y and lambda )
   *
   *  \param time current time
   */
  virtual int compute(double time) override;

  /** update dynamical systems state with the current global velocity values */
  virtual void update_dynamicalsystems_state();

  /** post-treatment of output from Numerics solver:
   *  set values of the unknowns of Interactions using (velocity,reaction)
   */
  virtual void postCompute() override;

  /* Check the compatibility fol the nslaw with the targeted OSNSP */
  bool checkCompatibleNSLaw(siconos::modeling::NonSmoothLaw& nslaw) override;

  void updateMu();

  /** print the data to the screen */
  void display() const override;
};
}  // namespace siconos::nonsmooth_formulations

#endif  // GlobalFrictionContact_H
