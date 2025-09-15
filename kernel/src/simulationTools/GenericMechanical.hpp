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
  Fricton-Contact Non-Smooth Problem
*/
#ifndef GENERICMECHANICAL_H
#define GENERICMECHANICAL_H

#include "Friction_cst.h" // from numerics, for solver id
#include "LinearOSNS.hpp"
#include "SiconosNumerics.h"

TYPEDEF_SPTR(GenericMechanicalProblem)

/** 
    Formalization and Resolution of a generic mechanical problem: It mixes
    bilateral equality, complementarity, impact and friction problems.
    
    This class is devoted to contains of a set of Non-Smooth Problem.
    
    \b Main functions:
    - formalization of the problem: computes M,q using the set of "active"
    Interactions from the simulation and the interactionBlock-matrices saved
    in the field interactionBlocks. 
    Functions: initialize(), computeInteractionBlock(), preCompute()
    - solving of the GenericMechanical problem: function compute(), used to call
    solvers from Numerics through the gmp_driver() interface of Numerics.
    - post-treatment of data: set values of y/lambda variables of the active
    Interaction (ie Interactions) using ouput results from the solver
    (velocity,reaction); function postCompute().
    
    For details regarding the available options, see Nonsmooth problems formulations and available solvers in users' guide.
*/
class GenericMechanical : public LinearOSNS {
protected:

  ACCEPT_SERIALIZATION(GenericMechanical);

  GenericMechanicalProblem *_pnumerics_GMP;

public:
  /** constructor from solver id
   *  
   *  \param numericsSolverId id of the internal friction solver of the generic
   *  problem default = SICONOS_FRICTION_3D_ONECONTACT_NSN
   */
  GenericMechanical(int FC3D_Solver_Id = SICONOS_FRICTION_3D_ONECONTACT_NSN);

  /** constructor from a pre-defined solver options set
   *
   *  \param options the options set
   */
  GenericMechanical(SP::SolverOptions options);

  /** destructor
   */
  ~GenericMechanical();

  // GETTERS/SETTERS

  // --- Others functions ---

  /** initialize the GenericMechanical problem(compute topology ...)
   *
   *  \param sim the simulation, owner of this OSNSPB
   */
  void initialize(SP::Simulation sim) override;

  /** Compute the unknown reaction and velocity and update the Interaction (y
   *  and lambda )
   *
   *  \param time double current time
   *  \return int information about the solver convergence (0: ok, >0 problem, see Numerics documentation)
   */
  int compute(double time) override;

  /** compute extra-diagonal interactionBlock-matrix
   *
   *  \param ed an edge descriptor
   */
  void
  computeInteractionBlock(const InteractionsGraph::EDescriptor &ed) override;

  /** compute diagonal Interaction block
   *
   *  \param vd  a vertex descriptor
   */
  void computeDiagonalInteractionBlock(
      const InteractionsGraph::VDescriptor &vd) override;

  /** print the data to the screen */
  void display() const override;

  /** 
      compute interactionBlocks if necessary (this depends on the type of
      OSNS, on the indexSets ...)
   */
  void updateInteractionBlocks() override;

  /** Check the compatibility fol the nslaw with the targeted OSNSP */
  bool checkCompatibleNSLaw(NonSmoothLaw &nslaw) override;

  ACCEPT_STD_VISITORS();
};

#endif // GenericMechanical_H
