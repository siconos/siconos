// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file
  Fricton-Contact Non-Smooth Problem
*/
#ifndef GENERICMECHANICAL_H
#define GENERICMECHANICAL_H

#include "FrictionContact_options.h"  // contains only enum. Ok. For SICONOS_FRICTION_3D_ONECONTACT_NSN
#include "GenericMechanical_cst.h"  // for SICONOS_GENERIC_MECHANICAL_NSGS
#include "LinearOSNS.hpp"

struct GenericMechanicalProblem;
struct SolverOptions;

namespace siconos::nonsmooth_formulations {
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

    For details regarding the available options, see Nonsmooth problems formulations and
   available solvers in users' guide.
*/
class GenericMechanical : public LinearOSNS {
 protected:
  ACCEPT_SERIALIZATION(GenericMechanical);

  GenericMechanicalProblem* _pnumerics_GMP{nullptr};

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
  GenericMechanical(std::shared_ptr<SolverOptions> options);

  /** destructor
   */
  ~GenericMechanical() noexcept;

  // GETTERS/SETTERS

  // --- Others functions ---

  /** initialize the GenericMechanical problem(compute topology ...)
   *
   *  \param sim the simulation, owner of this OSNSPB
   */
  void initialize(std::shared_ptr<siconos::simulation::Simulation> sim) override;

  /** Compute the unknown reaction and velocity and update the Interaction (y
   *  and lambda )
   *
   *  \param time double current time
   *  \return int information about the solver convergence (0: ok, >0 problem, see Numerics
   * documentation)
   */
  int compute(double time) override;

  /** compute extra-diagonal interactionBlock-matrix
   *
   *  \param ed an edge descriptor
   */
  void computeInteractionBlock(
      const siconos::graphs::InteractionsGraph::EDescriptor& ed) override;

  /** compute diagonal Interaction block
   *
   *  \param vd  a vertex descriptor
   */
  void computeDiagonalInteractionBlock(
      const siconos::graphs::InteractionsGraph::VDescriptor& vd) override;

  /** print the data to the screen */
  void display() const override;

  /**
      compute interactionBlocks if necessary (this depends on the type of
      OSNS, on the indexSets ...)
   */
  void updateInteractionBlocks(siconos::graphs::InteractionsGraph& indexSet) override;

  /** Check the compatibility fol the nslaw with the targeted OSNSP */
  bool checkCompatibleNSLaw(siconos::modeling::NonSmoothLaw& nslaw) override;
};
}  // namespace siconos::nonsmooth_formulations

#endif  // GenericMechanical_H
