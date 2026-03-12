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
/*! \file OneStepNSProblem.hpp
  \brief Interface to formalize and solve Non-Sooth problems
 */

#ifndef ONESTEPNSPROBLEM_H
#define ONESTEPNSPROBLEM_H

// #include <memory>
#include <set>

#include "SiconosSerialization.hpp"
#include "SimulationGraphs.hpp"
#include "TypeName.hpp"

struct SolverOptions;

namespace siconos::simulation {
class Simulation;
}

namespace siconos::nonsmooth_formulations {

/**
 Non Smooth Problem Formalization and Simulation

 This is an abstract class, that provides an interface to define a
 non smooth problem:
 - a formulation (ie the way the problem is written)
 - a solver (algorithm and solving formulation, that can be
 different from problem formulation)
 - routines to compute the problem solution.

 Two types of problem formulation are available :
 - Quadratic Problem
 - Linear Problem

 See derived classes (QP and LinearOSNS) for details.

 For Linear problems, the following formulations exists:
 - Linear Complementarity (LCP)
 - Mixed Linear Complementarity (MLCP)
 - Affine Variational Inequalities (AVI)
 - FrictionContact
 - Relay
 - Equality
 - GenericMechanical
 - MultipleImpact
 - GlobalFrictionContact

 The usual way to build and initialize a one-step nonsmooth problem is :

 - call constructor with the id of the required Numerics solver.
 (see Solver class or Numerics documentation for details on algorithm name and parameters).
 - initialize(simulation)
 Initialize process is usually done through model->initialize(simulation).
 See Examples for practical details.

 Options for Numerics and the driver for solvers

 When the Numerics driver is called a set of solver options (name, tolerance, max. number of
 iterations ...) is required --> SolverOptions.

 Default values are always set in solver options the OneStepNSProblem is built
 but if you need to set them yourself, please check Users'guide, Numerics solvers part.

*/
class OneStepNSProblem {
 protected:
  ACCEPT_SERIALIZATION(OneStepNSProblem);

  /** Numerics solver properties */
  std::shared_ptr<SolverOptions> _numerics_solver_options;

  /** size of the nonsmooth problem */
  siconos::algebra::Index _sizeOutput = 0;

  /** link to the simulation that owns the nonsmooth problem */
  std::shared_ptr<siconos::simulation::Simulation> _simulation;

  /** level of index sets that is considered by this osnsp */
  unsigned int _indexSetLevel = 0;

  /** level of input and output variables of the nonsmooth problems.
   *  We consider that the osnsp computes y[_inputOutputLevel] and lambda[_inputOutputLevel]
   */
  unsigned int _inputOutputLevel = 0;

  /** maximum value for sizeOutput. Set to the number of declared
      constraints by default (topology->getNumberOfConstraints());
      This value is used to allocate memory for M during initialize
      call. The best choice is to set maxSize to the estimated maximum
      dimension of the problem. It must not exceed ...
  */
  siconos::algebra::Index _maxSize = 0;

  /** Internal struc used to check that one and only one "type" of nslaw is used
  Same "type" means same class of nslaw AND same size.
  */
  struct NSLTypeKey {
    siconos::modeling::Type t;
    siconos::algebra::Index size;

    // for std::set
    bool operator<(const NSLTypeKey& other) const {
      return std::tie(t, size) < std::tie(other.t, other.size);
    }
  };

  /* set of nslaw types */
  std::set<NSLTypeKey> _nslawtype = {};

  /*During Newton it, this flag allows to update the numerics matrices only once if
   * necessary.*/
  bool _hasBeenUpdated = false;

  OneStepNSProblem() = default;

 private:
  OneStepNSProblem(const OneStepNSProblem&) = delete;
  OneStepNSProblem(OneStepNSProblem&&) = delete;
  OneStepNSProblem& operator=(const OneStepNSProblem&) = delete;
  OneStepNSProblem& operator=(OneStepNSProblem&&) = delete;

 public:
  /** constructor from a pre-defined solver options set.
   *
   *  \param options the options set
   */
  OneStepNSProblem(std::shared_ptr<SolverOptions> options)
      : _numerics_solver_options(options){};

  /** destructor
   */
  virtual ~OneStepNSProblem() noexcept = default;

  /** To get the SolverOptions structure
   *
   *  \return , the numerics structure used to save solver parameters
   */
  inline std::shared_ptr<SolverOptions> numericsSolverOptions() const {
    return _numerics_solver_options;
  };

  /** returns the dimension of the nonsmooth problem
   */
  inline siconos::algebra::Index getSizeOutput() const { return _sizeOutput; }

  /** get the simulation which owns this nonsmooth problem
   *
   *  \return a pointer on Simulation
   */
  inline std::shared_ptr<siconos::simulation::Simulation> simulation() const {
    return _simulation;
  }

  /** set the Simulation of the OneStepNSProblem
   *
   *  \param newS a pointer to Simulation
   */
  inline void setSimulationPtr(std::shared_ptr<siconos::simulation::Simulation> newS) {
    _simulation = newS;
  }

  /** \return the index set level */
  inline unsigned int indexSetLevel() const { return _indexSetLevel; }

  /** set the value of level min
   *
   *  \param newVal new value
   */
  inline void setIndexSetLevel(unsigned int newVal) { _indexSetLevel = newVal; }

  /** \return the Input/Output level */
  inline unsigned int inputOutputLevel() const { return _inputOutputLevel; }

  /** set the value of Input/Output level
   *
   *  \param newVal new value
   */
  inline void setInputOutputLevel(unsigned int newVal) { _inputOutputLevel = newVal; }

  /** \return the maximum value allowed for the dimension of the problem */
  inline siconos::algebra::Index maxSize() const { return _maxSize; }

  /** set the value of maxSize
   *
   *  \param newVal new value
   */
  inline void setMaxSize(const siconos::algebra::Index newVal) { _maxSize = newVal; }

  /** Turn on/off verbose mode in numerics solver*/
  void setNumericsVerboseMode(bool vMode);

  /**  set the verbose level in numerics solver*/
  void setNumericsVerboseLevel(int level);

  /**
      Check if the OSNSPb has interactions

      \return bool = true if the  osnsp has interactions, i.e. indexSet(_indexSetLevel)->size
     >0
   */
  bool hasInteractions() const;

  virtual void display() const {}

  /** Display the set of blocks for  a given indexSet
   *
   *  \param  indexSet  the concerned index set
   */
  virtual void displayBlocks(std::shared_ptr<siconos::graphs::InteractionsGraph> indexSet);

  /** compute interactionBlocks if necessary (this depends on the type of
   *  OSNS, on the indexSets ...)
   */
  virtual void updateInteractionBlocks();

  /** compute extra-diagonal interactionBlock-matrix
   *
   *  \param ed an edge descriptor
   */
  virtual void computeInteractionBlock(
      const siconos::graphs::InteractionsGraph::EDescriptor& ed) = 0;

  /** compute diagonal Interaction block
   *
   *  \param vd a vertex descriptor
   */
  virtual void computeDiagonalInteractionBlock(
      const siconos::graphs::InteractionsGraph::VDescriptor& vd) = 0;

  /** \return bool _hasBeenUpdated
   */
  bool hasBeenUpdated() { return _hasBeenUpdated; }

  /** turn activation flag
   *
   *  \param v to set _hasBeenUpdated.
   */
  void setHasBeenUpdated(bool v) { _hasBeenUpdated = v; }

  /**
     initialize the problem (topology and so on)

     \param sim the simulation that owns this OSNSPB
  */
  virtual void initialize(std::shared_ptr<siconos::simulation::Simulation> sim);

  /** prepare data of the osns for solving
   *
   *  \param time the current time
   *  \return true if the computation of the OSNS has to be carry on, false otherwise
   */
  virtual bool preCompute(double time) = 0;

  /** To run the solver for ns problem
   *
   *  \param time  current time
   *  \return int information about the solver convergence.
   */
  virtual int compute(double time) = 0;

  /** post treatment for output of the solver
   */
  virtual void postCompute() = 0;

  /**
      change the solver type and its default parameters

      - clear memory for the existing options set
      - create and initialize a new one

      \param solverId the new solver
   */
  void setSolverId(int solverId);

  /** \return the LU factorization of the OSI-related matrices used to compute the current
     InteractionBlock

      \param osi the integrator
      \param ds the concerned dynamical system
  */
  std::shared_ptr<siconos::algebra::SiconosDenseLUMatrix> getOSIMatrix(
      siconos::integrators::OneStepIntegrator& osi,
      std::shared_ptr<siconos::modeling::DynamicalSystem> ds);
};
}  // namespace siconos::nonsmooth_formulations
#endif  // ONESTEPNSPROBLEM_H
