/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
/*! \file OneStepNSProblem.hpp
  \brief Interface to formalize and solve Non-Sooth problems
 */

#ifndef ONESTEPNSPROBLEM_H
#define ONESTEPNSPROBLEM_H

#include "SiconosFwd.hpp"
#include "SimulationTypeDef.hpp"
#include "SimulationGraphs.hpp"

/** Non Smooth Problem Formalization and Simulation

   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) Apr 26, 2004

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
   - OSNSMultipleImpact
   - GlobalFrictionContact (unstable)

   The usual way to build and initialize a one-step nonsmooth problem is :
   - call constructor with the id of the required Numerics solver.
   (see Solver class or Numerics documentation for details on algorithm name and parameters).
   - initialize(simulation)
   Initialize process is usually done through model->initialize(simulation). 
   See Examples for practical details.

   \section osns_options Options for Numerics and the driver for solvers

   When the Numerics driver is called, two input arguments are
   required to set specific options:
   - the global Numerics options (verbose mode ...) --> NumericsOptions
   - the solver options (name, tolerance, max. number of iterations ...) --> _SolverOptions, \ref NumericsSolver.

   Default values are always set in solver options the OneStepNSProblem is built
   but if you need to set them yourself, please see \ref NumericsSolver. 

 */
class OneStepNSProblem
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(OneStepNSProblem);

  /** Numerics solver id */
  int _numerics_solver_id;

  /** Numerics structure used to solve solver options */
  SP::SolverOptions _numerics_solver_options;

  /** size of the problem to solve */
  unsigned int _sizeOutput;

  /** link to the simulation that owns the NSPb */
  SP::Simulation _simulation;

  /** level of index sets that is considered by this osnsp */
  unsigned int _indexSetLevel;

  /** level of input and output variables ofs osnsp.
   *  We consider that the osnsp computes y[_inputOutputLevel] and lambda[_inputOutputLevel]
   */
  unsigned int _inputOutputLevel;

  /** maximum value for sizeOutput. Set to the number of declared
      constraints by default (topology->getNumberOfConstraints());
      This value is used to allocate memory for M during initialize
      call. The best choice is to set maxSize to the estimated maximum
      dimension of the problem. It must not exceed ...
  */
  unsigned int _maxSize;

  /** Number of calls to the solver */
  unsigned int _nbIter;

  /** Numerics (C) structure used to define global options for
      Numerics functions calls */
  SP::NumericsOptions _numerics_options;

  /*During Newton it, this flag allows to update the numerics matrices only once if necessary.*/
  bool _hasBeenUpdated;

  // --- CONSTRUCTORS/DESTRUCTOR ---
  /** default constructor
   */
  OneStepNSProblem();

private:

  /** copy constructor (private => no copy nor pass-by value)
   */
  OneStepNSProblem(const OneStepNSProblem&) {};

  /** assignment (private => forbidden) 
   * \param osnsp
   * \return OneStepNSProblem&  osnsp
   */
  OneStepNSProblem& operator=(const OneStepNSProblem& osnsp);

public:
  /**  constructor with a solver from Numerics
   *  \param numericsSolverId id of numerics solver, see Numerics for the meaning
   */
  OneStepNSProblem(int numericsSolverId);

  /** destructor
   */
  virtual ~OneStepNSProblem(){};

  // --- GETTERS/SETTERS ---


  /** To get the SolverOptions structure
   *  \return , the numerics structure used to save solver parameters
   */
  inline SP::SolverOptions numericsSolverOptions() const
  {
    return _numerics_solver_options;
  };

  /** To get the NumericsOptions structure
   *  \return , the numerics structure used to save solver parameters
   */
  inline SP::NumericsOptions numericsOptions() const
  {
    return _numerics_options;
  };

  /** get dimension of the problem
   *  \return an unsigned ing
   */
  inline unsigned int getSizeOutput() const
  {
    return _sizeOutput;
  }

  /** set the value of sizeOutput
   *  \param newVal an unsigned int
   */
  inline void setSizeOutput(const unsigned int newVal)
  {
    _sizeOutput = newVal;
  }

  /** get the Simulation
   *  \return a pointer on Simulation
   */
  inline SP::Simulation simulation() const
  {
    return _simulation;
  }

  /** set the Simulation of the OneStepNSProblem
   *  \param newS a pointer to Simulation
   */
  inline void setSimulationPtr(SP::Simulation newS)
  {
    _simulation = newS;
  }

  /** get indexSetLevel
   *  \return an unsigned int
   */
  inline unsigned int indexSetLevel() const
  {
    return _indexSetLevel;
  }

  /** set the value of level min
   *  \param newVal an unsigned int
   */
  inline void setIndexSetLevel(unsigned int newVal)
  {
    _indexSetLevel = newVal;
  }



  /** get the Input/Output level
   *  \return an unsigned int
   */
  inline unsigned int inputOutputLevel() const
  {
    return _inputOutputLevel;
  }

  /** set the value of Input/Output level
   *  \param newVal an unsigned int
   */
  inline void setInputOutputLevel(unsigned int newVal)
  {
    _inputOutputLevel = newVal;
  }

  /** get maximum value allowed for the dimension of the problem
   *  \return an unsigned int
   */
  inline unsigned int maxSize() const
  {
    return _maxSize;
  }

  /** set the value of maxSize
   *  \param newVal an unsigned int
   */
  inline void setMaxSize(const unsigned int newVal)
  {
    _maxSize = newVal;
  }

  /** get the number of call to ns solver
   *  \return: an unsigned int
   */
  inline unsigned int getNumberOfIterations() const
  {
    return _nbIter;
  };

  /** set Numerics verbose mode
      \param vMode 1 for on, 0 for off
   */
  void setNumericsVerboseMode(bool vMode);

  /** reset stat (nbIter and CPUtime)
   */
  inline void resetStat()
  {
    _nbIter = 0;
  };

  /** Check if the OSNSPb has interactions.
      \return bool = true if the  osnsp has interactions, i.e. indexSet(_indexSetLevel)->size >0 
   */
  bool hasInteractions() const;

  /** display stat. info (CPU time and nb of iterations achieved)
   */
  void printStat();

  virtual void display() const
  {
    ;
  }

  /** Display the set of blocks for  a given indexSet
   * \param  indexSet  the concerned index set
   */
  virtual void displayBlocks(SP::InteractionsGraph indexSet);

  /** compute interactionBlocks if necessary (this depends on the type of
   * OSNS, on the indexSets ...)
   */
  virtual void updateInteractionBlocks();

  /** compute extra-diagonal interactionBlock-matrix
   *  \param ed an edge descriptor
   */
  virtual void computeInteractionBlock(const InteractionsGraph::EDescriptor& ed ) = 0;

  /** compute diagonal Interaction block
   * \param vd a vertex descriptor
   */
  virtual void computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd) = 0;

  /**
   * \return bool _hasBeenUpdated
   */
  bool hasBeenUpdated()
  {
    return _hasBeenUpdated;
  }

  /** 
   * \param v to set _hasBeenUpdated.
   */
  void setHasBeenUpdated(bool v)
  {
    _hasBeenUpdated = v;
  }

  /** initialize the problem(compute topology ...)
      \param sim the simulation, owner of this OSNSPB
    */
  virtual void initialize(SP::Simulation sim);

  /** prepare data of the osns for solving
   *  \param time the current time
   *  \return true if the computation of the OSNS has to be carry on, false otherwise
   */
  virtual bool preCompute(double time) = 0;

  /** To run the solver for ns problem
   *  \param time  current time
   *  \return int information about the solver convergence.
   */
  virtual int compute(double time) = 0;

  /** post treatment for output of the solver
   */
  virtual void postCompute() = 0;

  /** change the solver type. This requires a reset of the Solveroption struct
   * \param solverId the new solver
   */
  virtual void setSolverId(int solverId);

  /** get the OSI-related matrices used to compute the current InteractionBlock
      (Ex: for MoreauJeanOSI, W)
      \param osi the OSI of the concerned dynamical system
      \param ds the concerned dynamical system
      \return the required matrix.
  */
  SP::SimpleMatrix getOSIMatrix(SP::OneStepIntegrator osi, SP::DynamicalSystem ds);

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();

};

#endif // ONESTEPNSPROBLEM_H
