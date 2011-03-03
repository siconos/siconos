/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
/*! \file OneStepNSProblem.h
  \brief Interface to formalize and solve Non-Sooth problems
 */

#ifndef ONESTEPNSPROBLEM_H
#define ONESTEPNSPROBLEM_H

#include "InteractionsSet.hpp"
#include "SimulationTypeDef.hpp"
#include "NumericsOptions.h"
#include "NumericsMatrix.h"
#include "OSNSMatrix.hpp"

class Simulation;
class DynamicalSystem;
class UnitaryRelation;
class SiconosMatrix;
TYPEDEF_SPTR(NumericsOptions);
TYPEDEF_SPTR(SolverOptions);

/** default name for the OneStepNSProblem of the simulation */
const std::string DEFAULT_OSNS_NAME = "unamed";

/** Non Smooth Problem Formalization and Simulation
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * This is an abstract class, that provides an interface to define a
 * non smooth problem:
 * -> a formulation (ie the way the problem is written)
 * -> a solver (algorithm and solving formulation, that can be
        different from problem formulation)
 * -> routines to compute the problem solution.
 *
 * The available problem formulation, given by derived classes, are:
 *  - LCP
 *  - FrictionContact2D and 3D
 *  - QP
 *  - Relay
 *
 *  See Solver class or Numerics documentation for details on
 *  algorithm name and parameters.
 *
 *  Note: simulation is a required input for construction of a
 *  OneStepNSProblem.
 *
 * \section osns_options Options for Numerics and the driver for solvers
 *
 *  When the Numerics driver is called, two input arguments are
 *  required to set specific options:
 *  - one for general options in Numerics (verbose mode ...)
 * - the other to set the solver options (name, tolerance, max. number
      of iterations ...)
 *
 *  The second one is a member of the NonSmoothSolver, and thus filled
 *  during its construction. The general options are set thanks to
 *  specific functions (only setNumericsVerboseMode() at the time).
 *
 *
 * Remark: UnitaryMatrixRowIterator will be used to iterate through
 *  what corresponds to rows of unitaryBlocks (a row for a
 *  UnitaryRelation, named URrow) and for a row,
 *  UnitaryMatrixColumnIterator will be used to iterate through
 *  columns, ie through all the UnitaryRelations that are linked to
 *  URrow.
 *
 */
class OneStepNSProblem
{

protected:

  /** Algorithm/solver name
  std::string _numerics_solver_name;*/

  int _numerics_solver_id;
  /** Numerics structure used to solve solver options */
  SP::SolverOptions _numerics_solver_options;

  /** id/name of the problem */
  std::string _id;

  /** size of the problem to solve */
  unsigned int _sizeOutput;

  /** map that links each DynamicalSystem with the corresponding DSBlocks
      map < SP::DynamicalSystem , SiconosMatrix * > */
  MapOfDSMatrices _DSBlocks;



  /** link to the simulation that owns the NSPb */
  SP::Simulation _simulation;

  /** the XML object linked to the OneStepNSProblem to read XML data */
  SP::OneStepNSProblemXML _onestepnspbxml;

  /** set of Interactions: link to the Interactions of the Non Smooth
   * Dynamical System Note: no get or set functions for this object in
   * the class -> used only in OneStepNSProblem methods. */
  SP::InteractionsSet _OSNSInteractions;

  /** minimum index set number to be taken into account */
  unsigned int _levelMin;

  /** minimum index set number to be taken into account - For example,
   * if level_min = 1 and level_max = 2, first and second derivatives
   * of y and lambda will be taken into account in the non-smooth
   * problem. Usually, level_max depends only on the non-smooth law
   * and is given by the relative degree, whereas level_min can also
   * depends on the integrator. These values are computed by the
   * Simulation and given as input arguments for the OneStepNSProblem.
   * Classical values are (0,0) for electrical (degree 0) systems,
   * (1,1) for mechanical ones (degree 2).
   */
  unsigned int _levelMax;

  /** maximum value for sizeOutput. Set to the number of declared
      constraints by default (topology->getNumberOfConstraints());
      This value is used to allocate memory for M during initialize
      call. The best choice is to set maxSize to the estimated maximum
      dimension of the problem. It must not exceed ...
  */
  unsigned int _maxSize;

  /** Timer: cpu time spent in solver */
  clock_t _CPUtime;

  /** Number of calls to the solver */
  unsigned int _nbIter;

  /** Numerics (C) structure used to define global options for
      Numerics functions calls */
  SP::NumericsOptions _numerics_options;

  /*During Newton it, this flag allow to update the numerics matrices only once if necessary.*/
  bool _hasBeUpdated;


  // --- CONSTRUCTORS/DESTRUCTOR ---

  /** default constructor
   */
  OneStepNSProblem();

private:

  /** copy constructor (private => no copy nor pass-by value)
   */
  OneStepNSProblem(const OneStepNSProblem&) {};

  /** assignment (private => forbidden) */
  OneStepNSProblem& operator=(const OneStepNSProblem&);

public:
  OneStepNSProblem(const int newNumericsSolverId);
  /** depressed xml constructor
   *  \param string: problem type
   *  \param SP::OneStepNSProblemXML : the XML linked-object
   */
  OneStepNSProblem(const std::string&, SP::OneStepNSProblemXML);
  /**  xml constructor
   *  \param SP::OneStepNSProblemXML : the XML linked-object
   */
  OneStepNSProblem(SP::OneStepNSProblemXML);
  /** deressed constructor from data
   *  \param string: problem type
   *  \param string : id
   *  \param int : solver identifier
   */
  OneStepNSProblem(const std::string&, const std::string&, const int);
  /**  constructor from data
   *  \param string : id
   *  \param string : solver identifier
   */
  OneStepNSProblem(const std::string&, const int);
  /** destructor
   */
  virtual ~OneStepNSProblem();

  // --- GETTERS/SETTERS ---


  /** To get the SolverOptions structure
   *  \return , the numerics structure used to save solver parameters
   */
  inline SP::SolverOptions numericsSolverOptions() const
  {
    return _numerics_solver_options;
  };

  /** to get the id of the OneStepNSProblem
   *  \return string
   */
  inline std::string getId() const
  {
    return _id;
  }

  /** set the id of the OneStepNSProblem
   *  \param: string
   */
  inline void setId(const std::string& newVal)
  {
    _id = newVal;
  }

  /** get dimension of the problem
   *  \return an unsigned ing
   */
  inline unsigned int getSizeOutput() const
  {
    return _sizeOutput;
  }

  /** set the value of sizeOutput
   *  \param an unsigned int
   */
  inline void setSizeOutput(const unsigned int newVal)
  {
    _sizeOutput = newVal;
  }

  /** get the DSBlocks matrices map
   *  \return a MapOfDSMatrices
   */
  inline const MapOfDSMatrices getDSBlocks() const
  {
    return _DSBlocks;
  };

  /** get the DSBlock orresponding to DS1
   *  \param a pointer to DynamicalSystem, DS1
   *  \return a pointer to SiconosMatrix
   */
  SP::SiconosMatrix dSBlock(SP::DynamicalSystem) const ;

  /** set the map of DS matrices
   *  \param a MapOfDSMatrices
   */
  void setDSBlocks(const MapOfDSMatrices&);



  /** get the Simulation
   *  \return a pointer on Simulation
   */
  inline SP::Simulation simulation() const
  {
    return _simulation;
  }

  /** set the Simulation of the OneStepNSProblem
   *  \param a pointer to Simulation
   */
  inline void setSimulationPtr(SP::Simulation newS)
  {
    _simulation = newS;
  }

  /** get the Interactions set
   *  \return an InteractionsSet
   */
  inline SP::InteractionsSet interactions() const
  {
    return _OSNSInteractions;
  }

  /** get level min value
   *  \return an unsigned int
   */
  inline unsigned int levelMin() const
  {
    return _levelMin;
  }

  /** set the value of level min
   *  \param an unsigned int
   */
  inline void setLevelMin(unsigned int newVal)
  {
    _levelMin = newVal;
  }

  /** get level max value
   *  \return an unsigned int
   */
  inline unsigned int getLevelMax() const
  {
    return _levelMax;
  }

  /** set the value of level  max
   *  \param an unsigned int
   */
  inline void setLevelMax(unsigned int newVal)
  {
    _levelMax = newVal;
  }

  /** set the values of level min and max
   *  \param an unsigned int (levelMin value)
   *  \param an unsigned int (levelMax value)
   */
  inline void setLevels(unsigned int newMin, unsigned int newMax)
  {
    _levelMin = newMin;
    _levelMax = newMax;
  }

  /** get maximum value allowed for the dimension of the problem
   *  \return an unsigned int
   */
  inline unsigned int maxSize() const
  {
    return _maxSize;
  }

  /** set the value of maxSize
   *  \param an unsigned int
   */
  inline void setMaxSize(const unsigned int newVal)
  {
    _maxSize = newVal;
  }

  /** get the total (CPU) time spent in the solver
   *  \return: a double
   */
  inline double getCPUtime() const
  {
    return _CPUtime / (double)CLOCKS_PER_SEC;
  };

  /** get the number of call to ns solver
   *  \return: an unsigned int
   */
  inline unsigned int getNumberOfIterations() const
  {
    return _nbIter;
  };

  /** set Numerics verbose mode
      \param a boolean = 1 for on, = 0 for off
   */
  inline void setNumericsVerboseMode(bool vMode)
  {
    _numerics_options->verboseMode = vMode;
  };

  /** reset stat (nbIter and CPUtime)
   */
  inline void resetStat()
  {
    _CPUtime = 0;
    _nbIter = 0;
  };

  /** display stat. info (CPU time and nb of iterations achieved)
   */
  void printStat();

  virtual void display() const
  {
    ;
  }

  /** compute unitaryBlocks if necessary (this depends on the type of
   * OSNS, on the indexSets ...)
   */
  virtual void updateUnitaryBlocks();

  /** computes all diagonal and extra-diagonal unitaryBlock-matrices
   *  useless ?
   */
  virtual void computeAllUnitaryBlocks();

  /** compute extra-diagonal unitaryBlock-matrix
   *  \param an edge descriptor
   */
  virtual void computeUnitaryBlock(const UnitaryRelationsGraph::EDescriptor&) PURE_DEF;

  /** compute diagonal unitary block
   * \param a vertex descriptor
   */
  virtual void computeDiagonalUnitaryBlock(const UnitaryRelationsGraph::VDescriptor&) PURE_DEF;

  /** compute DSBlocks if necessary (this depends on the type of
      OSNS, on the indexSets ...)
  */
  void updateDSBlocks();

  /** computes all  DSBlock-matrices
   */
  void computeAllDSBlocks();

  /**
   * return _hasBeUpdated
   */
  bool hasBeUpdated()
  {
    return _hasBeUpdated;
  }
  /**
   * to set _hasBeUpdated.
   */
  void setHasBeUpdated(bool v)
  {
    _hasBeUpdated = v;
  }

  /** computes DSBlock-matrix that corresponds to DS1
   *  Move this to Unitary Relation class?
   *  \param a pointer to DynamicalSystem DS1
   */
  virtual void computeDSBlock(SP::DynamicalSystem);




  /** initialize the problem(compute topology ...)
      \param the simulation, owner of this OSNSPB
    */
  virtual void initialize(SP::Simulation);

  /** save Interactions states in Memory, called to save the current state of the Newton iteration.
   */
  virtual void saveInMemory();

  /** save y_k, called by TimeDiscretisation::process.
   */
  virtual void saveTimeStepInMemory();
  /** prepare data of the osns for solving
   *  param double : current time
   */
  virtual void preCompute(double) = 0;

  /** To run the solver for ns problem
   *  \param double : current time
   *  \return int, information about the solver convergence.
   */
  virtual int compute(double) = 0;

  /** post treatment for output of the solver
   */
  virtual void postCompute() = 0;

  /** copy the data of the OneStepNSProblem to the XML tree
   */
  virtual void saveNSProblemToXML() = 0;

  /** get the OSI-related matrices used to compute the current Unitary
   * Relation unitaryBlock (Ex: for Moreau, W)
   *  \param a pointer to UnitaryRelation
   *  \param a MapOfDSMatrices(in-out parameter)
   */
  virtual void getOSIMaps(SP::UnitaryRelation, MapOfDSMatrices&);

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();

  /** clear associated maps */
  void clear();

};

#endif // ONESTEPNSPROBLEM_H
