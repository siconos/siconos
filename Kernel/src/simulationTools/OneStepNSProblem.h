/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
/*! \file OneStepNSProblem.h
  \brief Interface to formalize and solve Non-Sooth problems
 */

#ifndef ONESTEPNSPROBLEM_H
#define ONESTEPNSPROBLEM_H

#include "InteractionsSet.h"
#include "SimulationTypeDef.h"
#include "Numerics_Options.h"
#include "NumericsMatrix.h"
#include "NonSmoothSolver.h"
#include "OSNSMatrix.h"

class Simulation;
class DynamicalSystem;
class UnitaryRelation;
class SiconosMatrix;
class OneStepNSProblemXML;

/** default name for the OneStepNSProblem of the simulation */
const std::string DEFAULT_OSNS_NAME = "unamed";

/** Non Smooth Problem Formalization and Simulation
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) Apr 26, 2004
 *
 * This is an abstract class, that provides an interface to define a non smooth problem:
 *   -> a formulation (ie the way the problem is written)
 *   -> a solver (algorithm and solving formulation, that can be different from problem formulation)
 *   -> routines to compute the problem solution.
 *
 * The available problem formulation, given by derived classes, are:
 *  - LCP
 *  - FrictionContact2D and 3D
 *  - QP
 *  - Relay
 *
 *  See Solver class or Numerics documentation for details on algorithm name and parameters.
 *
 *  Note: simulation is a required input for construction of a OneStepNSProblem.
 *
 * \section osns_options Options for Numerics and the driver for solvers
 *
 *  When the Numerics driver is called, two input arguments are required to set specific options:
 *  - one for general options in Numerics (verbose mode ...)
 *  - the other to set the solver options (name, tolerance, max. number of iterations ...)
 *
 *  The second one is a member of the NonSmoothSolver, and thus filled during its construction. \n
 *  The general options are set thanks to specific functions (only setNumericsVerboseMode() at the time).
 *
 *
 *
 * Remark: UnitaryMatrixRowIterator will be used to iterate through what corresponds to rows of blocks (a row for a UnitaryRelation, named URrow) and for
 *  a row, UnitaryMatrixColumnIterator will be used to iterate through columns, ie through all the UnitaryRelations that are linked to URrow.
 *
 */
class OneStepNSProblem
{

protected:

  /** type of the OneStepNSProblem (LCP ...) */
  std::string nspbType;

  /** id/name of the problem */
  std::string id;

  /** size of the problem to solve */
  unsigned int sizeOutput;

  /** map that links each UnitaryRelation with the corresponding blocks
      map < UnitaryRelationA * , map < UnitaryRelationB* , blockMatrixAB > >
      UnitaryRelations A and B are coupled through blockMatrixAB.  */
  MapOfMapOfUnitaryMatrices blocks;

  /** inside-class allocation flags. To each couple of Unitary Relations corresponds a bool, true if
      the block has been allocated in OneStepNSProblem, else false. */
  MapOfMapOfBool isBlockAllocatedIn;

  /** Solver for Non Smooth Problem*/
  NonSmoothSolver* solver;

  /** bool to check whether solver has been allocated inside the class or not */
  bool isSolverAllocatedIn;

  /** link to the simulation that owns the NSPb */
  Simulation *simulation;

  /** the XML object linked to the OneStepNSProblem to read XML data */
  OneStepNSProblemXML* onestepnspbxml;

  /** set of Interactions: link to the Interactions of the Non Smooth Dynamical System
   * Note: no get or set functions for this object in the class -> used only in OneStepNSProblem methods. */
  InteractionsSet * OSNSInteractions;

  /** minimum index set number to be taken into account */
  unsigned int levelMin;

  /** minimum index set number to be taken into account - For example, if level_min = 1 and level_max = 2, first and second derivatives of y and lambda will be
   * taken into account in the non-smooth problem. Usually, level_max depends only on the non-smooth law and is given by the relative degree, whereas level_min
   * can also depends on the integrator. These values are computed by the Simulation and given as input arguments for the OneStepNSProblem.
   * Classical values are (0,0) for electrical (degree 0) systems, (1,1) for mechanical ones (degree 2).
   */
  unsigned int levelMax;

  /** maximum value for sizeOutput. Set to the number of declared constraints by default (topology->getNumberOfConstraints());
      This value is used to allocate memory for M during initialize call. The best choice is to set maxSize to the estimated maximum
      dimension of the problem. It must not exceed ...
  */
  unsigned int maxSize;

  /** Timer: cpu time spent in solver */
  clock_t CPUtime;

  /** Number of calls to the solver */
  unsigned int nbIter;

  /** Numerics (C) structure used to define global options for Numerics functions calls */
  Numerics_Options * numerics_options;

  // --- CONSTRUCTORS/DESTRUCTOR ---

private:

  /** default constructor
   */
  OneStepNSProblem();

  /** copy constructor (private => no copy nor pass-by value)
   */
  OneStepNSProblem(const OneStepNSProblem&);

public:

  /** xml constructor
   *  \param string: problem type
   *  \param OneStepNSProblemXML* : the XML linked-object
   *  \param Simulation *: the simulation that owns the problem
   */
  OneStepNSProblem(const std::string&, OneStepNSProblemXML*, Simulation *);

  /** constructor from data
   *  \param string: problem type
   *  \param Simulation *: the simulation that owns this problem
   *  \param string : id
   *  \param Solver *: pointer on object that contains solver algorithm definition (optional)
   */
  OneStepNSProblem(const std::string&, Simulation *, const std::string&, NonSmoothSolver* = NULL);

  /** destructor
   */
  virtual ~OneStepNSProblem();

  // --- GETTERS/SETTERS ---

  /** to get the type of the OneStepNSProblem
   *  \return string
   */
  inline std::string getType() const
  {
    return nspbType;
  }

  /** set the type of the OneStepNSProblem
   *  \param: string
   */
  inline void setType(const std::string&  newVal)
  {
    nspbType = newVal;
  }

  /** to get the id of the OneStepNSProblem
   *  \return string
   */
  inline std::string getId() const
  {
    return id;
  }

  /** set the id of the OneStepNSProblem
   *  \param: string
   */
  inline void setId(const std::string& newVal)
  {
    id = newVal;
  }

  /** get dimension of the problem
   *  \return an unsigned ing
   */
  inline const unsigned int getSizeOutput() const
  {
    return sizeOutput;
  }

  /** set the value of sizeOutput
   *  \param an unsigned int
   */
  inline void setSizeOutput(const unsigned int newVal)
  {
    sizeOutput = newVal;
  }

  /** get the blocks matrices map
   *  \return a MapOfMapOfUnitaryMatrices
   */
  inline const MapOfMapOfUnitaryMatrices getBlocks() const
  {
    return blocks;
  };

  /** get the block orresponding to UR1 and UR2
   *  \param a pointer to UnitaryRelation, UR1
   *  \param a pointer to UnitaryRelation, optional, default value = NULL, in that case UR2 = UR1 (ie get "diagonal" block)
   *  \return a pointer to SiconosMatrix
   */
  SiconosMatrix* getBlockPtr(UnitaryRelation*, UnitaryRelation* = NULL) const ;

  /** set the map of unitary matrices
   *  \param a MapOfMapOfUnitaryMatrices
   */
  void setBlocks(const MapOfMapOfUnitaryMatrices&);

  /** clear the map of blocks (ie release memory)
   */
  void clearBlocks();

  /** get the NonSmoothSolver
   *  \return a pointer on NonSmoothSolver
   */
  inline NonSmoothSolver* getNonSmoothSolverPtr() const
  {
    return solver;
  }

  /** set the NonSmoothSolver of the OneStepNSProblem
   *  \param: a pointer on NonSmoothSolver
   */
  void setNonSmoothSolverPtr(NonSmoothSolver*);

  /** get the Simulation
   *  \return a pointer on Simulation
   */
  inline Simulation* getSimulationPtr() const
  {
    return simulation;
  }

  /** get the Interactions set
   *  \return an InteractionsSet
   */
  inline InteractionsSet * getInteractions() const
  {
    return OSNSInteractions;
  }

  /** get level min value
   *  \return an unsigned int
   */
  inline const unsigned int getLevelMin() const
  {
    return levelMin;
  }

  /** set the value of level min
   *  \param an unsigned int
   */
  inline void setLevelMin(unsigned int newVal)
  {
    levelMin = newVal;
  }

  /** get level max value
   *  \return an unsigned int
   */
  inline const unsigned int getLevelMax() const
  {
    return levelMax;
  }

  /** set the value of level  max
   *  \param an unsigned int
   */
  inline void setLevelMax(unsigned int newVal)
  {
    levelMax = newVal;
  }

  /** set the values of level min and max
   *  \param an unsigned int (levelMin value)
   *  \param an unsigned int (levelMax value)
   */
  inline void setLevels(unsigned int newMin, unsigned int newMax)
  {
    levelMin = newMin;
    levelMax = newMax;
  }

  /** get maximum value allowed for the dimension of the problem
   *  \return an unsigned int
   */
  inline const unsigned int getMaxSize() const
  {
    return maxSize;
  }

  /** set the value of maxSize
   *  \param an unsigned int
   */
  inline void setMaxSize(const unsigned int newVal)
  {
    maxSize = newVal;
  }

  /** get the total (CPU) time spent in the solver
   *  \return: a double
   */
  inline const double getCPUtime() const
  {
    return CPUtime / (double)CLOCKS_PER_SEC;
  };

  /** get the number of call to ns solver
   *  \return: an unsigned int
   */
  inline const unsigned int getNumberOfIterations() const
  {
    return nbIter;
  };

  /** set Numerics verbose mode
      \param a boolean = 1 for on, = 0 for off
   */
  inline void setNumericsVerboseMode(bool vMode)
  {
    numerics_options->verboseMode = vMode;
  };

  /** reset stat (nbIter and CPUtime)
   */
  inline void resetStat()
  {
    CPUtime = 0;
    nbIter = 0;
  };

  /** display stat. info (CPU time and nb of iterations achieved)
   */
  void printStat();

  /** compute blocks if necessary (this depends on the type of OSNS, on the indexSets ...)
   */
  void updateBlocks();

  /** computes all diagonal and extra-diagonal block-matrices
   */
  void computeAllBlocks();

  /** computes extra diagonal block-matrix that corresponds to UR1 and UR2
   *  Move this to Unitary Relation class?
   *  \param a pointer to UnitaryRelation
   *  \param a pointer to UnitaryRelation
   */
  virtual void computeBlock(UnitaryRelation*, UnitaryRelation*);

  /** initialize the problem(compute topology ...)
   */
  virtual void initialize();

  /** save Interactions states in Memory
   */
  void saveInMemory();

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

  /** get the OSI-related matrices used to compute the current Unitary Relation block (Ex: for Moreau, W and Theta)
   *  \param a pointer to UnitaryRelation
   *  \param a MapOfMatrices(in-out parameter)
   *  \param a MapOfDouble(in-out parameter)
   */
  virtual void getOSIMaps(UnitaryRelation*, MapOfMatrices&, MapOfDouble&);

};

#endif // ONESTEPNSPROBLEM_H
