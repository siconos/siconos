/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 *  \version 3.0.0.
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
 * Remark: UnitaryMatrixRowIterator will be used to iterate through what corresponds to rows of unitaryBlocks (a row for a UnitaryRelation, named URrow) and for
 *  a row, UnitaryMatrixColumnIterator will be used to iterate through columns, ie through all the UnitaryRelations that are linked to URrow.
 *
 */
class OneStepNSProblem : public boost::enable_shared_from_this<OneStepNSProblem>
{

protected:

  /** type of the OneStepNSProblem (LCP ...) */
  std::string nspbType;

  /** id/name of the problem */
  std::string id;

  /** size of the problem to solve */
  unsigned int sizeOutput;

  /** map that links each UnitaryRelation with the corresponding unitaryBlocks
      map < UnitaryRelationA * , map < UnitaryRelationB* , unitaryBlockMatrixAB > >
      UnitaryRelations A and B are coupled through unitaryBlockMatrixAB.  */
  MapOfMapOfUnitaryMatrices unitaryBlocks;

  /** inside-class allocation flags. To each couple of Unitary Relations corresponds a bool, true if
      the unitaryBlock has been allocated in OneStepNSProblem, else false. */
  MapOfMapOfUnitaryBool isUnitaryBlockAllocatedIn;

  /** map that links each DynamicalSystem with the corresponding DSBlocks
      map < SP::DynamicalSystem , SiconosMatrix * > */
  MapOfDSMatrices DSBlocks;

  /** map that links each UnitaryRelation and DynamicalSystem with the corresponding unitaryDSBlocks
      map < UnitaryRelationA * , map < DynamicalSystemB * , unitaryDSBlockMatrixAB > >
      UnitaryRelation A and DynamicalSystem B are coupled through unitaryDSBlockMatrixAB.  */
  MapOfUnitaryMapOfDSMatrices unitaryDSBlocks;

  /** map that links each DynamicalSystem and UnitaryRelation with the corresponding DSunitaryBlocks
      map < DynamicalSystemA * , map < UnitaryRelationB* , DSunitaryBlockMatrixAB > >
      Dynamical A and UnitaryRelation B are coupled through DSunitaryBlockMatrixAB.  */
  MapOfDSMapOfUnitaryMatrices DSUnitaryBlocks;

  /** Solver for Non Smooth Problem*/
  SP::NonSmoothSolver solver;

  /** link to the simulation that owns the NSPb */
  SP::Simulation simulation;

  /** the XML object linked to the OneStepNSProblem to read XML data */
  SP::OneStepNSProblemXML onestepnspbxml;

  /** set of Interactions: link to the Interactions of the Non Smooth Dynamical System
   * Note: no get or set functions for this object in the class -> used only in OneStepNSProblem methods. */
  SP::InteractionsSet OSNSInteractions;

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
  SP::Numerics_Options numerics_options;

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
   *  \param SP::OneStepNSProblemXML : the XML linked-object
   */
  OneStepNSProblem(const std::string&, SP::OneStepNSProblemXML);

  /** constructor from data
   *  \param string: problem type
   *  \param string : id
   *  \param Solver *: pointer on object that contains solver algorithm definition (optional)
   */
  OneStepNSProblem(const std::string&, const std::string&, SP::NonSmoothSolver = SP::NonSmoothSolver());

  /** destructor
   */
  virtual ~OneStepNSProblem() {};

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

  /** get the unitaryBlocks matrices map
   *  \return a MapOfMapOfUnitaryMatrices
   */
  inline const MapOfMapOfUnitaryMatrices getUnitaryBlocks() const
  {
    return unitaryBlocks;
  };

  /** get the unitaryBlock orresponding to UR1 and UR2
   *  \param a pointer to UnitaryRelation, UR1
   *  \param a pointer to UnitaryRelation, optional, default value = NULL, in that case UR2 = UR1 (ie get "diagonal" unitaryBlock)
   *  \return a pointer to SiconosMatrix
   */
  SP::SiconosMatrix getUnitaryBlockPtr(SP::UnitaryRelation,
                                       SP::UnitaryRelation = SP::UnitaryRelation()) const ;

  /** set the map of unitary matrices
   *  \param a MapOfMapOfUnitaryMatrices
   */
  void setUnitaryBlocks(const MapOfMapOfUnitaryMatrices&);

  /** clear the map of unitaryBlocks (ie release memory)
   */
  void clearUnitaryBlocks();

  /** get the DSBlocks matrices map
   *  \return a MapOfDSMatrices
   */
  inline const MapOfDSMatrices getDSBlocks() const
  {
    return DSBlocks;
  };

  /** get the DSBlock orresponding to DS1
   *  \param a pointer to DynamicalSystem, DS1
   *  \return a pointer to SiconosMatrix
   */
  SP::SiconosMatrix getDSBlockPtr(SP::DynamicalSystem) const ;

  /** set the map of DS matrices
   *  \param a MapOfDSMatrices
   */
  void setDSBlocks(const MapOfDSMatrices&);

  /** clear the map of DSBlocks (ie release memory)
   */
  void clearDSBlocks();


  /** get the unitaryDSBlocks matrices map
   *  \return a MapOfUnitaryMapOfDSMatrices
   */
  inline const MapOfUnitaryMapOfDSMatrices getUnitaryDSBlocks() const
  {
    return unitaryDSBlocks;
  };

  /** get the unitaryDSBlock corresponding to UR1 and DS2
   *  \param a pointer to UnitaryRelation, UR1
   *  \param a pointer to DynamicalSystem DS2
   *  \return a pointer to SiconosMatrix
   */
  SP::SiconosMatrix getUnitaryDSBlockPtr(SP::UnitaryRelation, SP::DynamicalSystem) const ;

  /** set the map of unitaryDS matrices
   *  \param a MapOfUnitaryMapOfDSMatrices
   */
  void setUnitaryDSBlocks(const MapOfUnitaryMapOfDSMatrices&);

  /** clear the map of unitaryBlocks (ie release memory)
   */
  void clearUnitaryDSBlocks();

  /** get the DSunitaryBlocks matrices map
    *  \return a MapOfDSMapOfUnitaryMatrices
    */
  inline const MapOfDSMapOfUnitaryMatrices getDSUnitaryBlocks() const
  {
    return DSUnitaryBlocks;
  };

  /** get the DSunitaryBlock corresponding to DS1 and UR2
   *  \param a pointer to UnitaryRelation, UR2
   *  \param a pointer to DynamicalSystem DS1
   *  \return a pointer to SiconosMatrix
   */
  SP::SiconosMatrix getDSUnitaryBlockPtr(SP::DynamicalSystem, SP::UnitaryRelation) const ;

  /** set the map of DSUnitary matrices
   *  \param a MapOfDSMapOfUnitaryMatrices
   */
  void setDSUnitaryBlocks(const MapOfDSMapOfUnitaryMatrices&);

  /** clear the map of unitaryBlocks (ie release memory)
   */
  void clearDSUnitaryBlocks();

  /** get the NonSmoothSolver
   *  \return a pointer on NonSmoothSolver
   */
  inline SP::NonSmoothSolver getNonSmoothSolverPtr() const
  {
    return solver;
  }

  /** set the NonSmoothSolver of the OneStepNSProblem
   *  \param: a pointer on NonSmoothSolver
   */
  void setNonSmoothSolverPtr(SP::NonSmoothSolver);

  /** get the Simulation
   *  \return a pointer on Simulation
   */
  inline SP::Simulation getSimulationPtr() const
  {
    return simulation;
  }

  /** set the Simulation of the OneStepNSProblem
   *  \param a pointer to Simulation
   */
  inline void setSimulationPtr(SP::Simulation newS)
  {
    simulation = newS;
  }

  /** get the Interactions set
   *  \return an InteractionsSet
   */
  inline SP::InteractionsSet getInteractions() const
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

  /** compute unitaryBlocks if necessary (this depends on the type of OSNS, on the indexSets ...)
   */
  void updateUnitaryBlocks();

  /** computes all diagonal and extra-diagonal unitaryBlock-matrices
   */
  void computeAllUnitaryBlocks();

  /** computes extra diagonal unitaryBlock-matrix that corresponds to UR1 and UR2
   *  Move this to Unitary Relation class?
   *  \param a pointer to UnitaryRelation
   *  \param a pointer to UnitaryRelation
   */
  virtual void computeUnitaryBlock(SP::UnitaryRelation, SP::UnitaryRelation);

  /** compute DSBlocks if necessary (this depends on the type of OSNS, on the indexSets ...)
  */
  void updateDSBlocks();

  /** computes all  DSBlock-matrices
   */
  void computeAllDSBlocks();

  /** computes DSBlock-matrix that corresponds to DS1
   *  Move this to Unitary Relation class?
   *  \param a pointer to DynamicalSystem DS1
   */
  virtual void computeDSBlock(SP::DynamicalSystem);


  /** compute UnitaryDSBlocks if necessary (this depends on the type of OSNS, on the indexSets ...)
  */
  void updateUnitaryDSBlocks();

  /** computes all unitaryDSBlock-matrices
   */
  void computeAllUnitaryDSBlocks();

  /** computes  unitaryBlock-matrix that corresponds to UR1 and DS2
   *  Move this to Unitary Relation class?
   *  \param a pointer to UnitaryRelation UR1
   *  \param a pointer to DynamicalSystems DS2
   */
  virtual void computeUnitaryDSBlock(SP::UnitaryRelation , SP::DynamicalSystem);


  /** compute DSUnitaryBlocks if necessary (this depends on the type of OSNS, on the indexSets ...)
  */
  void updateDSUnitaryBlocks();

  /** computes all DSunitaryBlock-matrices
   */
  void computeAllDSUnitaryBlocks();

  /** computes  DSUnitaryBlock-matrix that corresponds to UR1 and DS2
   *  Move this to Unitary Relation class?
   *  \param a pointer to UnitaryRelation UR1
   *  \param a pointer to DynamicalSystems DS2
   */
  virtual void computeDSUnitaryBlock(SP::DynamicalSystem, SP::UnitaryRelation);


  /** initialize the problem(compute topology ...)
      \param the simulation, owner of this OSNSPB
    */
  virtual void initialize(SP::Simulation);

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

  /** get the OSI-related matrices used to compute the current Unitary Relation unitaryBlock (Ex: for Moreau, W and Theta)
   *  \param a pointer to UnitaryRelation
   *  \param a MapOfDSMatrices(in-out parameter)
   *  \param a MapOfDouble(in-out parameter)
   */
  virtual void getOSIMaps(SP::UnitaryRelation, MapOfDSMatrices&, MapOfDouble&);

};

#endif // ONESTEPNSPROBLEM_H
