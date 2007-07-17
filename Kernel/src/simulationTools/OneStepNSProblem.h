/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
/*! \file
*/

#ifndef ONESTEPNSPROBLEM_H
#define ONESTEPNSPROBLEM_H

#include "InteractionsSet.h"
#include "RuntimeException.h"
#include "Solver.h"

const std::string DEFAULT_OSNSPB = "LCP";

class Simulation;
class SiconosVector;
class DynamicalSystem;
class UnitaryRelation;
class InteractionsSet;
class EqualityConstraint;
class SiconosMatrix;
class OneStepNSProblemXML;
class Solver;

/** map of SiconosMatrices with a UnitaryRelations as a key - Used for diagonal block-terms in assembled matrices of LCP etc ...*/
typedef std::map< UnitaryRelation* , SiconosMatrix*>  MapOfUnitaryMatrices;

/** corresponding iterator */
typedef MapOfUnitaryMatrices::iterator UnitaryMatrixColumnIterator ;
typedef MapOfUnitaryMatrices::const_iterator ConstUnitaryMatrixColumnIterator;

/** map of MapOfUnitaryMatrices with a UnitaryRelation as a key - Used for extra-diagonal block-terms in assembled matrices of LCP etc ..*/
typedef std::map< UnitaryRelation* , MapOfUnitaryMatrices >  MapOfMapOfUnitaryMatrices;

/** corresponding iterators */
typedef MapOfMapOfUnitaryMatrices::iterator UnitaryMatrixRowIterator ;
typedef MapOfMapOfUnitaryMatrices::const_iterator ConstUnitaryMatrixRowIterator ;

/** map of map of bools, with UnitaryRelations as keys */
typedef std::map< UnitaryRelation* , std::map<UnitaryRelation*, bool> >  MapOfMapOfBool;

// Remark: UnitaryMatrixRowIterator will be used to iterate through what corresponds to rows of blocks (a row for a UnitaryRelation, named URrow) and for
// a row, UnitaryMatrixColumnIterator will be used to iterate through columns, ie through all the UnitaryRelations that are linked to URrow.

/** map of SiconosMatrix; key = the related DS*/
typedef std::map<DynamicalSystem*, SiconosMatrix*> MapOfMatrices;

/** map of double; key = the related DS */
typedef std::map<DynamicalSystem*, double> MapOfDouble;

/** default name for One Step NS Problem of the simulation */
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

  /** map that links each UnitaryRelation with an int that gives the position of the corresponding block matrix
   *  in the full matrix (M in LCP case) in number of single elements
   */
  std::map< UnitaryRelation* , unsigned int > blocksPositions;

  /** map that links each UnitaryRelation with an int that gives the position of the corresponding block matrix
   *  in the full matrix (M in LCP case) in number of blocks (for use with block matrix storage)
   */
  std::map< UnitaryRelation* , unsigned int > blocksIndexes;

  /** Solver for Non Smooth Problem*/
  Solver* solver;

  /** bool to check whether solver has been allocated inside the class or not */
  bool isSolverAllocatedIn;

  /** link to the simulation that owns the NSPb */
  Simulation *simulation;

  /** the XML object linked to the OneStepNSProblem to read XML data */
  OneStepNSProblemXML* onestepnspbxml;

  /** set of Interactions: link to the Interactions of the Non Smooth Dynamical System
    * Note: no get or set functions for this object in the class -> used only in OneStepNSProblem methods. */
  InteractionsSet OSNSInteractions;

  /** minimum index set number to be taken into account */
  unsigned int levelMin;

  /** minimum index set number to be taken into account - For example, if level_min = 1 and level_max = 2, first and second derivatives of y and lambda will be
   * taken into account in the non-smooth problem. Usually, level_max depends only on the non-smooth law and is given by the relative degree, whereas level_min
   * can also depends on the integrator. These values are computed by the Simulation and given as input arguments for the OneStepNSProblem.
   * Classical values are (0,0) for electrical (degree 0) systems, (1,1) for mechanical ones (degree 2).
   */
  unsigned int levelMax;

  // --- CONSTRUCTORS/DESTRUCTOR ---

  /** default constructor
  *  \param string: problem type
  */
  OneStepNSProblem(const std::string = DEFAULT_OSNSPB);

public:

  /** xml constructor
  *  \param string: problem type
  *  \param OneStepNSProblemXML* : the XML linked-object
  *  \param Simulation *: the simulation that owns the problem
  */
  OneStepNSProblem(const std::string, OneStepNSProblemXML*, Simulation *);

  /** constructor from data
  *  \param string: problem type
  *  \param Simulation *: the simulation that owns this problem
  *  \param string : id
  *  \param Solver *: pointer on object that contains solver algorithm definition (optional)
  */
  OneStepNSProblem(const std::string, Simulation *, const std::string, Solver* = NULL);

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
  inline void setType(const std::string  newVal)
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
  inline void setId(const std::string newVal)
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

  /** get the Solver
  *  \return a pointer on Solver
  */
  inline Solver* getSolverPtr() const
  {
    return solver;
  }

  /** set the Solver of the OneStepNSProblem
  *  \param: a pointer on Solver
  */
  void setSolverPtr(Solver*);

  /** get the Simulation
  *  \return a pointer on Simulation
  */
  inline Simulation* getSimulationPtr() const
  {
    return simulation;
  }

  /** get the OneStepNSProblemXML
  *  \return a pointer on OneStepNSProblemXML
  */
  inline OneStepNSProblemXML* getOneStepNSProblemXML() const
  {
    return onestepnspbxml;
  }

  /** get the Interactions set
  *  \return an InteractionsSet
  */
  inline InteractionsSet getInteractions() const
  {
    return OSNSInteractions;
  }

  /** set the OneStepNSProblemXML
  *  \param a pointer on OneStepNSProblemXML
  */
  inline void setOneStepNSProblemXML(OneStepNSProblemXML* osnspb)
  {
    onestepnspbxml = osnspb;
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
  inline void setLevelMin(const unsigned int newVal)
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
  inline void setLevelMax(const unsigned int newVal)
  {
    levelMax = newVal;
  }

  /** set the values of level min and max
  *  \param an unsigned int (levelMin value)
  *  \param an unsigned int (levelMax value)
  */
  inline void setLevels(const unsigned int newMin, const unsigned int newMax)
  {
    levelMin = newMin;
    levelMax = newMax;
  }

  // --- OTHER FUNCTIONS ---

  /** fills in blocksPositions map, ie computes variables blocks positions in the full matrix (M in LCP case ...)
  */
  void computeUnitaryRelationsPositions();

  /** compute SizeOutput, ie count the total number
  *  of relations constrained
  */
  void computeSizeOutput();

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
  virtual void preCompute(const double) = 0;

  /** make the computation so solve the NS problem
  *  param double : current time
  */
  virtual void compute(const double) = 0;

  /** post treatment for output of the solver
   */
  virtual void postCompute() = 0;

  /** copy the data of the OneStepNSProblem to the XML tree
  *  \exception RuntimeException
  */
  virtual void saveNSProblemToXML() = 0;

  /** return exception and message if solver failed
  *  \param: output from solve_... (Numerics routine)
  */
  void check_solver(const int) const;

  /** get the OSI-related matrices used to compute the current Unitary Relation block (Ex: for Moreau, W and Theta)
  *  \param a pointer to UnitaryRelation
  *  \param a MapOfMatrices(in-out parameter)
  *  \param a MapOfDouble(in-out parameter)
  */
  virtual void getOSIMaps(UnitaryRelation*, MapOfMatrices&, MapOfDouble&);

};

#endif // ONESTEPNSPROBLEM_H
