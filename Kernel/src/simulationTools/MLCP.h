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
/*! \file MLCP.h
\brief Linear Complementarity Problem formulation and solving
*/

#ifndef MLCP_H
#define MLCP_H

#include "OneStepNSProblem.h"
#include "SimpleVector.h"
#include "SimpleMatrix.h"
#include "SparseBlockMatrix.h"
#include <sys/time.h>

class OneStepNSProblem;

/** Formalization and Resolution of a Mixed Linear Complementarity Problem (MLCP)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * \section MLCPintro Aim of the MLCP class
 *
 * This class is devoted to the formalization and the resolution of the
 * Mixed Linear Complementarity Problem (MLCP) defined by :
 *  * \f[
 * 0 =  Au + Cv + a
 * \f]
 * \f[
 * z =  Du + Bv + b
 * \f]
 * \f[
 * v \geq 0, z \geq 0,  z^{T} v =0
 * \f]
 * where
 *    - \f$ u \in R^{n} \f$ \f$ v \in R^{m} \f$  and \f$z \in R^{m} \f$ are the unknowns,
 *    - \f$ a \in R^{n} \f$ and \f$ b \in R^{m} \f$
 *    - \f$ A \in R^{n \times n } \f$
 *    - \f$ B \in R^{m \times m } \f$
 *    - \f$ C \in R^{n \times m } \f$
 *    - \f$ D \in R^{m \times n } \f$
 *
 *  The MLCP main components are:
 *  - a problem (variables A,B,C,D,a,b and size of the problem), which directly corresponds to the MixedLinearComplementarityProblem structure of Numerics
 *  - the unknowns u,v and z
 *  - a NonSmoothSolver, used to define a solver and its parameters (connected to Solver_Options structure of Numerics)
 *
 *  A MLCP is connected to a simulation that handles a NonSmoothDynamicalSystem and its Topology. \n
 *  IndexSets from simulation are used to know which constraints (UnitaryRelation) are active or not. \n
 *
 * \b Construction:
 *   - XML reading (inputs = xml node with tag "OneStepNSProblem" and a Simulation*)
 *   - Constructor from data (inputs = Simulations*, id, NonSmoothSolver*) - The solver is optional.
 * Main functions:
 *
 * \b Main functions:
 *  - formalization of the problem: computes A,B,C,D,a,b using the set of "active" UnitaryRelations from the simulation and \n
 *  the unitaryBlock-matrices saved in the field unitaryBlocks.\n
 *  Functions: initialize(), computeUnitaryBlock(), preCompute()
 *  - solving of the problem: function compute(), used to call solvers from Numerics through \n
 * the mlcp_driver() interface of Numerics.
 *  - post-treatment of data: set values of y/lambda variables of the active UR (ie Interactions) using \n
 *  ouput results from the solver (u,v,z); function postCompute().
 *
 *
 */
class MLCP : public OneStepNSProblem
{
protected:

  /** n is the number of equality */
  int n;
  /** m is the size of the complementarity conditions */
  int m;
  /** contains the matrix M of a MLCP system */
  OSNSMatrix *M;
  /** contains the vector q of a MLCP system */
  SiconosVector *q;
  /** contains the vector z of a MLCP system */
  SiconosVector *z;
  /** contains the vector w of a MLCP system */
  SiconosVector *w;
  /** The MLCP instance */
  MixedLinearComplementarity_Problem numerics_problem;


  /** Storage type for M - 0: SiconosMatrix (dense), 1: Sparse Storage (embedded into OSNSMatrix) */
  int MStorageType;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  bool isWAllocatedIn;
  bool isZAllocatedIn;
  bool isMAllocatedIn;
  bool isQAllocatedIn;


  /** default constructor (private)
   */
  MLCP();

public:

  /** xml constructor
  *  \param OneStepNSProblemXML* : the XML linked-object
  *  \param Simulation *: the simulation that owns the problem
  */
  MLCP(OneStepNSProblemXML*, Simulation*);

  /** constructor from data
  *  \param Simulation *: the simulation that owns this problem
  *  \param Solver* pointer to object that contains solver algorithm and formulation \n
  *  (optional, default = NULL => read .opt file in Numerics)
  *  \param String: id of the problem (default = "unamed")
  */
  MLCP(Simulation*, NonSmoothSolver* = NULL, const std::string& = "unamed_mlcp");

  /** destructor
  */
  ~MLCP();

  // --- n ---
  /** get the value of n,
  *  \return int
  */
  inline int getn() const
  {
    return n;
  }

  // --- numerics MLCP ---
  /** get the pointer on the Numerics MLCP,
  *  \return MixedLinearComplementarity_Problem*
  */
  inline MixedLinearComplementarity_Problem* getNumericsMLCP()
  {
    return &numerics_problem;
  }


  // --- W ---
  /** get the value of w, the initial state of the DynamicalSystem
  *  \return SimpleVector
  *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
  */
  inline const SimpleVector getW() const
  {
    return *w;
  }

  /** get w, the initial state of the DynamicalSystem
  *  \return pointer on a SimpleVector
  */
  inline SiconosVector* getWPtr() const
  {
    return w;
  }

  /** set the value of w to newValue
  *  \param SiconosVector newValue
  */
  void setW(const SiconosVector&);

  /** set w to pointer newPtr
  *  \param SiconosVector * newPtr
  */
  void setWPtr(SiconosVector*);

  // --- Z ---
  /** get the value of z, the initial state of the DynamicalSystem
  *  \return SimpleVector
  *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
  */
  inline const SimpleVector getZ() const
  {
    return *z;
  }

  /** get z, the initial state of the DynamicalSystem
  *  \return pointer on a SiconosVector
  */
  inline SiconosVector* getZPtr() const
  {
    return z;
  }

  /** set the value of z to newValue
  *  \param SiconosVector newValue
  */
  void setZ(const SiconosVector&);

  /** set z to pointer newPtr
  *  \param SiconosVector * newPtr
  */
  void setZPtr(SiconosVector*) ;

  // --- M ---

  /** get M
  *  \return pointer on a SiconosMatrix
  */
  inline OSNSMatrix* getMPtr() const
  {
    return M;
  }

  /** set the value of M to newValue
  *  \param newValue SiconosMatrix
  */
  void setM(const OSNSMatrix&);

  /** set M to pointer newPtr
   *  \param newPtr SiconosMatrix*
   */
  void setMPtr(OSNSMatrix *);

  // --- Q ---
  /** get the value of q, the initial state of the DynamicalSystem
  *  \return SimpleVector
  *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
  */
  inline const SimpleVector getQ() const
  {
    return *q;
  }

  /** get q, the initial state of the DynamicalSystem
  *  \return pointer on a SiconosVector
  */
  inline SiconosVector* getQPtr() const
  {
    return q;
  }

  /** set the value of q to newValue
  *  \param SiconosVector newValue
  */
  void setQ(const SiconosVector&);

  /** set q to pointer newPtr
  *  \param SiconosVector * newPtr
  */
  void setQPtr(SiconosVector*);

  /** get the type of storage for M */
  inline const int getMStorageType() const
  {
    return MStorageType;
  };

  /** set which type of storage will be used for M - Note that this function does not
      allocate any memory for M, it just sets an indicator for future use */
  inline void setMStorageType(int i)
  {
    MStorageType = i;
  };

  /** To initialize the MLCP problem(computes topology ...)
   */
  virtual void initialize();
  virtual void updateM();

  void reset();

  /** computes extra diagonal unitaryBlock-matrix that corresponds to UR1 and UR2
  *  Move this to Unitary Relation class?
  *  \param a pointer to UnitaryRelation
  *  \param a pointer to UnitaryRelation
  */
  void computeUnitaryBlock(UnitaryRelation*, UnitaryRelation*);

  /** compute vector q
  *  \param double : current time
  */
  void computeQ(double);

  /** pre-treatment for MLCP
  *  \param double : current time
  *  \return void
  */
  void preCompute(double);

  /** Compute the unknown z and w and update the Interaction (y and lambda )
  *  \param double : current time
  *  \return int, information about the solver convergence.
  */
  int compute(double);

  /** post-treatment for MLCP
  */
  void postCompute() ;

  /** print the data to the screen
  */
  void display() const;

  /** copy the data of the OneStepNSProblem to the XML tree
  *  \exception RuntimeException
  */
  void saveNSProblemToXML();

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
  *  \param OneStepNSProblem* : the one step problem which must be converted
  * \return a pointer on the problem if it is of the right type, NULL otherwise
  */
  static MLCP* convert(OneStepNSProblem* osnsp);

};

#endif // MLCP_H
