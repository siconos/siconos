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
/*! \file LCP.h
\brief Linear Complementarity Problem formulation and solving
*/

#ifndef LCP_H
#define LCP_H

#include "OneStepNSProblem.h"
#include "SimpleVector.h"
#include "SimpleMatrix.h"
#include "SparseBlockMatrix.h"
#include <sys/time.h>

#include "SiconosPointers.h"

class OneStepNSProblem;

/** Formalization and Resolution of a Linear Complementarity Problem (LCP)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 * \section LCPintro Aim of the LCP class
 *
 * This class is devoted to the formalization and the resolution of the
 * Linear Complementarity Problem (LCP) defined by :
 *  * \f[
 * w =  q + M z
 * \f]
 * \f[
 * w \geq 0, z \geq 0,  z^{T} w =0
 * \f]
 * where
 *    - \f$ w \in R^{n} \f$  and \f$z \in R^{n} \f$ are the unknowns,
 *    - \f$ M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
 *
 *  The LCP main components are:
 *  - a problem (variables M,q and size of the problem), which directly corresponds to the LinearComplementarityProblem structure of Numerics
 *  - the unknowns z and w
 *  - a NonSmoothSolver, used to define a solver and its parameters (connected to Solver_Options structure of Numerics)
 *
 *  A LCP is connected to a simulation that handles a NonSmoothDynamicalSystem and its Topology. \n
 *  IndexSets from simulation are used to know which constraints (UnitaryRelation) are active or not. \n
 *
 * \b Construction:
 *   - XML reading (inputs = xml node with tag "OneStepNSProblem" and a Simulation*)
 *   - Constructor from data (inputs = Simulations*, id, NonSmoothSolver*) - The solver is optional.
 * Main functions:
 *
 * \b Main functions:
 *  - formalization of the problem: computes M,q using the set of "active" UnitaryRelations from the simulation and \n
 *  the block-matrices saved in the field blocks.\n
 *  Functions: initialize(), computeBlock(), preCompute()
 *  - solving of the FrictionContact problem: function compute(), used to call solvers from Numerics through \n
 * the lcp_driver() interface of Numerics.
 *  - post-treatment of data: set values of y/lambda variables of the active UR (ie Interactions) using \n
 *  ouput results from the solver (z,w); function postCompute().
 *
 *
 * \todo : add "recover" function to start from old values of z and w.
 */
class LCP : public OneStepNSProblem
{
private:

  /** contains the vector w of a LCP system */
  SiconosVectorSPtr w;

  /** contains the vector z of a LCP system */
  SiconosVectorSPtr z;

  /** contains the matrix M of a LCP system */
  OSNSMatrixSPtr M;

  /** contains the vector q of a LCP system */
  SiconosVectorSPtr q;

#ifndef WithSmartPtr
  /** Flags to check whether pointers were allocated in class constructors or not */
  bool isWAllocatedIn;
  bool isZAllocatedIn;
  bool isMAllocatedIn;
  bool isQAllocatedIn;
#endif

  /** Storage type for M - 0: SiconosMatrix (dense), 1: Sparse Storage (embedded into OSNSMatrix) */
  int MStorageType;

  /** default constructor (private)
   */
  LCP();

public:

  /** xml constructor
  *  \param OneStepNSProblemXML* : the XML linked-object
  *  \param Simulation *: the simulation that owns the problem
  */
  LCP(OneStepNSProblemXML*, Simulation*);

  /** constructor from data
  *  \param Simulation *: the simulation that owns this problem
  *  \param Solver* pointer to object that contains solver algorithm and formulation \n
  *  (optional, default = NULL => read .opt file in Numerics)
  *  \param String: id of the problem (default = "unamed")
  */
  LCP(Simulation*, NonSmoothSolver* = NULL, const std::string& = "unamed_lcp");

  /** destructor
  */
  ~LCP();

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
  inline SiconosVectorSPtr getWPtr() const
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
  void setWPtr(SiconosVectorSPtr);

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
  inline SiconosVectorSPtr getZPtr() const
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
  void setZPtr(SiconosVectorSPtr) ;

  // --- M ---

  /** get M
  *  \return pointer on a OSNSMatrix
  */
  inline OSNSMatrixSPtr getMPtr() const
  {
    return M;
  }

  /** set the value of M to newValue
  *  \param newValue OSNSMatrix
  */
  void setM(const OSNSMatrix&);

  /** set M to pointer newPtr
   *  \param newPtr OSNSMatrix*
   */
  void setMPtr(OSNSMatrixSPtr);

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
  inline SiconosVectorSPtr getQPtr() const
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
  void setQPtr(SiconosVectorSPtr);

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

  /** To initialize the LCP problem(computes topology ...)
   */
  void initialize();

  /** computes extra diagonal block-matrix that corresponds to UR1 and UR2
  *  Move this to Unitary Relation class?
  *  \param a pointer to UnitaryRelation
  *  \param a pointer to UnitaryRelation
  */
  void computeBlock(UnitaryRelation*, UnitaryRelation*);

  /** compute vector q
  *  \param double : current time
  */
  void computeQ(double);

  /** pre-treatment for LCP
  *  \param double : current time
  *  \return void
  */
  void preCompute(double);

  /** Compute the unknown z and w and update the Interaction (y and lambda )
  *  \param double : current time
  *  \return int, information about the solver convergence.
  */
  int compute(double);

  /** post-treatment for LCP
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
  static LCP* convert(OneStepNSProblem* osnsp);

};

#endif // LCP_H
