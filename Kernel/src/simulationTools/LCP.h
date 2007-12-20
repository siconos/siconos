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
/*! \file LCP.h
Linear Complementarity Problem
*/

#ifndef LCP_H
#define LCP_H

#include "OneStepNSProblem.h"
#include "SimpleVector.h"
#include "SimpleMatrix.h"
#include "SparseBlockMatrix.h"
#include <sys/time.h>

class OneStepNSProblem;

/** Formalization and Resolution of a Linear Complementarity Problem (LCP)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) Apr 26, 2004
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
 * \todo Correct the computation of M with a correct concatenation process
 * \todo : add "recover" function to start from old values of z and w.
 * \todo Change the assembleM for a Sparse Block Matrix by a Boost Sparse Matrix of Dense Matrix. Creation of the Class of Sparse Block Matrix
 *
 */
class LCP : public OneStepNSProblem
{
private:

  /** contains the vector w of a LCP system */
  SiconosVector *w;

  /** contains the vector z of a LCP system */
  SiconosVector *z;

  /** contains the matrix M of a LCP system */
  SiconosMatrix *M;
  //ublas::matrix<double, ublas::column_major, ublas::bounded_array<double, 1000> * M;

  /** contains the vector q of a LCP system */
  SiconosVector *q;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  bool isWAllocatedIn;
  bool isZAllocatedIn;
  bool isMAllocatedIn;
  bool isQAllocatedIn;

  /** Sparse-Block Boost Matrix. Each block is a SiconosMatrix**/
  SparseBlockMatrix *MSparseBlock;

  /** pointer to function, called to prepare M matrix. If solver-block, connected to collectBlocks, else to assembleM*/
  void (LCP::*prepareM)();

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
   *  \param string : id
   *  \param string: solver name (optional)
   *  \param unsigned int   : MaxIter (optional) required if a solver is given
   *  \param double : Tolerance (optional) -> for NLGS, Gcp, Latin
   *  \param string : NormType (optional) -> never used at the time
   *  \param double : SearchDirection (optional) -> for Latin
   *  \param double : Rho (optional) -> for RPGS (regularization parameter)
   */
  LCP(Simulation *, const std::string&,  const std::string& = DEFAULT_SOLVER, unsigned int = DEFAULT_ITER, double = DEFAULT_TOL,
      unsigned int = DEFAULT_VERBOSE, const std::string& = DEFAULT_NORMTYPE, double = DEFAULT_SEARCHDIR,
      double = DEFAULT_RHO);

  /** constructor from data
  *  \param Solver* : pointer to object that contains solver algorithm and formulation
  *  \param Simulation *: the simulation that owns this problem
  *  \param String: id of the problem (default = DEFAULT_OSNS_NAME)
  */
  LCP(Solver*, Simulation*, const std::string& = DEFAULT_OSNS_NAME);

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

  /** get the value of M
  *  \return SimpleMatrix
  */
  inline const SimpleMatrix getM() const
  {
    return *M;
  }

  /** get M
  *  \return pointer on a SiconosMatrix
  */
  inline SiconosMatrix* getMPtr() const
  {
    return M;
  }

  /** set the value of M to newValue
  *  \param SiconosMatrix newValue
  */
  void setM(const SiconosMatrix&);

  /** set M to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setMPtr(SiconosMatrix *);

  /** get the structure used to save M as a list of blocks
   *  \return a SparseBlockMatrix
   */
  inline SparseBlockMatrix* getMSparsePtr() const
  {
    return MSparseBlock;
  }

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

  /** set the Solver of the OneStepNSProblem
   *  \param: a pointer on Solver
   */
  void setSolverPtr(Solver*);

  /** used (required) in case of a switch from a solver to a solver block
   */
  void initSolver();

  /** To initialize the LCP problem(computes topology ...)
   */
  void initialize();

  /** computes extra diagonal block-matrix that corresponds to UR1 and UR2
  *  Move this to Unitary Relation class?
  *  \param a pointer to UnitaryRelation
  *  \param a pointer to UnitaryRelation
  */
  void computeBlock(UnitaryRelation*, UnitaryRelation*);

  /** built sparse-block structured matrix M using already computed blocks
   */
  void collectMBlocks();

  /** built matrix M using already computed blocks
  */
  void assembleM();

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

  /** copy the matrix M of the OneStepNSProblem to the XML tree
  *  \exception RuntimeException
  */
  void saveMToXML();

  /** copy the vector q of the OneStepNSProblem to the XML tree
  *  \exception RuntimeException
  */
  void saveQToXML();

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
  *  \param OneStepNSProblem* : the one step problem which must be converted
  * \return a pointer on the problem if it is of the right type, NULL otherwise
  */
  static LCP* convert(OneStepNSProblem* osnsp);

  /** display stat. info (CPU time and nb of iterations achieved)
   */
  void printStat();

};

#endif // LCP_H
