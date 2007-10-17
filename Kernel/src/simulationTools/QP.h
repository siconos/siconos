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
/*! \file
 */
#ifndef QP_H
#define QP_H

#include "OneStepNSProblem.h"
#include "SimpleMatrix.h"
#include "SimpleVector.h"

class SimpleMatrix;
class SimpleVector;

/** Quadratic Problem
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) Apr 26, 2004
 *
 *
 *
 */
class QP : public OneStepNSProblem
{
private:

  /** contains the Q matrix of a QP problem */
  SiconosMatrix* Q;

  /** contains the p vector of a QP problem */
  SimpleVector* p;

  //  /** contains the data of the QP, according to siconos/numerics */
  //  QPStructure QPMethod;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  bool isQAllocatedIn;
  bool isPAllocatedIn;

public:

  /** xml constructor
   *  \param OneStepNSProblemXML* : the XML linked-object
   *  \param Simulation *: the simulation that owns the problem (optional)
   */
  QP(OneStepNSProblemXML*, Simulation * = NULL);

  /** Destructor */
  ~QP();

  // --- GETTERS/SETTERS ---

  // --- Q ---
  /** get the value of Q
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getQ() const
  {
    return *Q;
  }

  /** get Q
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getQPtr() const
  {
    return Q;
  }

  /** set the value of Q to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setQ(const SiconosMatrix& newValue)
  {
    *Q = newValue;
  }

  /** set Q to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setQPtr(SiconosMatrix *newPtr)
  {
    if (isQAllocatedIn) delete Q;
    Q = newPtr;
    isQAllocatedIn = false;
  }

  // --- P ---
  /** get the value of p, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getP() const
  {
    return *p;
  }

  /** get p, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getPPtr() const
  {
    return p;
  }

  /** set the value of p to newValue
   *  \param SimpleVector newValue
   */
  inline void setP(const SimpleVector& newValue)
  {
    *p = newValue;
  }

  /** set p to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setPPtr(SimpleVector* newPtr)
  {
    if (isPAllocatedIn) delete p;
    p = newPtr;
    isPAllocatedIn = false;
  }

  // --- OTHER FUNCTIONS ---

  /** make the computation so solve the NS problem
   */
  void compute(const double &);

  /** copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveNSProblemToXML();

  /** copy the matrix Q of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveQToXML();

  /** copy the vector p of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void savePToXML();

  /** print the data to the screen
   */
  void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static QP* convert(OneStepNSProblem* osnsp);

};

#endif // QP_H
