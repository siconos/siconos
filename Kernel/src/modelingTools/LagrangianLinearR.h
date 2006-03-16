/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
#ifndef LAGRANGIANLINEARRELATION_H
#define LAGRANGIANLINEARRELATION_H

#include "LagrangianR.h"
#include "LagrangianLinearRXML.h"

#include "SimpleVector.h"
#include "BlockVector.h"

/** \class LagrangianLinearR
 *  \brief Lagrangian Linear Relation, derived from class LagrangianR
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.3.
 *  \date (Creation) Apr 27, 2004
 *
 * This class defines and computes the Lagrangian Linear Relations defined by,
 * for the input \f$ y \f$,
 * \f[
 * y= H q + b + D\lambda
 * \f]
 * and for the output \f$ p\f$ defined by
 * \f[
 * p= H^t \lambda
 * \f]
 *  H is the only required input to built a LagrangianLinearR.
 *
 */

class LagrangianLinearR : public LagrangianR
{

private:

  /** \fn LagrangianLinearR()
   *  \brief Default constructor
   */
  LagrangianLinearR();

  /** H matrix such as y = Hq + ...*/
  SiconosMatrix* H;

  /** b vector such y = Hq + b + ...*/
  SimpleVector* b;

  /** D matrix, coefficient of lambda in y */
  SiconosMatrix * D;

  /* Flags to control where pointers have been allocated (in or out constructors)*/
  bool isHAllocatedIn;
  bool isBAllocatedIn;
  bool isDAllocatedIn;

public:

  /** \fn LagrangianLinearR(RelationXML*, Interaction* =NULL)
   *  \brief constructor with XML object of the parent class Relation
   *  \param RelationXML* : the XML object corresponding
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   */
  LagrangianLinearR(RelationXML*, Interaction* = NULL);

  /** \fn LagrangianLinearR(const SiconosMatrix& H, const SimpleVector& q, Interaction* = NULL);
   *  \brief constructor with in parameters, the data needed to build this Relation
   *  \param a SiconosMatrix to set H
   *  \param a SimpleVector to set b
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   */
  LagrangianLinearR(const SiconosMatrix&, const SimpleVector&, Interaction* = NULL);


  /** \fn LagrangianLinearR(const SiconosMatrix& H, Interaction* = NULL);
   *  \brief constructor with in parameters, the data needed to build this Relation
   *  \param a SiconosMatrix to set H
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   */
  LagrangianLinearR(const SiconosMatrix&, Interaction* = NULL);

  /** \fn LagrangianLinearR(const SiconosMatrix& H, const SimpleVector& q, const SiconosMatrix& D, Interaction* = NULL);
   *  \brief constructor with in parameters, the data needed to build this Relation
   *  \param a SiconosMatrix to set H
   *  \param a SimpleVector to set b
   *  \param a SiconosMatrix to set D
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   */
  LagrangianLinearR(const SiconosMatrix&, const SimpleVector&, const SiconosMatrix&, Interaction* = NULL);

  /** \fn LagrangianLinearR(const Relation&)
   *  \brief copy constructor
   *  \param a relation to copy
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   *  warning: the interaction link is not copied, set a new one!
   */
  LagrangianLinearR(const Relation &, Interaction* = NULL);

  ~LagrangianLinearR();

  // --- GETTERS/SETTERS
  // -- H --

  /** \fn  const SimpleMatrix getH() const
   *  \brief get the value of H
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getH() const
  {
    return *H;
  }

  /** \fn SiconosMatrix* getHPtr() const
   *  \brief get H
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getHPtr() const
  {
    return H;
  }

  /** \fn void setH (const SiconosMatrix& newValue)
   *  \brief set the value of H to newValue
   *  \param SiconosMatrix newValue
   */
  void setH(const SiconosMatrix&);

  /** \fn void setHPtr(SiconosMatrix* newPtr)
   *  \brief set H to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setHPtr(SiconosMatrix *);
  // -- b --

  /** \fn  const SimpleVector getB() const
   *  \brief get the value of b
   *  \return SimpleVector
   */
  inline const SimpleVector getB() const
  {
    return *b;
  }

  /** \fn SimpleVector* getBPtr() const
   *  \brief get b
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getBPtr() const
  {
    return b;
  }

  /** \fn void setB (const SimpleVector& newValue)
   *  \brief set the value of b to newValue
   *  \param SimpleVector newValue
   */
  void setB(const SimpleVector&);

  /** \fn void setBPtr(SimpleVector* newPtr)
   *  \brief set B to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setBPtr(SimpleVector *);

  // -- D --

  /** \fn  const SimpleMatrix getD() const
   *  \brief get the value of D
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getD() const
  {
    return *D;
  }

  /** \fn SiconosMatrix* getDPtr() const
   *  \brief get D
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getDPtr() const
  {
    return D;
  }

  /** \fn void setD (const SiconosMatrix& newValue)
   *  \brief set the value of D to newValue
   *  \param SiconosMatrix newValue
   */
  void setD(const SiconosMatrix&);

  /** \fn void setDPtr(SiconosMatrix* newPtr)
   *  \brief set D to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setDPtr(SiconosMatrix *);


  /** \fn void getHBlockDS(const int&,SiconosMatrix&) const
   *  \brief get in Matrix H the block corresponding to DS number int
   *  \param int, the ds number
   *  \param SiconosMatrix (in-out parameter): the resulting block matrix
   */
  void getHBlockDS(const int&, SiconosMatrix&) const;


  /** void getHBlockDS(DynamicalSystem * ds, SiconosMatrix&) const
   *  \brief get in Matrix H the block corresponding to ds
   *  \param a pointer to a dynamical system
   *  \param SiconosMatrix (in-out parameter): the resulting block matrix
   */
  void getHBlockDS(DynamicalSystem *, SiconosMatrix&) const;

  // --- OTHER FUNCTIONS ---

  /** \fn void computeFreeOutput(double time);
   *  \brief default function to compute y for the free state
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeFreeOutput(const double& time);

  /** \fn void computeOutput(double time);
   *  \brief default function to compute y
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeOutput(const double& time);

  /** \fn void computeInput(double time);
   *  \brief default function to compute lambda
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeInput(const double& time);

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the Relation to the XML tree
   *  \exception RuntimeException
   */
  void saveRelationToXML() const;

  /** \fn LagrangianLinearR* convert (Relation *r)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianLinearR* convert(Relation *r);

  /** \fn  void display() const
   * \brief main relation members display
   */
  void display() const;
};

#endif // LAGRANGIANLINEARRELATION_H
