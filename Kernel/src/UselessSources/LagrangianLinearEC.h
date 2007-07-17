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

#ifndef LAGRANGIANLINEAREC_H
#define LAGRANGIANLINEAREC_H

#include "LagrangianEC.h"
#include "LagrangianLinearECXML.h"
#include "SimpleVector.h"
#include "BlockVector.h"

//! Lagrangian Linear EqualityConstraint
/**  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) Apr 27, 2004
 *
 *
 */
class LagrangianLinearEC : public LagrangianEC
{
public:

  /** Default constructor
   */
  LagrangianLinearEC();

  /** constructor with XML object of the parent class EqualityConstraint
   *  \param EqualityConstraintXML* : the XML object corresponding
   */
  LagrangianLinearEC(EqualityConstraintXML*);

  /** constructor with in parameters, the data needed to build this Relation
   *  \param a SiconosMatrix to set h
   *  \param e SiconosVector to set b
   */
  LagrangianLinearEC(const SiconosMatrix&, const SimpleVector&);

  /** Destructor */
  ~LagrangianLinearEC();

  /** getter of the SimpleMatrix h
   *  \return a pointer on the SimpleMatrix h
   */
  inline SimpleMatrix getH(void) const
  {
    return *h;
  } ;

  /** getter of the SiconosVector b
   *  \return SimpleVector : value of b
   */
  inline /*SiconosVector*/SimpleVector getB(void) const
  {
    return *b;
  };

  /** getter of the SiconosMatrix* h
   *  \return a pointer on the SiconosMatrix* h
   */
  SiconosMatrix* getHPtr(void);

  /** getter of the SiconosVector* b
   *  \return a pointer on the SiconosVector b
   */
  SiconosVector* getBPtr(void);

  /** setter on the SiconosMatrix h
   *  \param a SiconosMatrix to set h
   */
  inline void setH(const SiconosMatrix &newH)
  {
    *h = newH;
  };

  /** set the vector b
   *  \param SimpleVector& : new value of b
   */
  inline void setB(/*SiconosVector*/SimpleVector& newB)
  {
    *b = newB;
  };

  /** default function to compute y for the free state
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeFreeOutput(double time);

  /** default function to compute y
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeOutput(double time);

  /** default function to compute lambda
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeInput(double time);

  /** copy the data of the Relation to the XML tree
   *  \exception RuntimeException
   */
  void saveEqualityConstraintToXML();

  /** print the data to the screen
   */
  void display() const;

  /** allows to create the Relation with an xml file, or the needed data
   *  \param LagrangianLinearECXML * : the XML object for this EqualityConstraint
   *  \param SiconosMatrix* : the matrix H of this EqualityConstraint
   *  \param SiconosVector* : the vector h of this EqualityConstraint
   *  \exception RuntimeException
   */
  void createEqualityConstraint(EqualityConstraintXML * ecXML,
                                SiconosMatrix* H = NULL, SiconosVector* b = NULL);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the EqualityConstraint which must be converted
   * \return a pointer on the EqualityConstraint if it is of the right type, NULL otherwise
   */
  static LagrangianLinearEC* convert(EqualityConstraint *ec);

protected:
  /** uses the EqualityConstraintXML of the LagrangianLinearEC to fill the fields of this EqualityConstraint
   *  \exception RuntimeException
   */
  void fillEqualityConstraintWithEqualityConstraintXML();


private:
  /** a specific matrix to the LagrangianLinearEC */
  SiconosMatrix *h;

  /** a specific vector to the LagrangianLinearEC */
  /*SiconosVector*/
  SimpleVector *b;
};

#endif // LAGRANGIANLINEAREC_H
