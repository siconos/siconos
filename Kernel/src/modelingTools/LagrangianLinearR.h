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
/*! \file LagrangianLinearR.h

*/
#ifndef LAGRANGIANLINEARRELATION_H
#define LAGRANGIANLINEARRELATION_H

#include "LagrangianR.h"

class SimpleVector;
class LagrangianR;
class LagrangianLinearRXML;

/**  Lagrangian Linear Relation.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) Apr 27, 2004
 *
 * Lagrangian Relation with:
 *
 * \f[
 * y= H q + b + D\lambda + Fz
 * \f]
 * \f[
 * \dot y= H \dot q + D\lambda
 * \f]
 *
 * \f[
 * p= H^t \lambda
 * \f]
 *
 *  H is the only required input to built a LagrangianLinearR.
 *
 */
class LagrangianLinearR : public LagrangianR
{

private:

  /** H matrix such as y = Hq + ...*/
  SiconosMatrix* H;

  /** b vector such y = Hq + b + ...*/
  SimpleVector* b;

  /** D matrix, coefficient of lambda in y */
  SiconosMatrix * D;

  /** F matrix, coefficient of z */
  SiconosMatrix * F;

  /** Default constructor
   */
  LagrangianLinearR();

  /** set all allocation flags to input value.
   *  \param bool
   */
  void initAllocationFlags(bool);

public:

  /** constructor with XML object of the parent class Relation
   *  \param RelationXML* : the XML object corresponding
   */
  LagrangianLinearR(RelationXML*);

  /** constructor with in parameters, the data needed to build this Relation
   *  \param a SiconosMatrix to set H
   *  \param a SimpleVector to set b
   */
  LagrangianLinearR(const SiconosMatrix&, const SimpleVector&);


  /** constructor with in parameters, the data needed to build this Relation
   *  \param a SiconosMatrix to set H
   */
  LagrangianLinearR(const SiconosMatrix&);

  /** constructor with in parameters, the data needed to build this Relation
   *  \param a SiconosMatrix to set H
   *  \param a SimpleVector to set b
   *  \param a SiconosMatrix to set D
   */
  LagrangianLinearR(const SiconosMatrix&, const SimpleVector&, const SiconosMatrix&);

  /** constructor with in parameters, the data needed to build this Relation
   *  \param a SiconosMatrix to set H
   *  \param a SimpleVector to set b
   *  \param a SiconosMatrix to set D
   *  \param a SiconosMatrix to set F
   */
  LagrangianLinearR(const SiconosMatrix&, const SimpleVector&, const SiconosMatrix&, const SiconosMatrix&);

  /** destructor
   */
  ~LagrangianLinearR();

  /** initialize LagrangianLinearR specific operators.
   */
  void initComponents();

  // --- GETTERS/SETTERS
  // -- H --

  /** get the value of H
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getH() const
  {
    return *H;
  }

  /** get H
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getHPtr() const
  {
    return H;
  }

  /** set the value of H to newValue
   *  \param SiconosMatrix newValue
   */
  void setH(const SiconosMatrix&);

  /** set H to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setHPtr(SiconosMatrix *);
  // -- b --

  /** get the value of b
   *  \return SimpleVector
   */
  inline const SimpleVector getB() const
  {
    return *b;
  }

  /** get b
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getBPtr() const
  {
    return b;
  }

  /** set the value of b to newValue
   *  \param SimpleVector newValue
   */
  void setB(const SimpleVector&);

  /** set B to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setBPtr(SimpleVector *);

  // -- D --

  /** get the value of D
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getD() const
  {
    return *D;
  }

  /** get D
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getDPtr() const
  {
    return D;
  }

  /** set the value of D to newValue
   *  \param SiconosMatrix newValue
   */
  void setD(const SiconosMatrix&);

  /** set D to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setDPtr(SiconosMatrix *);

  // -- F --

  /** get the value of F
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getF() const
  {
    return *F;
  }

  /** get F
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getFPtr() const
  {
    return F;
  }

  /** set the value of F to newValue
   *  \param SiconosMatrix newValue
   */
  void setF(const SiconosMatrix&);

  /** set F to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setFPtr(SiconosMatrix *);

  // --- OTHER FUNCTIONS ---

  /** default function to compute y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  void computeOutput(double, unsigned int = 0);

  /** default function to compute lambda
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  void computeInput(double, unsigned int);

  /** copy the data of the Relation to the XML tree
   *  \exception RuntimeException
   */
  void saveRelationToXML() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianLinearR* convert(Relation *r);

  /** main relation members display
   */
  void display() const;
};

#endif // LAGRANGIANLINEARRELATION_H
