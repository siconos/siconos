/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
/*! \file LagrangianLinearTIR.h

*/
#ifndef LAGRANGIANLINEARRELATION_H
#define LAGRANGIANLINEARRELATION_H

#include "LagrangianR.hpp"


/**  Lagrangian Linear Relation.

  \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) Apr 27, 2004

  Lagrangian Relation with:

  \f[
  y= Cq + e + D\lambda[0] + Fz
  \f]
  \f[
  \dot y= C \dot q + D\lambda[1]
  \f]

  \f[
  p= C^t \lambda[1]
  \f]

   C is the only required input to built a LagrangianLinearTIR.

*/
class LagrangianLinearTIR : public LagrangianR
{

protected:


  /** C*/
  //SP::SiconosMatrix C;

  /** D matrix, coefficient of lambda in y */
  //SP::SiconosMatrix D;

  /** F matrix, coefficient of z */
  SP::SiconosMatrix _F;

  /** e*/
  SP::SiconosVector _e;

  /** Default constructor
   */
  LagrangianLinearTIR(): LagrangianR(RELATION::LinearTIR) {};

public:

  /** constructor with XML object of the parent class Relation
   *  \param RelationXML* : the XML object corresponding
   */
  LagrangianLinearTIR(SP::RelationXML);

  /** create the Relation from a set of data
   *  \param SP::SiconosMatrix : the matrix C
   */
  LagrangianLinearTIR(SP::SiconosMatrix);

  /** create the Relation from a set of data
   *  \param SP::SiconosMatrix : C
   *  \param SP::SiconosMatrix : D
   *  \param SP::SiconosMatrix : F
   *  \param SP::SiconosVector : e
   */
  LagrangianLinearTIR(SP::SiconosMatrix, SP::SiconosMatrix, SP::SiconosMatrix, SP::SiconosVector);

  /** create the Relation from a set of data
   *  \param SP::SiconosMatrix : C
   *  \param SP::SiconosVector : e
   */
  LagrangianLinearTIR(SP::SiconosMatrix, SP::SiconosVector);

  /** create the Relation from a set of data
   *  \param SiconosMatrix : the matrix C
   */
  LagrangianLinearTIR(const SiconosMatrix&);

  /** create the Relation from a set of data
   *  \param SiconosMatrix : C
   *  \param SiconosMatrix : D
   *  \param SiconosMatrix : F
   *  \param SiconosVector : e
   */
  LagrangianLinearTIR(const SiconosMatrix&, const SiconosMatrix&, const SiconosMatrix&, const SiconosVector&);

  /** create the Relation from a set of data
   *  \param SiconosMatrix : C
   *  \param SiconosVector : e
   */
  LagrangianLinearTIR(const SiconosMatrix&, const SiconosVector&);

  /** destructor
   */
  virtual ~LagrangianLinearTIR() {};

  /** initialize LagrangianLinearTIR specific operators.
   */
  void initComponents();

  /** default function to compute h
   *  \param double : current time
   */
  void computeh(double);
  void computeg(double time);


  /** default function to compute y
   *  \param double: not used
   *  \param unsigned int: not used
   */
  void computeOutput(double, unsigned int = 0);

  /** default function to compute r
   *  \param double : not used
   *  \param unsigned int: not used
   */
  void computeInput(double, unsigned int = 0);

  // GETTERS/SETTERS

  // -- C --
  /** get the value of C
   *  \return plugged matrix

  inline const SimpleMatrix getC() const { return *C; }
  */
  /** get C
   *  \return pointer on a plugged matrix
   */
  inline SP::SiconosMatrix C() const
  {
    return _jachq;
  }

  /** set the value of C to newValue
   *  \param a plugged matrix

  void setC(const SiconosMatrix& newValue)
    {
      setObject<SimpleMatrix,SP::SiconosMatrix,SiconosMatrix>(C,newValue);
    }
  */

  /** set C to pointer newPtr
   *  \param a SP to plugged matrix
   */
  inline void setCPtr(SP::SiconosMatrix newPtr)
  {
    _jachq = newPtr;
  }

  // -- D --

  /** get the value of D
   *  \return plugged matrix

  inline const SimpleMatrix getD() const { return *D; }
  */

  /** get D
   *  \return pointer on a plugged matrix
   */
  inline SP::SiconosMatrix D() const
  {
    return _jachlambda;
  }

  /** set the value of D to newValue
   *  \param a plugged matrix

  void setD(const SiconosMatrix& newValue)
    {
      setObject<SimpleMatrix,SP::SiconosMatrix,SiconosMatrix>(D,newValue);
    }
  */

  /** set D to pointer newPtr
   *  \param a SP to plugged matrix
   */
  inline void setDPtr(SP::SiconosMatrix newPtr)
  {
    _jachlambda = newPtr;
  }

  // -- F --

  /** get the value of F
   *  \return plugged matrix

  inline const SimpleMatrix getF() const { return *_F; }
  */

  /** get F
   *  \return pointer on a plugged matrix
   */
  inline SP::SiconosMatrix F() const
  {
    return _F;
  }

  /** set the value of F to newValue
   *  \param a plugged matrix

  void setF(const SiconosMatrix& newValue)
    {
      setObject<SimpleMatrix,SP::SiconosMatrix,SiconosMatrix>(F,newValue);
    }
  */

  /** set F to pointer newPtr
   *  \param a SP to plugged matrix
   */
  inline void setFPtr(SP::SiconosMatrix newPtr)
  {
    _F = newPtr;
  }

  // -- e --
  /** get the value of e
   *  \return plugged vector

  inline const SimpleVector getE() const { return *e; }
  */

  /** get e
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector e() const
  {
    return _e;
  }

  /** set the value of e to newValue
   *  \param a plugged vector

  void setE(const SiconosVector& newValue)
    {
      setObject<SimpleVector, SP::SiconosVector,SiconosVector>(e,newValue);
    }
  */

  /** set e to pointer newPtr
   *  \param a SP to plugged vector
   */
  inline void setEPtr(SP::SiconosVector newPtr)
  {
    _e = newPtr;
  }

  /** get matrix Jach[index]
   *  \return a SimpleMatrix
  const SimpleMatrix getJachx() const;
  const SimpleMatrix getJachlambda() const;
   */

  /** get a pointer on matrix Jach[index]
   *  \return a pointer on a SiconosMatrix
   */

  /** copy the data of the Relation to the XML tree
   */
  void saveRelationToXML() const;

  /** print the data to the screen
   */
  void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianLinearTIR* convert(Relation *r);

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(LagrangianLinearTIR);

#endif
