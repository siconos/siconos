/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
/*! \file FirstOrderLinearTIR.hpp

 */
#ifndef FirstOrderLinearTIR_H
#define FirstOrderLinearTIR_H

#include "FirstOrderR.hpp"

/** Linear Time Invariant Relation, derived from class FirstOrderR

\author SICONOS Development Team - copyright INRIA
\version 3.0.0.
\date Apr 15, 2007

Linear Relation for First Order Dynamical Systems:

\f{eqnarray}
y &=& Cx(t) + Fz + D\lambda + e\\

R &=& B\lambda
\f}

 */
class FirstOrderLinearTIR : public FirstOrderR
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderLinearTIR);





  /** F matrix, coefficient of z */
  SP::SiconosMatrix _F;

  /** e*/
  SP::SiconosVector _e;


  /** default constructor, protected
  */
  //FirstOrderLinearTIR():FirstOrderR(RELATION::LinearTIR) {};

public:

  /** constructor with XML object of the parent class Relation
  *  \param RelationXML* : the XML object corresponding
  */
  FirstOrderLinearTIR(SP::RelationXML);

  FirstOrderLinearTIR();
  /** create the Relation from a set of data
  *  \param SP::SiconosMatrix : the matrix C
  *  \param SP::SiconosMatrix : the matrix B
  */
  FirstOrderLinearTIR(SP::SiconosMatrix, SP::SiconosMatrix);

  /** create the Relation from a set of data
  *  \param SP::SiconosMatrix : C
  *  \param SP::SiconosMatrix : D
  *  \param SP::SiconosMatrix : F
  *  \param SP::SiconosVector : e
  *  \param SP::SiconosMatrix : B
  */
  FirstOrderLinearTIR(SP::SiconosMatrix, SP::SiconosMatrix, SP::SiconosMatrix, SP::SiconosVector, SP::SiconosMatrix);

  /** create the Relation from a set of data
  *  \param SiconosMatrix : the matrix C
  *  \param SiconosMatrix : the matrix B
  */
  FirstOrderLinearTIR(const SiconosMatrix&, const SiconosMatrix&);

  /** create the Relation from a set of data
  *  \param SiconosMatrix : C
  *  \param SiconosMatrix : D
  *  \param SiconosMatrix : F
  *  \param SiconosVector : e
  *  \param SiconosMatrix : B
  */
  FirstOrderLinearTIR(const SiconosMatrix&, const SiconosMatrix&, const SiconosMatrix&, const SiconosVector&, const SiconosMatrix&);



  /** destructor
  */
  virtual ~FirstOrderLinearTIR() {};

  /** initialize the relation (check sizes, memory allocation ...)
  \param SP to Interaction: the interaction that owns this relation
  */
  void initialize(Interaction& inter);

  /** Gets the number of computed jacobians for h
  \return an unsigned int.
  inline unsigned int getNumberOfJacobiansForH() const
  {
  if(D) return 2;
  else return 1;
  }
  */


  /** Gets the number of computed jacobian for g
  \return an unsigned int.
  inline unsigned int getNumberOfJacobiansForG() const { return 1;}
  */

  /** default function to compute h
  *  \param double : current time
  */
  void computeh(const double time, Interaction& inter);

  /** default function to compute g
  *  \param double : current time
  */
  void computeg(const double time, Interaction& inter);

  /** default function to compute y
  *  \param double: not used
  *  \param unsigned int: not used
  */
  void computeOutput(const double time, Interaction& inter, unsigned int = 0);

  /** default function to compute r
  *  \param double : not used
  *  \param unsigned int: not used
  */
  void computeInput(const double time, Interaction& inter, unsigned int = 0);

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
    return _jachx;
  }

  /** set the value of C to newValue
  *  \param a plugged matrix

  void setC(const SiconosMatrix& newValue){
  setObject<SimpleMatrix,SP::SiconosMatrix,SiconosMatrix>(C,newValue);
  }
  */
  /** set C to pointer newPtr
  *  \param a SP to plugged matrix
  */
  inline void setCPtr(SP::SiconosMatrix newPtr)
  {
    _jachx = newPtr;
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

  void setD(const SiconosMatrix& newValue){
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
  */

  inline SP::SiconosMatrix F() const
  {
    return _F;
  }
  /** get F
  *  \return pointer on a plugged matrix
  */

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

  inline const SiconosVector getE() const { return *e; }
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
  setObject<SiconosVector, SP::SiconosVector,SiconosVector>(e,newValue);
  }
  */

  /** set e to pointer newPtr
  *  \param a SP to plugged vector
  */
  inline void setEPtr(SP::SiconosVector newPtr)
  {
    _e = newPtr;
  }

  // -- B --
  /** get the value of B
  *  \return plugged matrix

  inline const SimpleMatrix getB() const { return *B; }
  */
  /** get B
  *  \return pointer on a plugged matrix
  */
  inline SP::SiconosMatrix B() const
  {
    return _jacglambda;
  }

  /** set the value of B to newValue
  *  \param a plugged matrix

  void setB(const SiconosMatrix& newValue)
  {
  setObject<SimpleMatrix,SP::SiconosMatrix,SiconosMatrix>(B,newValue);
  }
  */

  /** set B to pointer newPtr
  *  \param a SP to plugged matrix
  */
  inline void setBPtr(SP::SiconosMatrix newPtr)
  {
    _jacglambda = newPtr;
  }

  /** get matrix Jach[index]
  *  \return a SimpleMatrix
  const SimpleMatrix getJachx() const;
  */

  /** get a pointer on matrix Jach[index]
  *  \return a pointer on a SiconosMatrix
  */

  /** get matrix Jacg[index]
  *  \return a SimpleMatrix
  const SimpleMatrix getJacg(unsigned int  index = 0) const;
  */

  /** get a pointer on matrix Jacg[index]
  *  \return a pointer on a SiconosMatrix
  */

  /** copy the data of the Relation to the XML tree
  */
  void saveRelationToXML() const;

  /** print the data to the screen
  */
  void display() const;

  virtual void computeJachx(const double time, Interaction& inter) {};
  virtual void computeJachlambda(const double time, Interaction& inter) {};
  virtual void computeJacglambda(const double time, Interaction& inter) {};
  /**
  * return true if the relation is linear.
  */

  virtual bool isLinear()
  {
    return true;
  }
  /** encapsulates an operation of dynamic casting. Needed by Python interface.
  *  \param Relation * : the relation which must be converted
  * \return a pointer on the relation if it is of the right type, NULL otherwise
  */
  static FirstOrderLinearTIR* convert(Relation *r);

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderLinearTIR);

#endif
