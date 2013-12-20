/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

public:

  /** default constructor, protected
  */
  FirstOrderLinearTIR();

  /** constructor with XML object of the parent class Relation
  *  \param relxml the XML corresponding object
  */
  FirstOrderLinearTIR(SP::RelationXML relxml);

  /** create the Relation from a set of data
  *  \param C the matrix C
  *  \param B the matrix B
  */
  FirstOrderLinearTIR(SP::SiconosMatrix C, SP::SiconosMatrix B);

  /** create the Relation from a set of data
  *  \param C the C matrix
  *  \param D the D matrix
  *  \param F the F matrix
  *  \param e the e matrix
  *  \param B the B matrix
  */
  FirstOrderLinearTIR(SP::SiconosMatrix C, SP::SiconosMatrix D, SP::SiconosMatrix F, SP::SiconosVector e, SP::SiconosMatrix B);

  /** destructor
  */
  virtual ~FirstOrderLinearTIR() {};

  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction that owns this relation
  */
  void initialize(Interaction& inter);

  /** default function to compute h
  *  \param time current time
  *  \param inter the interaction that owns this relation
  */
  void computeh(double time, Interaction& inter);

  /** default function to compute g
  *  \param time current time
  *  \param inter the interaction that owns this relation
  */
  void computeg(double time, Interaction& inter);

  /** default function to compute y
  *  \param double: not used
  *  \param unsigned int: not used
  */
  void computeOutput(double time, Interaction& inter, unsigned int = 0);

  /** default function to compute r
  *  \param double : not used
  *  \param unsigned int: not used
  */
  void computeInput(double time, Interaction& inter, unsigned int = 0);

  // GETTERS/SETTERS

  /** get C
  *  \return pointer on a plugged matrix
  */
  inline SP::SiconosMatrix C() const
  {
    return _jachx;
  }

  /** set C to pointer newPtr
  *  \param a SP to plugged matrix
  */
  inline void setCPtr(SP::SiconosMatrix newPtr)
  {
    _jachx = newPtr;
  }

  /** get D
  *  \return pointer on a plugged matrix
  */
  inline SP::SiconosMatrix D() const
  {
    return _jachlambda;
  }

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

  /** set F to pointer newPtr
  *  \param a SP to plugged matrix
  */
  inline void setFPtr(SP::SiconosMatrix newPtr)
  {
    _F = newPtr;
  }

  /** get e
  *  \return pointer on a plugged vector
  */
  inline SP::SiconosVector e() const
  {
    return _e;
  }


  /** set e to pointer newPtr
  *  \param a SP to plugged vector
  */
  inline void setEPtr(SP::SiconosVector newPtr)
  {
    _e = newPtr;
  }

  /** get B
  *  \return pointer on a plugged matrix
  */
  inline SP::SiconosMatrix B() const
  {
    return _jacglambda;
  }

  /** set B to pointer newPtr
  *  \param a SP to plugged matrix
  */
  inline void setBPtr(SP::SiconosMatrix newPtr)
  {
    _jacglambda = newPtr;
  }

  /** copy the data of the Relation to the XML tree
  */
  void saveRelationToXML() const;

  /** print the data to the screen
  */
  void display() const;

  virtual void computeJachx(double time, Interaction& inter) {};
  virtual void computeJachlambda(double time, Interaction& inter) {};
  virtual void computeJacglambda(double time, Interaction& inter) {};

  /**
  * return true if the relation is linear.
  */
  virtual bool isLinear()
  {
    return true;
  }

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderLinearTIR)

#endif
