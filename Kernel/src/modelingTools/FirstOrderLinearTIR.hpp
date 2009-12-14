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
/*! \file FirstOrderLinearTIR.h

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




  /** F matrix, coefficient of z */
  SP::SiconosMatrix F;

  /** e*/
  SP::SiconosVector e;


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
  void initialize(SP::Interaction);

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
  void computeh(double);

  /** default function to compute g
   *  \param double : current time
   */
  void computeg(double);

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
  inline SP::SiconosMatrix getCPtr() const
  {
    return JacXH;
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
    JacXH = newPtr;
  }

  // -- D --

  /** get the value of D
   *  \return plugged matrix

  inline const SimpleMatrix getD() const { return *D; }
  */
  /** get D
   *  \return pointer on a plugged matrix
   */
  inline SP::SiconosMatrix getDPtr() const
  {
    return JacLH;
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
    JacLH = newPtr;
  }

  // -- F --

  /** get the value of F
   *  \return plugged matrix

  inline const SimpleMatrix getF() const { return *F; }
  */
  /** get F
   *  \return pointer on a plugged matrix
   */
  inline SP::SiconosMatrix getFPtr() const
  {
    return F;
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
    F = newPtr;
  }

  // -- e --
  /** get the value of e
   *  \return plugged vector

  inline const SimpleVector getE() const { return *e; }
  */
  /** get e
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector getEPtr() const
  {
    return e;
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
    e = newPtr;
  }

  // -- B --
  /** get the value of B
   *  \return plugged matrix

  inline const SimpleMatrix getB() const { return *B; }
  */
  /** get B
   *  \return pointer on a plugged matrix
   */
  inline SP::SiconosMatrix getBPtr() const
  {
    return JacLG;
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
    JacLG = newPtr;
  }

  /** get matrix JacH[index]
   *  \return a SimpleMatrix
  const SimpleMatrix getJacXH() const;
   */

  /** get a pointer on matrix JacH[index]
   *  \return a pointer on a SiconosMatrix
   */

  /** get matrix JacG[index]
   *  \return a SimpleMatrix
  const SimpleMatrix getJacG(unsigned int  index = 0) const;
   */

  /** get a pointer on matrix JacG[index]
   *  \return a pointer on a SiconosMatrix
   */

  /** copy the data of the Relation to the XML tree
   */
  void saveRelationToXML() const;

  /** print the data to the screen
   */
  void display() const;

  virtual void computeJacXH(double) {};
  virtual void computeJacLH(double) {};
  virtual void computeJacLG(double) {};

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static FirstOrderLinearTIR* convert(Relation *r);
};

TYPEDEF_SPTR(FirstOrderLinearTIR);

#endif
