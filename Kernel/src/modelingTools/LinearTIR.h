/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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
/*! \file LinearTIR.h

*/
#ifndef LINEARTIRELATION_H
#define LINEARTIRELATION_H

#include "Relation.h"
#include "LinearTIRXML.h"

/** Linear Time Invariant Relation, derived from class Relation
 *
 * \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date Apr 27, 2004
 *
 * This class defines and computes the Linear Time Invariant Relation defined by:
 *
 * \f[
 * y = h(x,t,\lambda,u,...) = C x + Fu + D \lambda + e \\
 *
 * R = g(\lambda,t) = B \lambda
 * \f]
 *
 * Warning: For this class, h and g plug-in are not used but after all connected to default plug-in functions.
 * To use plug-in rather than C, F ... matrices, call setComputeOutput and/or setComputeInput functions.
 * This will connect h and g to the given plug-in and when calling computeOutput/Input functions, h and g will
 * be used rather than C, F ... This means that, in this case, C, F and other matrices or vectors may be useless.
 * It's up to user to be coherent and to provide good connections between plug-in and matrices if necessary.
 */
class LinearTIR : public Relation
{

private:
  /** Relation is given by: \f$ y= C x + D \lambda + Fu + e \f$*/
  /** and \f$ r = B\lambda\f$ */

  /** \var C */
  SiconosMatrix* C;
  /** \var D*/
  SiconosMatrix* D;
  /** \var F*/
  SiconosMatrix* F;
  /** \var e*/
  SimpleVector* e;
  /** \var B*/
  SiconosMatrix* B;

  /** \var isAllocatedIn
   * Flags to know if pointers have been allocated in constructors or not*/
  /* the order is the one of members list above (C,D,F,e,B)  */
  std::vector<bool> isAllocatedIn;

  /** the XML object linked to the LinearTIR to read XML data */
  // LinearTIRXML * lTIRxml;

  /** Default (private) constructor
   */
  LinearTIR();

public:

  /** xml constructor
   *  \param LinearTIRXML* : the XML object corresponding
   */
  LinearTIR(RelationXML*);

  /** create the Relation from a set of data
   *  \param SiconosMatrix : the matrix C
   *  \param SiconosMatrix : the matrix B
   *  \exception RuntimeException
   */
  LinearTIR(const SiconosMatrix& , const SiconosMatrix&);

  /** create the Relation from a set of data
   *  \param SiconosMatrix : C
   *  \param SiconosMatrix : D
   *  \param SiconosMatrix : F
   *  \param SimpleVector  : e
   *  \param SiconosMatrix : B
   *  \exception RuntimeException
   */
  LinearTIR(const SiconosMatrix& , const SiconosMatrix& ,
            const SiconosMatrix& , const SimpleVector& ,
            const SiconosMatrix&);

  /** copy constructor
   *  \param a relation to copy
   *  warning: the interaction link is not copied, set a new one!
   */
  LinearTIR(const Relation &);

  /** destructor
   */
  ~LinearTIR();

  /** initialize the relation (check sizes, memory allocation ...)
   */
  void initialize();

  // GETTERS/SETTERS

  // -- C --

  /** get the value of C
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getC() const
  {
    return *C;
  }

  /** get C
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getCPtr() const
  {
    return C;
  }

  /** set the value of C to newValue
   *  \param SiconosMatrix newValue
   */
  void setC(const SiconosMatrix&);

  /** set C to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setCPtr(SiconosMatrix *);

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
  void setFPtr(SiconosMatrix *) ;

  // -- e --

  /** get the value of e
   *  \return SimpleVector
   */
  inline const SimpleVector getE() const
  {
    return *e;
  }

  /** get e
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getEPtr() const
  {
    return e;
  }

  /** set the value of e to newValue
   *  \param SimpleVector newValue
   */
  void setE(const SimpleVector&);

  /** set E to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setEPtr(SimpleVector *);

  // -- B --

  /** get the value of B
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getB() const
  {
    return *B;
  }

  /** get B
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getBPtr() const
  {
    return B;
  }

  /** set the value of B to newValue
   *  \param SiconosMatrix newValue
   */
  void setB(const SiconosMatrix&);

  /** set B to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setBPtr(SiconosMatrix *) ;

  /** get block of C corresponding to ds
   *  \param a pointer to dynamical system and a SiconosMatrix (in-out parameter)
   */
  void getCBlockDSPtr(DynamicalSystem*, SiconosMatrix&) const;

  /** get block of C corresponding to DS number n
   *  \param an int and a SiconosMatrix (in-out parameter)
   */
  void getCBlockDSPtr(const int , SiconosMatrix&) const;

  /** get block of B corresponding to ds
   *  \param a pointer to dynamical system and a SiconosMatrix (in-out parameter)
   */
  void getBBlockDSPtr(DynamicalSystem* , SiconosMatrix&) const;

  /** get block of B corresponding to DS number n
   *  \param an int and a SiconosMatrix (in-out parameter)
   */
  void getBBlockDSPtr(const int , SiconosMatrix&) const;

  // --- OTHER FUNCTIONS ---

  /** computes y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  void computeOutput(const double, const unsigned int = 0);

  /** computes yFree AND save it into y
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  void computeFreeOutput(const double = 0, const unsigned int = 0);

  /** default function to compute lambda
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  void computeInput(const double, const unsigned int);

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
  static LinearTIR* convert(Relation *r);
};

#endif // LINEARTIRELATION_H
