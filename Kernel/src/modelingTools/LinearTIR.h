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
#ifndef LINEARTIRELATION_H
#define LINEARTIRELATION_H

#include "Relation.h"
#include "LinearTIRXML.h"

/** \class LinearTIR
 *  \brief Linear Time Invariant Relation, derived from class Relation
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.3.
 *  \date Apr 27, 2004
 *
 * This class defines and computes the Linear Time Invariant Relation defined by:
 *
 * \f[
 * y = h(x,t,\lambda,u,...) = C x + Fu + D \lambda + e
 *
 * R = g(\lambda,t) = B \lambda +a
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
public:

  /** \fn LinearTIR(RelationXML*, Interaction* =NULL)
   *  \brief xml constructor
   *  \param LinearTIRXML* : the XML object corresponding
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   */
  LinearTIR(RelationXML*, Interaction* = NULL);

  /** \fn void LinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newB)
   *  \brief create the Relation from a set of data
   *  \param SiconosMatrix : the matrix C
   *  \param SiconosMatrix : the matrix B
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   *  \exception RuntimeException
      */
  LinearTIR(const SiconosMatrix& , const SiconosMatrix&, Interaction* = NULL);

  /** \fn void LinearTIR(const SiconosMatrix& newC, const SiconosMatrix& newD,
   *                     const SiconosMatrix& newF, const SimpleVector& newE,
   *                     const SiconosMatrix& newB, const SimpleVector& newA)
   *  \brief create the Relation from a set of data
   *  \param SiconosMatrix : C
   *  \param SiconosMatrix : D
   *  \param SiconosMatrix : F
   *  \param SimpleVectorx : e
   *  \param SiconosMatrix : B
   *  \param SimpleVector : a
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   *  \exception RuntimeException
   */
  LinearTIR(const SiconosMatrix& , const SiconosMatrix& ,
            const SiconosMatrix& , const SimpleVector& ,
            const SiconosMatrix& , const SimpleVector& , Interaction* = NULL);

  /** \fn LinearTIR(const Relation&)
   *  \brief copy constructor
   *  \param a relation to copy
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   *  warning: the interaction link is not copied, set a new one!
   */
  LinearTIR(const Relation &, Interaction* = NULL);

  /** \fn ~LinearTIR()
   *  \brief destructor
   */
  ~LinearTIR();

  // GETTERS/SETTERS

  // -- C --

  /** \fn  const SimpleMatrix getC() const
   *  \brief get the value of C
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getC() const
  {
    return *C;
  }

  /** \fn SiconosMatrix* getCPtr() const
   *  \brief get C
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getCPtr() const
  {
    return C;
  }

  /** \fn void setC (const SiconosMatrix& newValue)
   *  \brief set the value of C to newValue
   *  \param SiconosMatrix newValue
   */
  void setC(const SiconosMatrix&);

  /** \fn void setCPtr(SiconosMatrix* newPtr)
   *  \brief set C to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setCPtr(SiconosMatrix *);

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

  // -- F --

  /** \fn  const SimpleMatrix getF() const
   *  \brief get the value of F
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getF() const
  {
    return *F;
  }

  /** \fn SiconosMatrix* getFPtr() const
   *  \brief get F
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getFPtr() const
  {
    return F;
  }

  /** \fn void setF (const SiconosMatrix& newValue)
   *  \brief set the value of F to newValue
   *  \param SiconosMatrix newValue
   */
  void setF(const SiconosMatrix&);

  /** \fn void setFPtr(SiconosMatrix* newPtr)
   *  \brief set F to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setFPtr(SiconosMatrix *) ;

  // -- e --

  /** \fn  const SimpleVector getE() const
   *  \brief get the value of e
   *  \return SimpleVector
   */
  inline const SimpleVector getE() const
  {
    return *e;
  }

  /** \fn SimpleVector* getEPtr() const
   *  \brief get e
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getEPtr() const
  {
    return e;
  }

  /** \fn void setE (const SimpleVector& newValue)
   *  \brief set the value of e to newValue
   *  \param SimpleVector newValue
   */
  void setE(const SimpleVector&);

  /** \fn void setEPtr(SimpleVector* newPtr)
   *  \brief set E to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setEPtr(SimpleVector *);

  // -- B --

  /** \fn  const SimpleMatrix getB() const
   *  \brief get the value of B
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getB() const
  {
    return *B;
  }

  /** \fn SiconosMatrix* getBPtr() const
   *  \brief get B
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getBPtr() const
  {
    return B;
  }

  /** \fn void setB (const SiconosMatrix& newValue)
   *  \brief set the value of B to newValue
   *  \param SiconosMatrix newValue
   */
  void setB(const SiconosMatrix&);

  /** \fn void setBPtr(SiconosMatrix* newPtr)
   *  \brief set B to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setBPtr(SiconosMatrix *) ;

  // -- a --

  /** \fn  const SimpleVector getA() const
   *  \brief get the value of a
   *  \return SimpleVector
   */
  inline const SimpleVector getA() const
  {
    return *a;
  }

  /** \fn SimpleVector* getAPtr() const
   *  \brief get a
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getAPtr() const
  {
    return a;
  }

  /** \fn void setA (const SimpleVector& newValue)
   *  \brief set the value of a to newValue
   *  \param SimpleVector newValue
   */
  void setA(const SimpleVector&);

  /** \fn void setAPtr(SimpleVector* newPtr)
   *  \brief set a to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setAPtr(SimpleVector *);

  /** \fn void getCBlockDSPtr(DynamicalSystem* ds, SiconosMatrix&) const
   *  \brief get block of C corresponding to ds
   *  \param a pointer to dynamical system and a SiconosMatrix (in-out parameter)
   */
  void getCBlockDSPtr(DynamicalSystem*, SiconosMatrix&) const;

  /** \fn void getCBlockDSPtr(const int & n, SiconosMatrix&) const
   *  \brief get block of C corresponding to DS number n
   *  \param an int and a SiconosMatrix (in-out parameter)
   */
  void getCBlockDSPtr(const int &, SiconosMatrix&) const;

  /** \fn void getBBlockDSPtr(DynamicalSystem* ds, SiconosMatrix&) const
   *  \brief get block of B corresponding to ds
   *  \param a pointer to dynamical system and a SiconosMatrix (in-out parameter)
   */
  void getBBlockDSPtr(DynamicalSystem* , SiconosMatrix&) const;

  /** \fn void getBBlockDSPtr(const int & n , SiconosMatrix&) const
   *  \brief get block of B corresponding to DS number n
   *  \param an int and a SiconosMatrix (in-out parameter)
   */
  void getBBlockDSPtr(const int &, SiconosMatrix&) const;

  // --- OTHER FUNCTIONS ---

  /** \fn void computeOutput(const double& time)
   *  \brief computes y
   *  \param double : current time
   */
  void computeOutput(const double&);

  /** \fn void computeFreeOutput(const double& time)
   *  \brief computes yFree AND save it into y
   *  \param double : current time
   */
  void computeFreeOutput(const double& = 0);

  /** \fn void computeInput(double time);
   *  \brief default function to compute lambda
   *  \param double : current time
   *  \exception RuntimeException
   */
  void computeInput(const double&);

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the Relation to the XML tree
   */
  void saveRelationToXML() const;

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn LinearTIR* convert (Relation *r)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation * : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LinearTIR* convert(Relation *r);

private:
  /** \fn LinearTIR();
   *  \brief Default (private) constructor
   */
  LinearTIR();

  /** Relation is given by: \f$ y= C x + D \lambda + Fu + e \f$*/
  /** and \f$ r = B\lambda\f$ */
  /** */
  SiconosMatrix* C;
  /** */
  SiconosMatrix* D;
  /** */
  SiconosMatrix* F;
  /** */
  SimpleVector* e;
  /** */
  SiconosMatrix* B;
  /** */
  SimpleVector* a;
  /** Flags to know if pointers have been allocated in constructors or not*/
  /* the order is the one of members list above (C,D,F,e,B,a)  */
  std::vector<bool> isAllocatedIn;

  /** the XML object linked to the LinearTIR to read XML data */
  // LinearTIRXML * lTIRxml;
};

#endif // LINEARTIRELATION_H
