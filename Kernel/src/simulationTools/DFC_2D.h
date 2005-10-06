/* Siconos version 1.0, Copyright INRIA 2005.
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
#ifndef DFC_2D_H
#define DFC_2D_H

#include "OneStepNSProblem.h"
#include "DFC_2DXML.h"

// hazardous dependency
#include "Moreau.h"
#include "LagrangianLinearR.h"
#include "NewtonImpactLawNSL.h"

class OneStepNSProbem;

/** \class DFC_2D
 *  \brief This class is devoted to the formalization and the resolution of the
 * Contact Friction Dual problem (DFC_2D)
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 12, 2005
 */
class DFC_2D : public OneStepNSProblem
{
public:
  /** \fn DFC_2D()
   *  \brief default constructor
   */
  DFC_2D();

  /** \fn DFC_2D(OneStepNSProblemXML*, Strategy*=NULL)
   *  \brief xml constructor
   *  \param OneStepNSProblemXML* : the XML linked-object
   *  \param Strategy *: the strategy that owns the problem (optional)
   */
  DFC_2D(OneStepNSProblemXML*, Strategy* = NULL);

  ~DFC_2D();

  // GETTERS/SETTERS

  /** \fn int getNDfc_2D(void)
   *  \brief get the size nDfc_2D of the DFC_2D
   *  \return an int
   */
  inline const int getNDfc_2D(void) const
  {
    return nDfc_2D;
  }

  /** \fn void setNLcp(const int&)
   *  \brief set the size of the DFC_2D
   *  \param the size
   */
  inline void setNDfc_2D(const int& newValue)
  {
    nDfc_2D = newValue;
  }

  // --- W ---
  /** \fn  const SimpleVector getW(void) const
   *  \brief get the value of w, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getW() const
  {
    return *w;
  }

  /** \fn SimpleVector* getWPtr(void) const
   *  \brief get w, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getWPtr() const
  {
    return w;
  }

  /** \fn void setW(const SimpleVector& newValue)
   *  \brief set the value of w to newValue
   *  \param SimpleVector newValue
   */
  inline void setW(const SimpleVector& newValue)
  {
    *w = newValue;
  }

  /** \fn void setWPtr(SimpleVector* newPtr)
   *  \brief set w to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setWPtr(SimpleVector* newPtr)
  {
    if (isWAllocatedIn) delete w;
    w = newPtr;
    isWAllocatedIn = false;
  }

  // --- Z ---
  /** \fn  const SimpleVector getZ(void) const
   *  \brief get the value of z, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getZ() const
  {
    return *z;
  }

  /** \fn SimpleVector* getZPtr(void) const
   *  \brief get z, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getZPtr() const
  {
    return z;
  }

  /** \fn void setZ(const SimpleVector& newValue)
   *  \brief set the value of z to newValue
   *  \param SimpleVector newValue
   */
  inline void setZ(const SimpleVector& newValue)
  {
    *z = newValue;
  }

  /** \fn void setZPtr(SimpleVector* newPtr)
   *  \brief set z to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setZPtr(SimpleVector* newPtr)
  {
    if (isZAllocatedIn) delete z;
    z = newPtr;
    isZAllocatedIn = false;
  }

  // --- M ---

  /** \fn  const SiconosMatrix getM(void) const
   *  \brief get the value of M
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getM() const
  {
    return *M;
  }

  /** \fn SiconosMatrix* getMPtr(void) const
   *  \brief get M
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getMPtr() const
  {
    return M;
  }

  /** \fn void setM (const SiconosMatrix& newValue)
   *  \brief set the value of M to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setM(const SiconosMatrix& newValue)
  {
    *M = newValue;
  }

  /** \fn void setMPtr(SiconosMatrix* newPtr)
   *  \brief set M to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setMPtr(SiconosMatrix *newPtr)
  {
    if (isMAllocatedIn) delete M;
    M = newPtr;
    isMAllocatedIn = false;
  }

  // --- Q ---
  /** \fn  const SimpleVector getQ(void) const
   *  \brief get the value of q, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getQ() const
  {
    return *q;
  }

  /** \fn SimpleVector* getQPtr(void) const
   *  \brief get q, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getQPtr() const
  {
    return q;
  }

  /** \fn void setQ(const SimpleVector& newValue)
   *  \brief set the value of q to newValue
   *  \param SimpleVector newValue
   */
  inline void setQ(const SimpleVector& newValue)
  {
    *q = newValue;
  }

  /** \fn void setQPtr(SimpleVector* newPtr)
   *  \brief set q to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setQPtr(SimpleVector* newPtr)
  {
    if (isQAllocatedIn) delete q;
    q = newPtr;
    isQAllocatedIn = false;
  }

  /** \fn void compute(const double&)
   *  \brief Compute the unknown z and w and update the Interaction (y and lambda )
   *  \return void
   */
  void preDFC_2D(const double&);

  /** \fn void compute(void)
   *  \brief Compute the unknown z and w and update the Interaction (y and lambda )
   *  \return void
   */
  void compute(const double&);

  /** \fn void computeM (void)
   *  \brief make the computation of matrix M
   */
  void computeM(void);

  /** \fn void computeQ (void)
   *  \brief make the computation of the SiconosVector q
   *  \param double : current time
   */
  void computeQ(const double& time);

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveNSProblemToXML();

  /** \fn void saveMToXML()
   *  \brief copy the matrix M of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveMToXML();

  /** \fn void saveQToXML()
   *  \brief copy the vector q of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveQToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;

  /** \fn DFC_2D* convert (OneStepNSProblem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static DFC_2D* convert(OneStepNSProblem* osnsp);

private:

  /** Size of the DFC_2D */
  int nDfc_2D;

  /** contains the vector w of a DFC_2D system */
  SimpleVector* w;

  /** contains the vector z of a DFC_2D system */
  SimpleVector* z;

  /** contains the matrix M of a DFC_2D system */
  SiconosMatrix* M;

  /** contains the vector q of a DFC_2D system */
  SimpleVector* q;
  /** Flags to check wheter pointers were allocated in class constructors or not */
  bool isWAllocatedIn;
  bool isZAllocatedIn;
  bool isMAllocatedIn;
  bool isQAllocatedIn;
};

#endif // DFC_2D_H
