/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
/*! \file LinearBC.h

*/
#ifndef LINEARBC_H
#define LINEARBC_H

#include "BoundaryCondition.h"
#include "LinearBCXML.h"
#include "check.h"

//! Linear Boundary Conditions
/**  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.0.
 *  \date (Creation) May 6, 2004
 *
 *
 */
class LinearBC : public BoundaryCondition
{
public:

  /** Basic constructor
   */
  LinearBC();

  /** constructor with XML object
   *  \param The object XML which contains data of the boundary condition.
   */
  LinearBC(BoundaryConditionXML*);

  ~LinearBC();

  /** allow to get the SimpleMatrix omega0
   *  \return the SimpleMatrix omega0
   */
  inline SimpleMatrix getOmega0(void) const
  {
    return *omega0;
  };

  /** allow to get the SimpleMatrix omegaT
   *  \return the SimpleMatrix omegaT
   */
  inline SimpleMatrix getOmegaT(void) const
  {
    return *omegaT;
  };

  /** get vector omega
   *  \return SimpleVector : value of omega
   */
  inline SimpleVector getOmega(void) const
  {
    return *omega;
  };


  /** allow to set the SiconosMatrix omega0
   *  \param the SiconosMatrix to set omega0
   */
  inline void setOmega0(SiconosMatrix &M)
  {
    *omega0 = M;
  };

  /** allow to set the SiconosMatrix omegaT
   *  \param the SiconosMatrix to set omegaT
   */
  inline void setOmegaT(SiconosMatrix &M)
  {
    *omegaT = M;
  };

  /** set vector omega
   *  \param SimpleVector& : new value of omega
   */
  inline void setOmega(/*SiconosVector*/SimpleVector& v)
  {
    *omega = v;
  };

  /** copy the data of the BoundaryCondition to the XML tree
   *  \exception RuntimeException
   */
  void saveBCToXML();

  /** allows to create the BoundaryCondition with an xml file, or the needed data
   *  \param BoundaryConditionXML* : the XML object for this BoundaryCondition
   *  \param SiconosVector* : the omega vector of this BoundaryCondition
   *  \param SiconosVector* : the omega0 matrix of this BoundaryCondition
   *  \param SiconosMatrix* : the omegaT matrix of this BoundaryCondition
   *  \exception RuntimeException
   */
  void createBoundaryCondition(BoundaryConditionXML * bcXML,
                               SiconosVector* omega = NULL,
                               SiconosMatrix* omega0 = NULL, SiconosMatrix* omegaT = NULL); //,DynamicalSystem* ds=NULL);

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param BoundaryCondition* : the boundary condition which must be converted
   * \return a pointer on the boundary condition if it is of the right type, NULL otherwise
   */
  static LinearBC* convert(BoundaryCondition* bc);

protected:
  /** uses the BoundaryConditionXML of the BoundaryCondition to fill the fields of this BoundaryCondition
   *  \exception RuntimeException
   */
  void fillBCWithBCXML();

private:
  /** Initial matrix of boundary conditions */
  SiconosMatrix* omega0;
  /** Current matrix of boundary conditions */
  SiconosMatrix* omegaT;
  /** Vector omega of the BoundaryCondition */
  SimpleVector* omega;
};

#endif // LINEARBC_H

