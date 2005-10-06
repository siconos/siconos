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
#ifndef LINEARBC_H
#define LINEARBC_H

#include "BoundaryCondition.h"
#include "LinearBCXML.h"
#include "check.h"

/** \class LinearBC
 *  \brief kind of BoundaryCondition
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) May 6, 2004
 *
 *
 */
class LinearBC : public BoundaryCondition
{
public:

  /** \fn LinearBC();
   *  \brief Basic constructor
   */
  LinearBC();

  /** \fn LinearBC(BoundaryConditionXML*);
   *  \brief constructor with XML object
   *  \param The object XML which contains data of the boundary condition.
   */
  LinearBC(BoundaryConditionXML*);

  ~LinearBC();

  /** \fn SiconosMatrix getOmega0(void)
   *  \brief allow to get the SiconosMatrix omega0
   *  \return the SiconosMatrix omega0
   */
  inline SiconosMatrix getOmega0(void) const
  {
    return this->omega0;
  };

  /** \fn SiconosMatrix getOmegaT(void)
   *  \brief allow to get the SiconosMatrix omegaT
   *  \return the SiconosMatrix omegaT
   */
  inline SiconosMatrix getOmegaT(void) const
  {
    return this->omegaT;
  };

  /** \fn SimpleVector getOmega(void)
   *  \brief get vector omega
   *  \return SimpleVector : value of omega
   */
  inline /*SiconosVector*/SimpleVector getOmega(void) const
  {
    return this->omega;
  };


  /** \fn void setOmega0(SiconosMatrix)
   *  \brief allow to set the SiconosMatrix omega0
   *  \param the SiconosMatrix to set omega0
   */
  inline void setOmega0(SiconosMatrix &M)
  {
    this->omega0 = M;
  };

  /** \fn void setOmegaT(SiconosMatrix)
   *  \brief allow to set the SiconosMatrix omegaT
   *  \param the SiconosMatrix to set omegaT
   */
  inline void setOmegaT(SiconosMatrix &M)
  {
    this->omegaT = M;
  };

  /** \fn void setOmega(SimpleVector&)
   *  \brief set vector omega
   *  \param SimpleVector& : new value of omega
   */
  inline void setOmega(/*SiconosVector*/SimpleVector& v)
  {
    this->omega = v;
  };


  //////////////////////

  /** \fn void saveBCToXML()
   *  \brief copy the data of the BoundaryCondition to the XML tree
   *  \exception RuntimeException
   */
  void saveBCToXML();

  /** \fn void createBoundaryCondition(BoundaryConditionXML * bcXML,
        SiconosVector* omega, SiconosMatrix* omega0, SiconosMatrix* omegaT)
   *  \brief allows to create the BoundaryCondition with an xml file, or the needed data
   *  \param BoundaryConditionXML* : the XML object for this BoundaryCondition
   *  \param SiconosVector* : the omega vector of this BoundaryCondition
   *  \param SiconosVector* : the omega0 matrix of this BoundaryCondition
   *  \param SiconosMatrix* : the omegaT matrix of this BoundaryCondition
   *  \exception RuntimeException
   */
  void createBoundaryCondition(BoundaryConditionXML * bcXML,
                               SiconosVector* omega = NULL,
                               SiconosMatrix* omega0 = NULL, SiconosMatrix* omegaT = NULL); //,DynamicalSystem* ds=NULL);

  /** \fn LinearBC* convert (BoundaryCondition* bc)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param BoundaryCondition* : the boundary condition which must be converted
   * \return a pointer on the boundary condition if it is of the right type, NULL otherwise
   */
  static LinearBC* convert(BoundaryCondition* bc);

protected:
  /** \fn void fillBCWithBCXML()
   *  \brief uses the BoundaryConditionXML of the BoundaryCondition to fill the fields of this BoundaryCondition
   *  \exception RuntimeException
   */
  void fillBCWithBCXML();


private:
  /** Initial matrix of boundary conditions */
  SiconosMatrix omega0;
  /** Current matrix of boundary conditions */
  SiconosMatrix omegaT;
  /** Vector omega of the BoundaryCondition */
  /*SiconosVector*/
  SimpleVector omega;
};

#endif // LINEARBC_H

