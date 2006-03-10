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
#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H

#include "DynamicalSystem.h"
#include "SiconosMatrix.h"
#include "SiconosVector.h"
#include "BoundaryConditionXML.h"
#include "RuntimeException.h"
#include "SiconosConst.h"

#include "check.h"
#include <string>
#include <vector>

const std::string LINEARBC = "LinearBC";
const std::string NLINEARBC = "NonLinearBC";
const std::string PERIODICBC = "PeriodicBC";

class DynamicalSystem;

/** \class BoundaryCondition
 *  \brief represents the boundary conditions for a NSDS BVP
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.3.
 *  \date (Creation) May 24, 2004
 *
 *
 */
class BoundaryCondition
{
public:

  /** \fn BoundaryCondition();
   *  \brief Basic constructor
   */
  BoundaryCondition();

  /** \fn BoundaryCondition(BoundaryConditionXML*);
   *  \brief constructor with XML object
   *  \param The object XML which contains data of the boundary condition.
   */
  BoundaryCondition(BoundaryConditionXML*);

  virtual ~BoundaryCondition();

  /** \fn inline string getType()
   *  \brief allows to get the type of the BoundaryCondition
   *  \return string : the type of the BoundaryCondition
   */
  inline std::string  getType() const
  {
    return this->boundaryType;
  }

  /** \fn inline void setBoundaryConditionXML( BoundaryConditionXML* bcxml )
   *  \brief allows to set the BoundaryConditionXML
   *  \param BoundaryConditionXML* : the BoundaryConditionXML of the BoundaryCondition
   */
  inline void setBoundaryConditionXML(BoundaryConditionXML* bcxml)
  {
    this->bcXML = bcxml;
  }

protected:
  /** \fn void fillBCWithBCXML()
   *  \brief uses the BoundaryConditionXML of the BoundaryCondition to fill the fields of this BoundaryCondition
   */
  virtual void fillBCWithBCXML();

  /** type of condition : linear, non linear, etc. */
  std::string  boundaryType;
  BoundaryConditionXML* bcXML;

};

#endif // BOUNDARYCONDITION_H
