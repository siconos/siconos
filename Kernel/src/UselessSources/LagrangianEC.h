/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
/*! \file
*/

#ifndef LAGRANGIANEC_H
#define LAGRANGIANEC_H

#include "EqualityConstraint.h"
#include "LagrangianECXML.h"

//! Lagrangian EqualityConstraint
/**  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date 17/01/2005
 *
 *
 */
class LagrangianEC : public EqualityConstraint
{
public:

  /** default constructor
  */
  LagrangianEC();
  ~LagrangianEC();

  LagrangianEC(EqualityConstraintXML*);


  /** allows to create the EqualityConstraint with an xml file, or the needed data
  *  \param LagrangianECXML * : the XML object for this EqualityConstraint
  *  \exception RuntimeException
  */
  void createEqualityConstraint(EqualityConstraintXML * ecXML , int number = -1,
                                SiconosMatrix *G = NULL, std::vector<DSInputOutput*> *dsioVector = NULL);

private:

};

#endif //LAGRANGIANEC_H
