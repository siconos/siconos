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
/*! \file BoundaryConditionXML.h

*/

#ifndef BOUNDARYCONDITIONXML_H
#define BOUNDARYCONDITIONXML_H

#include "SiconosDOMTreeTools.h"

//! describes the boundary conditions for a BVP NSDS
/**  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date May 24, 2004
 *
 */
class BoundaryConditionXML
{
public:

  BoundaryConditionXML();
  BoundaryConditionXML(xmlNode *);
  ~BoundaryConditionXML();

  /** Return the type of the BoundaryConditionXML
   *  \return The string type of the BoundaryConditionXML
   */
  inline std::string getType()
  {
    std::string type((char*)rootBCNode->name);
    return type;
    //return SiconosDOMTreeTools::getStringAttributeValue(this->rootBCNode, BC_TYPE);
  }

protected:
  /* node of the DOM tree containing the attributes of the BouindaryCondition
   * The BoundaryCondition is formed like that : <BoundaryCondition> <[type]> <[type]/> <BoundaryCondition/>
   * the node here corrrespond to the tag <[type]>
   */
  xmlNode *rootBCNode;
};

#endif // BOUNDARYCONDITIONXML_H
