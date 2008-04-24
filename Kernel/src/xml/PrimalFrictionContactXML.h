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
/*! \file PrimalFrictionContactXML.h
  \brief XML management for PrimalFrictionContact problems
 */

#ifndef __PrimalFrictionContactXML__
#define __PrimalFrictionContactXML__

#include "OneStepNSProblemXML.h"
#include "SimpleMatrix.h"
#include "SimpleVector.h"

class OneStepNSProblem;

/** XML management for PrimalFrictionContact
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 05/18/2004
 *
 *  XML tag for PrimalFrictionContact (example for 3D case)

 \code
 <OneStepNSProblems_List>
 <PrimalFrictionContact Type="3D" StorageType="0">
 <NonSmoothSolver name="Lemke">
 <iparam  vectorSize='3'>1 3 4</iparam>
 <dparam  vectorSize='3'>1.3 1e-17 1e-12</dparam>
 </NonSmoothSolver>
 </PrimalFrictionContact>
 </OneStepNSProblems_List>

 \endcode

 attributes type is required. \n
 storageType is optional.

 Compared to base-class OneStepNSProblem, the present one just add the attribute "type".

*/
class PrimalFrictionContactXML : public OneStepNSProblemXML
{
public:
  PrimalFrictionContactXML() : OneStepNSProblemXML() {}

  /** Build a PrimalFrictionContactXML object from a DOM tree describing a PrimalFrictionContact
   *   \param PrimalFrictionContactNode : the PrimalFrictionContact DOM tree
   *   \exception XMLException : if a property of the PrimalFrictionContact lacks in the DOM tree
   */
  PrimalFrictionContactXML(xmlNodePtr PrimalFrictionContactNode) : OneStepNSProblemXML(PrimalFrictionContactNode) {}

  /** Destructor */
  ~PrimalFrictionContactXML() {}

  /** Checks if attribute "type" is given in PrimalFrictionContact tag
   *  \return a bool
   */
  inline bool hasProblemDim() const
  {
    return (SiconosDOMTreeTools::hasAttributeValue(rootNode, "Type"));
  }

  /** Returns the value of attribute "type" in PrimalFrictionContact tag (ie dim of the problem, 2D or 3D)
   *  \return an integer (2 or 3)
   */
  inline int getProblemDim() const
  {
    if (!hasProblemDim())
      XMLException::selfThrow("PrimalFrictionContactXML::getDimNSProblem - Attribute named type does not exists in tag PrimalFrictionContact.");
    return SiconosDOMTreeTools::getAttributeValue<int>(rootNode, "Type");
  }

};


#endif
