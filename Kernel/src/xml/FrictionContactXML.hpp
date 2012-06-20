/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
/*! \file FrictionContactXML.hpp
  \brief XML management for FrictionContact problems
 */

#ifndef __FrictionContactXML__
#define __FrictionContactXML__

#include "OneStepNSProblemXML.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosVector.hpp"

class OneStepNSProblem;

/** XML management for FrictionContact
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 05/18/2004
 *
 *  XML tag for FrictionContact (example for 3D case)

 \code
 <OneStepNSProblems_List>
 <FrictionContact Type="3D" StorageType="0">
 <NonSmoothSolver name="Lemke">
 <iparam  vectorSize='3'>1 3 4</iparam>
 <dparam  vectorSize='3'>1.3 1e-17 1e-12</dparam>
 </NonSmoothSolver>
 </FrictionContact>
 </OneStepNSProblems_List>

 \endcode

 attributes type is required. \n
 storageType is optional.

 Compared to base-class OneStepNSProblem, the present one just add the attribute "type".

*/
class FrictionContactXML : public OneStepNSProblemXML
{
public:
  FrictionContactXML() : OneStepNSProblemXML() {}

  /** Build a FrictionContactXML object from a DOM tree describing a FrictionContact
   *   \param FrictionContactNode : the FrictionContact DOM tree
   *   \exception XMLException : if a property of the FrictionContact lacks in the DOM tree
   */
  FrictionContactXML(xmlNodePtr FrictionContactNode) : OneStepNSProblemXML(FrictionContactNode) {}

  /** Destructor */
  ~FrictionContactXML() {}

  /** Checks if attribute "type" is given in FrictionContact tag
   *  \return a bool
   */
  inline bool hasProblemDim() const
  {
    return (SiconosDOMTreeTools::hasAttributeValue(rootNode, "Type"));
  }

  /** Returns the value of attribute "type" in FrictionContact tag (ie dim of the problem, 2D or 3D)
   *  \return an integer (2 or 3)
   */
  inline int getProblemDim() const
  {
    if (!hasProblemDim())
      XMLException::selfThrow("FrictionContactXML::getDimNSProblem - Attribute named type does not exists in tag FrictionContact.");
    return SiconosDOMTreeTools::getAttributeValue<int>(rootNode, "Type");
  }

};

DEFINE_SPTR(FrictionContactXML);
#endif
