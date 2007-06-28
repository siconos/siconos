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

/*! \file
*/

/** XML management for OneStepIntegrator
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.0.
 *   \date 05/17/2004
 *
 */

#ifndef __OneStepIntegratorXML__
#define __OneStepIntegratorXML__

#include "SiconosDOMTreeTools.h"

// Tags
const std::string INTERACTIONS_CONCERNED = "Interactions_Concerned";

class OneStepIntegratorXML
{
protected:

  /** root node (tag OneStepIntegrator) */
  xmlNodePtr rootNode;

  /** DSConcerned node (tag DS_Concerned)
   * Child of rootIntegratorXMLNode */
  xmlNodePtr DSConcernedNode;

  /** interactionConcerned node (tag Interactions_Concerned)
   * Child of rootIntegratorXMLNode */
  xmlNodePtr interactionsConcernedNode;

public:

  OneStepIntegratorXML();

  /** Build a OneStepIntegratorXML object from the DOM tree
  *   \param OneStepIntegratorNode : the OneStepIntegrator DOM tree
  */
  OneStepIntegratorXML(xmlNodePtr OneStepIntegratorNode);

  /** destructor
  */
  virtual ~OneStepIntegratorXML() {};

  /** get the root node of the OneStepIntegratorXML
  *   \return xmlNodePtr
  */
  inline xmlNodePtr getRootNode() const
  {
    return rootNode;
  }

  /** get the DSConcerned node of the OneStepIntegratorXML
  *   \return xmlNodePtr
  */
  inline xmlNodePtr getDSConcernedNode() const
  {
    return DSConcernedNode;
  }

  /** get the interactionsConcerned node of the OneStepIntegratorXML
  *   \return xmlNodePtr
  */
  inline xmlNodePtr getInteractionsConcernedNode() const
  {
    return interactionsConcernedNode;
  }

  /** Return the type of the OneStepIntegratorXML
  *   \return a string
  */
  inline const std::string getType() const
  {
    return (char*)rootNode->name;
  }

  /** true if tag DSConcerned is present
  *  \return a bool
  */
  inline bool hasDSList()
  {
    return !(DSConcernedNode == NULL);
  }

  /** true if tag interactionsConcernedNode is present
  *  \return a bool
  */
  inline bool hasInteractionsList()
  {
    return !(interactionsConcernedNode == NULL);
  }

  /** attribute of the DSConcerned tag - True if all DS of the NSDS are concerned
  *  \return a bool
  */
  inline bool hasAllDS() const
  {
    if (SiconosDOMTreeTools::hasAttributeValue(DSConcernedNode, ALL_ATTRIBUTE))
      return SiconosDOMTreeTools::getAttributeValue<bool>(DSConcernedNode, ALL_ATTRIBUTE);
    else return false;
  }

  /** attribute of the interactionsConcerned tag - True if all the interactions of the nsds are concerned
  *  \return a bool
  */
  inline bool hasAllInteractions() const
  {
    if (SiconosDOMTreeTools::hasAttributeValue(interactionsConcernedNode, ALL_ATTRIBUTE))
      return SiconosDOMTreeTools::getAttributeValue<bool>(interactionsConcernedNode, ALL_ATTRIBUTE);
    else return false;
  }

  /** return a vector<int> of ds numbers related to the OSI
  *  \param in-out vector<int>
  */
  void getDSNumbers(std::vector<int>&);

  /** return a vector<int> of interactions numbers related to the OSI
  *  \param in-out vector<int>
  */
  void getInteractionsNumbers(std::vector<int>&);
};


#endif
