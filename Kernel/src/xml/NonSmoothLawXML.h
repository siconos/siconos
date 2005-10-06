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
/** \class NSLawXML
*   \brief This class manages NSLaw data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/13/2004
*
*
*
* NSLawXML allows to manage data of a NSLaw DOM tree.
*/


#ifndef __NSLawXML__
#define __NSLawXML__

#include "SiconosDOMTreeTools.h"
#include "NonSmoothLaw.h"

class NonSmoothLaw;

//const string NSLAW_TYPE = "type";
class NonSmoothLawXML
{
public:
  NonSmoothLawXML();
  NonSmoothLawXML(xmlNode*);
  virtual ~NonSmoothLawXML();

  /** \fn int getType()
  *   \brief Return the type of the NSLawXML
  *   \return The string type of the NSLawXML
  */
  inline std::string getType()
  {
    //return SiconosDOMTreeTools::getStringAttributeValue(this->rootNSLawXMLNode, NSLAW_TYPE);
    std::string type((char*) rootNSLawXMLNode->name);
    return type;
  }

  /** \fn xmlNode* getNode()
  *   \brief Return the node of the NonSmoothLawXML
  *   \return xmlNode* : the node of the NonSmoothLawXML in the DOM tree
  */
  inline xmlNode* getNode() const
  {
    return rootNSLawXMLNode;
  }

  /** \fn void updateNonSmoothLawXML( xmlNode* node, NonSmoothLaw* nsl )
  *   \brief makes the operations to create the NonSmoothLaw of the Interaction
  *   \param xmlNode* : the root node of the NonSmoothLawXML
  *   \param Relation* : the NonSmoothLaw of this NonSmoothLawXML
  */
  void updateNonSmoothLawXML(xmlNode* node, NonSmoothLaw* nsl);

protected:
  xmlNode * rootNSLawXMLNode;
};


#endif
