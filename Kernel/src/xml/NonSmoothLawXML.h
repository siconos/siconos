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
/** \class NSLawXML
 *   \brief This class manages NSLaw data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.1.3.
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
protected:
  xmlNodePtr rootNSLawXMLNode;

public:

  /** \fn NonSmoothLawXML()
   *   \brief default constructor
   */
  NonSmoothLawXML();

  /** \fn NonSmoothLawXML(xmlNodePtr);
   *   \brief constructor from DOM tree data
   */
  NonSmoothLawXML(xmlNodePtr);

  /** \fn ~NonSmoothLawXML()
   *   \brief destructor
   */
  virtual ~NonSmoothLawXML();

  /** \fn const string getType() const
   *   \brief get the type of non smooth law
   *   \return a string
   */
  inline const std::string getType() const
  {
    std::string type((char*) rootNSLawXMLNode->name);
    return type;
  }

  /** \fn xmlNodePtr getNode()
   *   \brief Return node named "NonSmoothLaw"
   *   \return an xmlNodePtr
   */
  inline xmlNodePtr getNode() const
  {
    return rootNSLawXMLNode;
  }

  /** \fn void updateNonSmoothLawXML( xmlNodePtr node, NonSmoothLaw* nsl )
   *   \brief makes the operations to create the NonSmoothLaw of the Interaction
   *   \param xmlNodePtr : the root node of the NonSmoothLawXML
   *   \param Relation* : the NonSmoothLaw of this NonSmoothLawXML
   */
  void updateNonSmoothLawXML(xmlNodePtr node, NonSmoothLaw* nsl);
};

#endif
