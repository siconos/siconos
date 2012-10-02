/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

/*! \file NonSmoothLawXML.hpp
  \brief

*/

#ifndef __NSLawXML__
#define __NSLawXML__

#include "SiconosDOMTreeTools.hpp"

class NonSmoothLaw;

/** XML management for NonSmoothLaw
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 05/13/2004
 *
 *
 *
 * NSLawXML allows to manage data of a NSLaw DOM tree.
 */
class NonSmoothLawXML
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NonSmoothLawXML);


  /** rootNode */
  xmlNodePtr rootNode;

  /** size */
  xmlNodePtr sizeNode;

public:

  /** default constructor
  */
  NonSmoothLawXML();

  /** constructor from DOM tree data
  */
  NonSmoothLawXML(xmlNodePtr);

  /** destructor
  */
  virtual ~NonSmoothLawXML();

  /** get the type of non smooth law
  *   \return a string
  */
  inline const std::string getType() const
  {
    std::string type((char*) rootNode->name);
    return type;
  }

  /** Return node named "NonSmoothLaw"
  *   \return an xmlNodePtr
  */
  inline xmlNodePtr getNode() const
  {
    return rootNode;
  }

  /** return true if size node is present
  *  \return a bool
  */
  bool hasSize() const
  {
    return (sizeNode);
  };

  /** Return the size of the InteractionXML
  *   \return an unsigned int
  */
  inline unsigned int getSize() const
  {
    if (!hasSize())
      XMLException::selfThrow("NonSmoothLawXML::getSize() : sizeNode == NULL");
    return SiconosDOMTreeTools::getContentValue<int>(sizeNode);
  }

  /** to save the size of the Interaction
  *  \return an unsigned int
  */
  void setSize(const unsigned int);

  /** makes the operations to create the NonSmoothLaw of the Interaction
  *   \param xmlNodePtr : the root node of the NonSmoothLawXML
  *   \param Relation* : the NonSmoothLaw of this NonSmoothLawXML
  */
  void updateNonSmoothLawXML(xmlNodePtr node, NonSmoothLaw* nsl);
};
#endif
