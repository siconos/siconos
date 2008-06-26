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

/*! \file RelationXML.h
  \brief XML data reading for Relation.
*/

#ifndef RelationXML_H
#define RelationXML_H

#include "SiconosDOMTreeTools.h"
#include "RelationTypes.hpp"

class Relation;

/** XML management for Relation
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 04/13/2004
 *
 */
class RelationXML
{
protected:

  /** Relation node.*/
  xmlNodePtr rootNode;

  /**Default construcor. */
  RelationXML();

  /** Copy => private: no copy nor pass-by value. */
  RelationXML(const RelationXML&);

  /**Basic construcor.
   \param a pointer to the relation xml node.
  */
  RelationXML(xmlNodePtr);

public:

  /** Destructor */
  virtual ~RelationXML();

  /** Returns the type of the RelationXML
   *  \return a string.
   */
  inline const RELATIONTYPES getType() const
  {
    std::string type((char*)rootNode->name);
    if (type == "LagrangianRelation")
      return Lagrangian;
    else if (type == "FirstOrderRelation")
      return FirstOrder;
    else
    {
      XMLException::selfThrow("RelationXML - getType: unknown type of Relation.");
      return Lagrangian;
    }
  }

  /** Returns the sub-type of the Relation
   *  \return a string.
   */
  const RELATIONSUBTYPES getSubType() const
  {
    std::string res = SiconosDOMTreeTools::getStringAttributeValue(rootNode, "type");
    if (res == "NonLinear")
      return NonLinearR;
    else if (res == "Linear")
      return LinearR;
    else if (res == "Type1")
      return Type1R;
    else if (res == "LinearTI")
      return LinearTIR;
    else if (res == "Scleronomous")
      return ScleronomousR;
    else if (res == "Rheonomous")
      return RheonomousR;
    else if (res == "Compliant")
      return CompliantR;
    else
    {
      XMLException::selfThrow("RelationXML - getType: unknown type of Relation.");
      return NonLinearR;
    }
  }

  /** Returns the node of the RelationXML
   *   \return an xmlNodePtr.
   */
  inline xmlNodePtr getRootNode()const
  {
    return rootNode;
  }

  /** To create the Relation of the Interaction. Useless?
   *   \param xmlNodePtr : the root node of the RelationXML
   *   \param Relation* : the Relation of this RelationXML
   */
  void updateRelationXML(xmlNodePtr, Relation*);
};

#endif
