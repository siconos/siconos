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

using RELATION::TYPES;
using RELATION::SUBTYPES;

class RelationXML
{
private:
  /** Copy => private: no copy nor pass-by value. */
  RelationXML(const RelationXML&);

protected:

  /** Relation node.*/
  xmlNodePtr rootNode;

  /**Default constructor. */
  RelationXML(): rootNode(NULL) {};

  /**Basic construcor.
   \param a pointer to the relation xml node.
  */
  RelationXML(xmlNodePtr node): rootNode(node) {};

public:

  /** Destructor */
  virtual ~RelationXML() {};

  /** Returns the type of the RelationXML
   */
  const RELATION::TYPES getType() const ;

  /** Returns the sub-type of the Relation
   *  \return a string.
   */
  const RELATION::SUBTYPES getSubType() const;

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
  inline void updateRelationXML(xmlNodePtr node, Relation*)
  {
    rootNode = node;
  }

};

#endif
