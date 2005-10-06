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

/** \class LagrangianRXML
 *   \brief This class manages Lagrangian Relation data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 08/12/2004
 *
 *
 *
 * LagrangianRXML allows to manage data of a LNLRelation DOM tree.
 */

#ifndef __LAGRANGIANRelationXML__
#define __LAGRANGIANRelationXML__

#include "RelationXML.h"

class LagrangianRXML : public RelationXML
{
public:

  LagrangianRXML();

  /** \fn LagrangianRXML(xmlNode * LNLRelationNode)
   *   \brief Build a LagrangianRXML object from a DOM tree describing a Relation with LNL type
   *   \param LagrangianRXML : the LagrangianR DOM tree
   *   \exception XMLException : if a property of the LagrangianLinear Relation lacks in the DOM tree
   */
  LagrangianRXML(xmlNode * LNLRelationNode);

  ~LagrangianRXML();

  /** \fn int getComputeInputPlugin()
   *   \brief Returns the computeInput plugin of the Relation
   *   \return string which defines the plugin
   */
  std::string  getComputeInputPlugin() const;

  /** \fn int getComputeOutputPlugin()
   *   \brief Returns the computeOutput plugin of the Relation
   *   \return string which defines the plugin
   */
  std::string  getComputeOutputPlugin() const;

private:
  //Nodes
};


#endif
