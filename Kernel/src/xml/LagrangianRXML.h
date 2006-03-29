/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
 *   \version 1.1.4.
 *   \date 08/12/2004
 *
 *
 *
 * LagrangianRXML allows to manage data of a LNLRelation DOM tree.
 */

#ifndef __LAGRANGIANRelationXML__
#define __LAGRANGIANRelationXML__

#include "RelationXML.h"
#include <vector>

const std::string LAGRANGIANR_MATRIXPLUGIN = "matrixPlugin";

class LagrangianRXML : public RelationXML
{

protected:
  //Nodes

  xmlNodePtr hNode; // node for h function (plug-in)
  std::vector<xmlNodePtr> GNode; // node for G function(s) (plug-in or matrix)

public:

  LagrangianRXML();

  /** \fn LagrangianRXML(xmlNode * LNLRelationNode)
   *   \brief Build a LagrangianRXML object from a DOM tree describing a Relation with LNL type
   *   \param LagrangianRXML : the LagrangianR DOM tree
   *   \exception XMLException : if a property of the LagrangianLinear Relation lacks in the DOM tree
   */
  LagrangianRXML(xmlNode * LNLRelationNode);

  ~LagrangianRXML();

  /** \fn std::string  getLagrangianType() const;
   *   \brief Returns the type of constraints in the Lagrangian relation
   *   \return a string
   */
  inline std::string getLagrangianType() const
  {
    if (SiconosDOMTreeTools::hasAttributeValue(rootRelationXMLNode, "type"))
      return SiconosDOMTreeTools::getStringAttributeValue(rootRelationXMLNode, "type");
    else return "scleronomic";
  }

  // === h matrix/Plug-in ===

  /** \fn bool hasH()
   *  \brief checks if h is defined in the DOM tree
   *  \return a bool
   */
  inline bool hasH()
  {
    return (hNode != NULL);
  } const

  /** \fn inline string getHPlugin()
   *   \brief Return h plug-in name, if it exists
   *   \return a string
   *  \exception XMLException
   */
  inline std::string getHPlugin() const
  {
    return  SiconosDOMTreeTools::getStringAttributeValue(hNode, PLUGIN_ATTRIBUTE);
  }

  /** \fn void setHPlugin(const string& plugin)
   *   \brief to save the h plug-in name
   *   \param a string
   */
  inline void setHPlugin(const std::string& plugin)
  {
    if (hNode == NULL)
    {
      hNode = SiconosDOMTreeTools::createSingleNode(rootRelationXMLNode, "h");
      xmlNewProp(hNode, (xmlChar*)(LAGRANGIANR_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else
      SiconosDOMTreeTools::setStringAttributeValue(hNode, PLUGIN_ATTRIBUTE, plugin);
  }

  // === G matrix/Plug-in ===

  /** \fn bool isGPlugin(const unsigned int & = 0)
   *   \brief Return true if G is computed with a plug-in
   *   \return a bool
   */
  bool isGPlugin(const unsigned int & = 0) const;

  /** \fn bool hasG(const unsigned int & = 0)
   *  \brief checks if G is defined in the DOM tree
   *  \return a bool
   */
  bool hasG(const unsigned int & = 0) const;

  /** \fn inline string getGPlugin(const unsigned int & = 0)
   *   \brief Return G plug-in name, if it exists
   *   \return a string
   *  \exception XMLException
   */
  std::string getGPlugin(const unsigned int & = 0) const ;

  /** \fn SimpleMatrix getGMatrix(const unsigned int & = 0)
   *   \brief Return the G matrix of the LagrangianRXML
   *   \return a SimpleMatrix
   *  \exception XMLException
   */
  SimpleMatrix getGMatrix(const unsigned int & = 0) const;

  /** \fn void setGPlugin(const string& plugin,const unsigned int & = 0)
   *   \brief to save the G plug-in name
   *   \param a string
   */
  void setGPlugin(const std::string& plugin, const unsigned int & = 0);

  /** \fn void setGMatrix(SiconosMatrix *newMat,const unsigned int & = 0)
   *   \brief to save the G matrix of the LagrangianRXML
   *   \return a pointer to SiconosMatrix
   */
  void setGMatrix(SiconosMatrix *newMat, const unsigned int & = 0);

};


#endif
