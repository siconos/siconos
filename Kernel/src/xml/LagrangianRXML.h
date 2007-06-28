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

/*! \file LagrangianRXML.h
\brief XML data for Lagrangian Relations.
*/

#ifndef LagrangianRXML_H
#define LagrangianRXML_H

#include "RelationXML.h"
#include "SimpleMatrix.h"
#include "SimpleVector.h"

/** XML management for LagrangianR
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.0.
 *   \date 08/12/2004
 *
 */
class LagrangianRXML : public RelationXML
{

protected:

  /** plug-in name for h.*/
  xmlNodePtr hNode;

  /** plug-in names for gradients of h.*/
  std::vector<xmlNodePtr> GNode;

  /** plug-in name for hDot.*/
  xmlNodePtr hDotNode;

public:

  /** Default constructor. */
  LagrangianRXML();

  /** Build a LagrangianRXML object from a DOM tree describing a Relation with LNL type
   *   \param LagrangianRXML : the LagrangianR DOM tree
   */
  LagrangianRXML(xmlNodePtr);

  /** Destructor */
  ~LagrangianRXML();

  /** Gets the number of G defined in XML file (ie number of gradients for h).
      \return an unsigned int.
  */
  unsigned int getNumberOfGradients() const;

  /** checks if h is defined in the DOM tree
   *  \return a bool
   */
  inline bool hasH() const
  {
    return (hNode != NULL);
  };

  /** Return h plug-in name, if it exists
   *   \return a string
   */
  inline std::string getHPlugin() const
  {
    return  SiconosDOMTreeTools::getStringAttributeValue(hNode, "plugin");
  }

  /** to save the h plug-in name
   *   \param a string
   */
  inline void setHPlugin(const std::string& plugin)
  {
    if (hNode == NULL)
    {
      hNode = SiconosDOMTreeTools::createSingleNode(rootNode, "h");
      xmlNewProp(hNode, (xmlChar*)"matrixPlugin", (xmlChar*)plugin.c_str());
    }
    else
      SiconosDOMTreeTools::setStringAttributeValue(hNode, "plugin", plugin);
  }

  /** checks if hDot is defined in the DOM tree
   *  \return a bool
   */
  inline bool hasHDot() const
  {
    return (hDotNode != NULL);
  };

  /** Return true if hDot is calculated with a plugin.
   *  \return a bool
   */
  inline const bool isHDotPlugin() const
  {
    return xmlHasProp(hDotNode, (xmlChar*)"plugin");
  };

  /** Return hDot plug-in name, if it exists
   *   \return a string
   */
  inline std::string getHDotPlugin() const
  {
    if (!isHDotPlugin())
      XMLException::selfThrow("LagrangianRXML - getHDotPlugin : hDot is not plugged.");
    return  SiconosDOMTreeTools::getStringAttributeValue(hDotNode, "plugin");
  }

  /** Return hDot vector.
   *   \return a SimpleVector.
   */
  inline const SimpleVector getHDotVector() const
  {
    if (isHDotPlugin())
      XMLException::selfThrow("LagrangianRXML - getHDotVector : hDot is plugged, not set with a vector.");
    return  SiconosDOMTreeTools::getSiconosVectorValue(hDotNode);
  }

  /** to save the hDot plug-in name
   *   \param a string
   */
  void setHDotPlugin(const std::string&);

  /** to save the hDot vector.
   *   \param a SiconosVector
   */
  void setHDotVector(const SiconosVector&);

  /** Return true if G[i] is computed with a plug-in
   * \param an int: the required index for G.
   * \return a bool
   */
  bool isGPlugin(unsigned int = 0) const;

  /** checks if G[i] is defined in the DOM tree
   * \param an int: the required index for G.
   * \return a bool
   */
  bool hasG(unsigned int = 0) const;

  /** Return G plug-in name, if it exists
  * \param an int: the required index for G.
  * \return a string
  */
  std::string getGPlugin(unsigned int = 0) const ;

  /** Return the G[i] matrix.
   * \param an int: the required index for G.
   * \return a SimpleMatrix
   */
  SimpleMatrix getGMatrix(unsigned int = 0) const;

  /** to save G[i] plug-in name.
   * \param a string, the plug-in name.
   * \param an int: the required index for G.
   */
  void setGPlugin(const std::string& , unsigned int = 0);

  /** to save G[i] matrix.
   * \param a pointer to a SiconosMatrix.
   * \param an int: the required index for G.
   */
  void setGMatrix(SiconosMatrix *, unsigned int = 0);

};


#endif
