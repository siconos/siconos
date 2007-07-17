/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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

/*! \file FirstOrderLinearDSXML.h

*/


#ifndef __LINEARSYSTEMDSXML__
#define __LINEARSYSTEMDSXML__

#include "FirstOrderNonLinearDSXML.h"

const std::string LDS_A = "A";
const std::string LDS_B = "b";

class SimpleMatrix;

/** XML management for FirstOrderLinearDS
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.1.
 *   \date 05/11/2004
 *
 *
 *
 * FirstOrderLinearDSXML allows to manage data of a FirstOrderLinearDS DOM tree.
 */
class FirstOrderLinearDSXML : public FirstOrderNonLinearDSXML
{
private:

  //Nodes
  xmlNodePtr ANode; /**< A in \f$ M \dot x = Ax+b \f$ */
  xmlNodePtr bNode; /**< b in \f$ M \dot x = Ax+b \f$ */

  /** Default constructor
   */
  FirstOrderLinearDSXML();

public:

  /** Build a FirstOrderLinearDSXML object from a DOM tree describing a DS
  *   \param xmlNode * linearSystemDSNode : the linearSystemDS DOM tree
  *   \param bool isBVP : if NonSmoothDynamicalSystem is BVP, linearSystemDS has boundary condition
  */
  FirstOrderLinearDSXML(xmlNodePtr, bool);

  /** destructor
  */
  inline ~FirstOrderLinearDSXML() {};

  /** return the A of the FirstOrderLinearDSXML
  *   \return a SimpleMatrix
  */
  inline const SimpleMatrix getA() const
  {
    if (isAPlugin())
      XMLException::selfThrow("FirstOrderLinearDSXML - getA: A is not given but calculated with a plugin");
    return  SiconosDOMTreeTools::getSiconosMatrixValue(ANode);
  }

  /** return the A plug-in name of the FirstOrderLinearDSXML
  *   \return a string
  */
  inline const std::string getAPlugin() const
  {
    if (!isAPlugin())
      XMLException::selfThrow("FirstOrderLinearDSXML - getAPlugin : A is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(ANode, MATRIXPLUGIN);
  }

  /** to save the A of the FirstOrderLinearDSXML
  *   \param The A SiconosMatrix to save
  */
  void setA(const SiconosMatrix& m);

  /** to save the A plugin
  *   \param a string (name of the plug-in)
  */
  void setAPlugin(const std::string& plugin);

  /** Return the b plug-in name of the FirstOrderLinearDSXML
  *   \return a string
  */
  inline const std::string getBPlugin() const
  {
    if (!isBPlugin())
      XMLException::selfThrow("FirstOrderLinearDSXML - getBPlugin : b is not calculated from a plugin ; b vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(bNode, VECTORPLUGIN);
  }

  /** Return b vector of the FirstOrderLinearDSXML
  *   \return a SimpleVector
  */
  inline const SimpleVector getBVector() const
  {
    if (isBPlugin())
      XMLException::selfThrow("FirstOrderLinearDSXML - getBVector : b vector is not given ; b is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(bNode);
  }

  /** to save the b vector of the FirstOrderLinearDSXML
  *   \param The b SimpleVector to save
  */
  void setB(const SiconosVector& v);

  /** to save the B plugin
  *   \param a string (name of the plug-in)
  */
  void setBPlugin(const std::string& plugin);

  /** Return true if A is calculated from a plugin
  */
  inline bool isAPlugin() const
  {
    return xmlHasProp(ANode, (xmlChar *) MATRIXPLUGIN.c_str());
  }

  /** Return true if b is calculated from a plugin
  */
  inline bool isBPlugin() const
  {
    return xmlHasProp(bNode, (xmlChar *) VECTORPLUGIN.c_str());
  }

  /** returns true if ANode is defined
  *  \return true if ANode is defined
  */
  inline bool hasA() const
  {
    return (ANode != NULL);
  }

  /** returns true if bNode is defined
  *  \return true if bNode is defined
  */
  inline bool hasB() const
  {
    return (bNode != NULL);
  }

};

#endif
