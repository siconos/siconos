/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

/*! \file LinearDSXML.h

*/


#ifndef __LINEARSYSTEMDSXML__
#define __LINEARSYSTEMDSXML__

#include "DynamicalSystemXML.h"

const std::string LDS_A = "A";
const std::string LDS_Mxdot = "Mxdot";
const std::string LDS_B = "b";

class SimpleMatrix;
/** XML management for LinearDS
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.0.1.
 *   \date 05/11/2004
 *
 *
 *
 * LinearDSXML allows to manage data of a LinearDS DOM tree.
 */
class LinearDSXML : public DynamicalSystemXML
{
private:

  //Nodes
  xmlNode * ANode; /**< A in \f$ Mxdot \dot x = Ax+b \f$ */
  xmlNode * MxdotNode; /**< Mxdot in \f$ Mxdot \dot x = Ax+b \f$ */
  xmlNode * bNode; /**< b in \f$ Mxdot \dot x = Ax+b \f$ */

public:

  /** Default constructor
  */
  LinearDSXML();

  /** Build a LinearDSXML object from a DOM tree describing a DS
  *   \param xmlNode * linearSystemDSNode : the linearSystemDS DOM tree
  *   \param bool isBVP : if NonSmoothDynamicalSystem is BVP, linearSystemDS has boundary condition
  */
  LinearDSXML(xmlNode * linearSystemDSNode, const bool& isBVP);

  /** destructor
  */
  ~LinearDSXML();

  /** return the A of the LinearDSXML
  *   \return a SimpleMatrix
  */
  inline const SimpleMatrix getA() const
  {
    if (isAPlugin())
      XMLException::selfThrow("LinearDSXML - getA: A is not given but calculated with a plugin");
    return  SiconosDOMTreeTools::getSiconosMatrixValue(ANode);
  }

  /** return the A plug-in name of the LinearDSXML
  *   \return a string
  */
  inline const std::string getAPlugin() const
  {
    if (!isAPlugin())
      XMLException::selfThrow("LinearDSXML - getAPlugin : A is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(ANode, DS_MATRIXPLUGIN);
  }

  /** to save the A of the LinearDSXML
  *   \param The A SiconosMatrix to save
  */
  void setA(const SiconosMatrix& m);

  /** to save the A plugin
  *   \param a string (name of the plug-in)
  */
  void setAPlugin(const std::string& plugin);

  /** return the optional matrix Mxdot of the LinearDSXML
  *   \return a SimpleMatrix
  */
  inline const SimpleMatrix getMxdot() const
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(MxdotNode);
  }

  /** to save the optional Mxdot of the LinearDSXML
  *   \param The Mxdot SiconosMatrix to save
  */
  void setMxdot(const SiconosMatrix& m);

  /** Return the b plug-in name of the LinearDSXML
  *   \return a string
  */
  inline const std::string getBPlugin() const
  {
    if (!isBPlugin())
      XMLException::selfThrow("LinearDSXML - getUPlugin : b is not calculated from a plugin ; b vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(bNode, DS_VECTORPLUGIN);
  }

  /** Return b vector of the LinearDSXML
  *   \return a SimpleVector
  */
  inline const SimpleVector getBVector() const
  {
    if (isBPlugin())
      XMLException::selfThrow("LinearDSXML - getBVector : b vector is not given ; b is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(bNode);
  }

  /** to save the b vector of the LinearDSXML
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
    return xmlHasProp((xmlNodePtr)ANode, (xmlChar *) DS_MATRIXPLUGIN.c_str());
  }

  /** Return true if b is calculated from a plugin
  */
  inline bool isBPlugin() const
  {
    return xmlHasProp((xmlNodePtr)bNode, (xmlChar *) DS_VECTORPLUGIN.c_str());
  }

  /** returns true if ANode is defined
  *  \return true if ANode is defined
  */
  inline bool hasA() const
  {
    return (ANode != NULL);
  }

  /** returns true if MxdotNode is defined
  *  \return true if MxdotNode is defined
  */
  inline bool hasMxdot() const
  {
    return (MxdotNode != NULL);
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
