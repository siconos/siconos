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

/** \class LinearDSXML
 *   \brief This class manages LinearSystem DS data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.1.3.
 *   \date 05/11/2004
 *
 *
 *
 * LinearDSXML allows to manage data of a LinearDS DOM tree.
 */

#ifndef __LINEARSYSTEMDSXML__
#define __LINEARSYSTEMDSXML__

#include "DynamicalSystemXML.h"

const std::string LDS_A = "A";
const std::string LDS_B = "b";
const std::string LDS_E = "E";

class LinearDSXML : public DynamicalSystemXML
{
public:
  LinearDSXML();

  /** \fn LinearDSXML(xmlNode * LagrangianDSNode, bool isBVP)
   *   \brief Build a LinearDSXML object from a DOM tree describing a DS
   *   \param xmlNode * linearSystemDSNode : the linearSystemDS DOM tree
   *   \param bool isBVP : if NonSmoothDynamicalSystem is BVP, linearSystemDS has boundary condition
   */
  LinearDSXML(xmlNode * linearSystemDSNode, const bool& isBVP);

  ~LinearDSXML();

  /** \fn SiconosMatrix getA()
   *   \brief Return the A of the LinearDSXML
   *   \return The A SiconosMatrix of the LinearDSXML
   */
  inline const SiconosMatrix getA() const
  {
    if (isAPlugin())
      XMLException::selfThrow("LinearDSXML - getA: A is not given but calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosMatrixValue(ANode);
  }

  /** \fn inline string getAPlugin()
   *   \brief Return the A Plugin name of the LinearDSXML
   *  \exception XMLException
   */
  inline const std::string getAPlugin() const
  {
    if (!isAPlugin())
      XMLException::selfThrow("LinearDSXML - getAPlugin : A is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(ANode, DS_MATRIXPLUGIN);
  }

  /** \fn void setA(SiconosMatrix *m)
   *   \brief allows to save the A of the LinearDSXML
   *   \return The A SiconosMatrix to save
   */
  inline void setA(const SiconosMatrix& m)
  {
    if (ANode != NULL)
      SiconosDOMTreeTools::setSiconosMatrixNodeValue(ANode, m);
    else
      ANode = SiconosDOMTreeTools::createMatrixNode(rootDynamicalSystemXMLNode, LDS_A, m);
  }

  /** \fn SiconosMatrix getE()
   *   \brief Return the E of the LinearDSXML
   *   \return The E SiconosMatrix of the LinearDSXML
   */
  inline const SiconosMatrix getE() const
  {
    if (isEPlugin())
      XMLException::selfThrow("LinearDSXML - getE: E is not given but calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosMatrixValue(ENode);
  }

  /** \fn inline string getEPlugin()
   *   \brief Return the E Plugin name of the LinearDSXML
   *  \exception XMLException
   */
  inline const std::string getEPlugin() const
  {
    if (!isEPlugin())
      XMLException::selfThrow("LinearDSXML - getEPlugin : E is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(ENode, DS_MATRIXPLUGIN);
  }

  /** \fn void setE(SiconosMatrix *m)
   *   \brief allows to save the E of the LinearDSXML
   *   \param The SiconosMatrix to save
   */
  inline void setE(const SiconosMatrix &m)
  {
    if (ENode != NULL)
      SiconosDOMTreeTools::setSiconosMatrixNodeValue(ENode, m);
    else ENode = SiconosDOMTreeTools::createMatrixNode(rootDynamicalSystemXMLNode, LDS_E, m);
  }

  /** \fn inline string getBPlugin()
   *   \brief Return the b Plugin name of the LinearDSXML
   *   \return The b Plugin name of the LinearDSXML
   *  \exception XMLException
   */
  inline const std::string getBPlugin() const
  {
    if (!isBPlugin())
      XMLException::selfThrow("LinearDSXML - getUPlugin : b is not calculated from a plugin ; b vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(bNode, DS_VECTORPLUGIN);
  }

  /** \fn inline SimpleVector getBVector()
   *   \brief Return b vector of the LinearDSXML
   *   \return SimpleVector : value of b of LinearDSXML
   *  \exception XMLException
   */
  inline const SimpleVector getBVector() const
  {
    if (isBPlugin())
      XMLException::selfThrow("LinearDSXML - getBVector : b vector is not given ; b is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(bNode);
  }

  /** \fn inline void setBVector(SiconosVector *v)
   *   \brief allows to save the b vector of the LinearDSXML
   *   \return The b SimpleVector to save
   */
  inline void setBVector(const SiconosVector& v)
  {
    if (bNode != NULL)
      SiconosDOMTreeTools::setSiconosVectorNodeValue(bNode, v);
    else bNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LDS_B, v);
  }

  /** \fn bool isAPlugin()
   *   \brief Return true if A is calculated from a plugin
   */
  inline bool isAPlugin() const
  {
    return xmlHasProp((xmlNodePtr)ANode, (xmlChar *) DS_MATRIXPLUGIN.c_str());
  }

  /** \fn bool isBPlugin()
   *   \brief Return true if b is calculated from a plugin
   */
  inline bool isBPlugin() const
  {
    return xmlHasProp((xmlNodePtr)bNode, (xmlChar *) DS_VECTORPLUGIN.c_str());
  }

  /** \fn bool isEPlugin()
   *   \brief Return true if E is calculated from a plugin
   */
  inline bool isEPlugin() const
  {
    return xmlHasProp((xmlNodePtr)ENode, (xmlChar *) DS_MATRIXPLUGIN.c_str());
  }

  /** \fn bool hasXX()
   * \brief return true if XXnode exists */
  inline bool hasA() const
  {
    return (ANode != NULL);
  }
  inline bool hasB() const
  {
    return (bNode != NULL);
  }
  inline bool hasE() const
  {
    return (ENode != NULL);
  }


  /** \fn void updateDynamicalSystemXML(xmlNode*, DynamicalSystem*, BoundaryCondition*)
   *   \brief makes the operations to add a DynamicalSystem to the NonSmoothDynamicalSystemXML
   *   \param xmlNode* : the root node of this DynamicalSystem
   *   \param DynamicalSystem* : the DynamicalSystem of this DynamicalSystemXML
   *   \param BoundaryCondition* : the BoundaryCondition of the DS if the NonSmoothDynamicalSystem is BVP (optional)
   */
  void updateDynamicalSystemXML(xmlNode*, DynamicalSystem*, BoundaryCondition* bc = NULL);


private:

  //Nodes
  xmlNode * ANode;
  xmlNode * bNode;
  xmlNode * ENode;
};

#endif
