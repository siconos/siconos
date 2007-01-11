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

/*! \file DynamicalSystemXML.h

*/

#ifndef __DynamicalSystemXML__
#define __DynamicalSystemXML__

#include "SiconosDOMTreeTools.h"

#include "SimpleVector.h"
#include "SimpleMatrix.h"
#include "SiconosMemoryXML.h"

class SiconosMemory;
class DynamicalSystem;
class SiconosMemoryXML;

//Tags
const std::string DS_X0 = "x0";
const std::string DS_X = "x";
const std::string DS_R = "R";
const std::string DS_U = "u";
const std::string DS_T = "T";
const std::string DS_XMEMORY = "xMemory";
const std::string DS_RMEMORY = "RMemory";
const std::string DS_STEPSINMEMORY = "StepsInMemory";
const std::string DS_F = "f";
const std::string DS_M = "M";
const std::string DS_JACOBIANXF = "jacobianXF";
const std::string DS_MATRIXPLUGIN = "matrixPlugin";
const std::string DS_VECTORPLUGIN = "vectorPlugin";

class SimpleMatrix;
class SiconosMatrix;
class SimpleVector;
class SiconosVector;

/** XML management for DynamicalSystem
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.0.1.
 *   \date 04/04/2004
 *
 *
 * DynamicalSystemXML allows to manage data of a DynamicalSystem DOM tree.
 */
class DynamicalSystemXML
{

protected:

  xmlNodePtr rootDynamicalSystemXMLNode;/**< root node  */
  xmlNodePtr parentNode; /**< node that owns root node (NonSmoothDynamicalSystemNode)  */
  xmlNodePtr x0Node;/**< initial state */
  xmlNodePtr xNode;/**<  state (usefull is start from recovery xml-file*/
  xmlNodePtr stepsInMemoryNode; /**< size of memory */
  xmlNodePtr xMemoryNode;/**<  memory vector for x*/
  xmlNodePtr rMemoryNode;/**< memory vector for r */
  xmlNodePtr MNode; /**< M in \f$ M \dot x = Ax+b \f$ */
  xmlNodePtr fNode;/**< f(x,t) */
  xmlNodePtr jacobianXFNode;/**< jacobian of f according to x */
  xmlNodePtr uSizeNode;/**< size of control vector */
  xmlNodePtr uNode;/**< control term */
  xmlNodePtr TNode;/**< coefficient of control term */

  SiconosMemoryXML * xMemoryXML;/** <xml object for xMemory*/
  SiconosMemoryXML * rMemoryXML;/**<xml object for rMemory*/

public:

  /** Default constructor */
  DynamicalSystemXML();

  /** Destructor */
  virtual ~DynamicalSystemXML();

  /** Build a DynamicalSystemXML object from the DynamicalSystem node of the xml DOMtree
   *   \param an xmlNodePtr DynamicalSystemNode
   *   \param bool isBVP : if true, NonSmoothDynamicalSystem is a Boundary value problem
   */
  DynamicalSystemXML(xmlNodePtr DynamicalSystemNode, const bool&);

  /** Return the number of the DynamicalSystem
   *   \return an integer
   */
  inline const int getNumber() const
  {
    return SiconosDOMTreeTools::getAttributeValue<int>(rootDynamicalSystemXMLNode, NUMBER_ATTRIBUTE);
  }

  /** Return the type of the DynamicalSystem
   *   \return a string
   */
  inline const std::string getType() const
  {
    std::string res((char*)rootDynamicalSystemXMLNode->name);
    return res;
  }

  /** Return the id of the DynamicalSystem
   *   \return The string id of the DynamicalSystem
   */
  inline const std::string getId() const
  {
    return SiconosDOMTreeTools::getStringAttributeValue(rootDynamicalSystemXMLNode, "Id");
  }

  /** return true if id is given
   *  \return a bool
   */
  inline bool hasId() const
  {
    return (xmlHasProp(rootDynamicalSystemXMLNode, (xmlChar*)"Id"));
  }

  /** to save the id of the DynamicalSystem
   *   \param a string
   */
  inline void setId(const std::string& s)
  {
    SiconosDOMTreeTools::setStringAttributeValue(rootDynamicalSystemXMLNode, "Id", s);
  }

  /** Returns initial state of the DynamicalSystem (x0)
   *   \return a SimpleVector
   */
  inline const SimpleVector getX0() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(x0Node);
  }

  /** save x0 of the DynamicalSystem
   *   \param a SiconosVector
   */
  inline void setX0(const SiconosVector& v)
  {
    if (!hasX0())
      x0Node = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, DS_X0, v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(x0Node, v);
  }

  /** Returns the x state-vector of the DynamicalSystem
   *   \return SimpleVector
   */
  inline const SimpleVector getX() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(xNode);
  }

  /** save x of the DynamicalSystem
   *   \param a SiconosVector
   */
  inline void setX(const SiconosVector& v)
  {
    if (!hasX())
      xNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, DS_X, v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(xNode, v);
  }

  /** return the optional matrix M of the LinearDSXML
   *   \return a SimpleMatrix
   */
  inline const SimpleMatrix getM() const
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(MNode);
  }

  /** to save the optional M of the LinearDSXML
   *   \param The M SiconosMatrix to save
   */
  void setM(const SiconosMatrix& m);

  /** Returns the steps in memory for the DynamicalSystemXML
   *   \return The integer number of steps in memory for the DynamicalSystemXML
   */
  inline const unsigned int getStepsInMemory() const
  {
    return  SiconosDOMTreeTools::getContentValue<unsigned int>(stepsInMemoryNode);
  }

  /** to save the steps in memory for the DynamicalSystemXML
   *   \param an unsigned int
   */
  void setStepsInMemory(const unsigned int&);

  /** Returns the xMemoryXML* of the DynamicalSystemXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXML* getXMemoryXML() const
  {
    return xMemoryXML;
  }

  /** to save the XMemory of the DynamicalSystemXML
   *   \param SiconosMemory* smem : SiconosMemory to save
   */
  void setXMemory(const SiconosMemory&);

  /** Returns the rMemoryXML* of the DynamicalSystemXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXML* getRMemoryXML() const
  {
    return rMemoryXML;
  }

  /** to save the rMemory of the DynamicalSystemXML
   *   \param SiconosMemory*
   */
  void setRMemory(const SiconosMemory&);

  // === f ===
  /** Return the name of the f plug-in
   *  \return a string
   */
  inline const std::string getFPlugin() const
  {
    if (!isFPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getFPlugin : f is not calculated from a plugin since a f vector is given.");
    return  SiconosDOMTreeTools::getStringAttributeValue(fNode, "vectorPlugin");
  }

  /** return f vector
   *   \return SimpleVector
   */
  inline const SimpleVector getFVector() const
  {
    if (isFPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getFVector : f vector is not given since f is calculated using a plug-in");
    return  SiconosDOMTreeTools::getSiconosVectorValue(fNode);
  }

  /** to save the f vector
   *   \param a SiconosVector
   */
  void setFVector(const SiconosVector&v);

  /** to save the F plugin
   *   \param a string (name of the plug-in)
   */
  void setFPlugin(const std::string& plugin);

  // === JacobianXF ===
  /** Return the name of the jacobianXF plug-in
   *  \return a string
   */
  inline const std::string getJacobianXFPlugin() const
  {
    if (!isJacobianXFPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getJacobianXFPlugin : jacobianXF is not calculated from a plugin since a jacobianXF matrix is given.");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianXFNode, "matrixPlugin");
  }

  /** return jacobianXF matrix
   *   \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianXFMatrix() const
  {
    if (isJacobianXFPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getJacobianXFMatrix : jacobianXF matrix is not given since jacobianXF is calculated using a plug-in");
    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianXFNode);
  }

  /** to save the jacobianXF matrix
   *   \param a SiconosMatrix
   */
  void setJacobianXFMatrix(const SiconosMatrix&v);

  /** to save the JacobianXF plugin of the LagrangianDSXML
   *   \param a string (name of the plug-in)
   */
  void setJacobianXFPlugin(const std::string& plugin);

  /** get size of vector u
   */
  inline const unsigned int getUSize() const
  {
    if (!hasUSize())
      XMLException::selfThrow("DynamicalSystemXML - getUSize: this node does not exist");
    return  SiconosDOMTreeTools::getContentValue<unsigned int>(uSizeNode);
  }

  /** to save the size of vector u
   *   \param an unsigned int
   */
  inline void setUSize(const unsigned int& us)
  {
    if (!hasUSize())
      uSizeNode = SiconosDOMTreeTools::createIntegerNode(rootDynamicalSystemXMLNode, "uSize", us);
    else SiconosDOMTreeTools::setIntegerContentValue(uSizeNode, us);
  }

  // --- u ---

  /** Return the u plug-in name
   *   \return a string
   */
  inline const std::string getUPlugin() const
  {
    if (!isUPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getUPlugin : u is not calculated from a plugin ; u vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(uNode, DS_VECTORPLUGIN);
  }

  /** Return u vector
   *   \return a SimpleVector
   *  \exception XMLException
   */
  inline const SimpleVector getUVector() const
  {
    if (isUPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getUVector : u vector is not given ; u is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(uNode);
  }

  /** to save the u vector
   *   \param a SiconosVector
   */
  void setUVector(const SiconosVector&);

  /** to save the u plugin
   *   \param a string (name of the plug-in)
   */
  void setUPlugin(const std::string& plugin);

  // --- T ---

  /** get T Matrix
   *   \return a SimpleMatrix
   */
  inline const SimpleMatrix getTMatrix() const
  {
    if (isTPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getT: T is not given but calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosMatrixValue(TNode);
  }

  /** get Plugin name to compute T
   *  \exception XMLException
   */
  inline const std::string getTPlugin() const
  {
    if (!isTPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getTPlugin : T is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(TNode, DS_MATRIXPLUGIN);
  }

  /** save T
   *   \param The SiconosMatrix to save
   */
  void setT(const SiconosMatrix &m);

  /** to save the T plugin
   *   \param a string (name of the plug-in)
   */
  void setTPlugin(const std::string& plugin);

  /** returns true if x0Node is defined
   *  \return true if x0Node is defined
   */
  inline bool hasX0() const
  {
    return (x0Node != NULL);
  }

  /** returns true if xNode is defined
   *  \return true if xNode is defined
   */
  inline bool hasX() const
  {
    return (xNode != NULL);
  }

  /** returns true if stepsInMemoryNode is defined
   *  \return true if stepsInMemoryNode is defined
   */
  inline bool hasStepsInMemory() const
  {
    return (stepsInMemoryNode != NULL);
  }

  /** returns true if xMemoryNode is defined
   *  \return true if xMemoryNode is defined
   */
  inline bool hasXMemory()const
  {
    return (xMemoryNode != NULL);
  }

  /** returns true if RMemory is defined
   *  \return true if RMemory is defined
   */
  inline bool hasRMemory() const
  {
    return (rMemoryNode != NULL);
  }

  /** returns true if MNode is defined
   *  \return true if MNode is defined
   */
  inline bool hasM() const
  {
    return (MNode != NULL);
  }

  /** returns true if fNode is defined
   *  \return true if fNode is defined
   */
  inline bool hasF() const
  {
    return (fNode != NULL);
  }

  /** returns true if jacobianXFNode is defined
   *  \return true if jacobianXFNode is defined
   */
  inline bool hasJacobianXF() const
  {
    return (jacobianXFNode != NULL);
  }

  /** returns true if uSizeNode is defined
   *  \return true if jacobianXFNode is defined
   */
  inline bool hasUSize() const
  {
    return (uSizeNode != NULL);
  }

  /** returns true if uNode is defined
   *  \return true if jacobianXFNode is defined
   */
  inline bool hasU() const
  {
    return (uNode != NULL);
  }

  /** returns true if TNode is defined
   *  \return true if jacobianXFNode is defined
   */
  inline bool hasT() const
  {
    return (TNode != NULL);
  }

  /** Return true if f is calculated from a plugin
   */
  inline bool isFPlugin() const
  {
    return xmlHasProp((xmlNodePtr)fNode, (xmlChar *) DS_VECTORPLUGIN.c_str());
  }

  /** Return true if jacobianXF is calculated from a plugin
   */
  inline bool isJacobianXFPlugin() const
  {
    return xmlHasProp((xmlNodePtr)jacobianXFNode, (xmlChar *) DS_MATRIXPLUGIN.c_str());
  }

  /** Return true if u is calculated from a plugin
   */
  inline bool isUPlugin() const
  {
    return xmlHasProp((xmlNodePtr)uNode, (xmlChar *) DS_VECTORPLUGIN.c_str());
  }

  /** Return true if T is calculated from a plugin
   */
  inline bool isTPlugin() const
  {
    return xmlHasProp((xmlNodePtr)TNode, (xmlChar *) DS_MATRIXPLUGIN.c_str());
  }
};

#endif
