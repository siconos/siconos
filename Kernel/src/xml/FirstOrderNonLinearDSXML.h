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

/*! \file FirstOrderNonLinearDSXML.h

*/

#ifndef __FirstOrderNonLinearDSXML__
#define __FirstOrderNonLinearDSXML__

#include "DynamicalSystemXML.h"
#include "SimpleVector.h"
#include "SimpleMatrix.h"
#include "SiconosMemory.h"

//Tags
const std::string DS_X0 = "x0";
const std::string DS_X = "x";
const std::string DS_R = "R";
const std::string DS_XMEMORY = "xMemory";
const std::string DS_M = "M";
const std::string DS_F = "f";
const std::string DS_JACOBIANXF = "jacobianXF";

class SiconosMemory;
class SiconosMemoryXML;
class SimpleMatrix;
class SiconosMatrix;
class SimpleVector;
class SiconosVector;

/** XML management for FirstOrderNonLinearDS
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 04/04/2004
 *
 * Reading/Writing of xml data for the FirstOrderNonLinearDS class and its derived classes.
 *
 */
class FirstOrderNonLinearDSXML: public DynamicalSystemXML
{
protected:

  xmlNodePtr x0Node;/**< initial state */
  xmlNodePtr xNode;/**<  state (usefull is start from recovery xml-file*/
  xmlNodePtr MNode; /**< M in \f$ M \dot x = f(x,t,z) \f$ */
  xmlNodePtr fNode;/**< f(x,t) */
  xmlNodePtr jacobianXFNode;/**< jacobian of f according to x */
  xmlNodePtr xMemoryNode;/**<  memory vector for x*/
  SiconosMemoryXML * xMemoryXML;/** <xml object for xMemory*/

  /** Default constructor */
  FirstOrderNonLinearDSXML();

public:

  /** Build a FirstOrderNonLinearDSXML object from the DynamicalSystem node of the xml DOMtree
   *   \param an xmlNodePtr DynamicalSystemNode
   *   \param bool isBVP : if true, NonSmoothDynamicalSystem is a Boundary value problem
   */
  FirstOrderNonLinearDSXML(xmlNodePtr,  bool);

  /** Destructor */
  virtual ~FirstOrderNonLinearDSXML();

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
      x0Node = SiconosDOMTreeTools::createVectorNode(rootNode, DS_X0, v);
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
      xNode = SiconosDOMTreeTools::createVectorNode(rootNode, DS_X, v);
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

  // === f ===
  /** Return the name of the f plug-in
   *  \return a string
   */
  inline const std::string getFPlugin() const
  {
    if (!isFPlugin())
      XMLException::selfThrow("FirstOrderNonLinearDSXML - getFPlugin : f is not calculated from a plugin since a f vector is given.");
    return  SiconosDOMTreeTools::getStringAttributeValue(fNode, "vectorPlugin");
  }

  /** return f vector
   *   \return SimpleVector
   */
  inline const SimpleVector getFVector() const
  {
    if (isFPlugin())
      XMLException::selfThrow("FirstOrderNonLinearDSXML - getFVector : f vector is not given since f is calculated using a plug-in");
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
      XMLException::selfThrow("FirstOrderNonLinearDSXML - getJacobianXFPlugin : jacobianXF is not calculated from a plugin since a jacobianXF matrix is given.");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianXFNode, "matrixPlugin");
  }

  /** return jacobianXF matrix
   *   \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianXFMatrix() const
  {
    if (isJacobianXFPlugin())
      XMLException::selfThrow("FirstOrderNonLinearDSXML - getJacobianXFMatrix : jacobianXF matrix is not given since jacobianXF is calculated using a plug-in");
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

  /** returns true if x0Node is defined
   *  \return a bool
   */
  inline bool hasX0() const
  {
    return (x0Node != NULL);
  }

  /** returns true if xNode is defined
   *  \return a bool
   */
  inline bool hasX() const
  {
    return (xNode != NULL);
  }

  /** returns true if xMemoryNode is defined
   *  \return a bool
   */
  inline bool hasXMemory()const
  {
    return (xMemoryNode != NULL);
  }

  /** returns true if MNode is defined
   *  \return a bool
   */
  inline bool hasM() const
  {
    return (MNode != NULL);
  }

  /** returns true if fNode is defined
   *  \return a bool
   */
  inline bool hasF() const
  {
    return (fNode != NULL);
  }

  /** returns true if jacobianXFNode is defined
   *  \return a bool
   */
  inline bool hasJacobianXF() const
  {
    return (jacobianXFNode != NULL);
  }

  /** Return true if f is calculated from a plugin
   *  \return a bool
   */
  inline bool isFPlugin() const
  {
    return xmlHasProp(fNode, (xmlChar *) VECTORPLUGIN.c_str());
  }

  /** Return true if jacobianXF is calculated from a plugin
   *  \return a bool
   */
  inline bool isJacobianXFPlugin() const
  {
    return xmlHasProp(jacobianXFNode, (xmlChar *) MATRIXPLUGIN.c_str());
  }
};

#endif
