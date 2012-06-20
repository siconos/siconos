/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*! \file FirstOrderNonLinearDSXML.hpp

*/

#ifndef __FirstOrderNonLinearDSXML__
#define __FirstOrderNonLinearDSXML__

#include "DynamicalSystemXML.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosMemory.hpp"

//Tags
const std::string DS_X0 = "x0";
const std::string DS_X = "x";
const std::string DS_R = "R";
const std::string DS_XMEMORY = "xMemory";
const std::string DS_M = "M";
const std::string DS_F = "f";
const std::string DS_JACOBIANXF = "jacobianfx";

class SiconosMemory;
class SiconosMemoryXML;
class SimpleMatrix;
class SiconosMatrix;
class SiconosVector;
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
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderNonLinearDSXML);


  xmlNodePtr x0Node;/**< initial state */
  xmlNodePtr xNode;/**<  state (usefull is start from recovery xml-file*/
  xmlNodePtr MNode; /**< M in \f$ M \dot x = f(x,t,z) \f$ */
  xmlNodePtr fNode;/**< f(x,t) */
  xmlNodePtr jacobianfxNode;/**< jacobian of f according to x */
  xmlNodePtr xMemoryNode;/**<  memory vector for x*/
  SP::SiconosMemoryXML xMemoryXML;/** <xml object for xMemory*/

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
   *   \return a SiconosVector
   */
  inline const SiconosVector getX0() const
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
   *   \return SiconosVector
   */
  inline const SiconosVector getx() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(xNode);
  }

  /** save x of the DynamicalSystem
   *   \param a SiconosVector
   */
  inline void setX(const SiconosVector& v)
  {
    if (!hasx())
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

  /** Return the name of the M plug-in
   *  \return a string
   */
  inline const std::string getMPlugin() const
  {
    if (!isMPlugin())
      XMLException::selfThrow("FirstOrderNonLinearDSXML - getMPlugin : jacobianfx is not calculated from a plugin since a jacobianfx matrix is given.");
    return  SiconosDOMTreeTools::getStringAttributeValue(MNode, "matrixPlugin");
  }

  /** return M matrix
   *   \return SimpleMatrix
   */
  inline const SimpleMatrix getMMatrix() const
  {
    if (isMPlugin())
      XMLException::selfThrow("FirstOrderNonLinearDSXML - getMMatrix : jacobianfx matrix is not given since M is calculated using a plug-in");
    return  SiconosDOMTreeTools::getSiconosMatrixValue(MNode);
  }

  /** to save the M matrix
   *   \param a SiconosMatrix
   */
  void setMMatrix(const SiconosMatrix&v);

  /** to save the Jacobianfx plugin of the LagrangianDSXML
   *   \param a string (name of the plug-in)
   */
  void setMPlugin(const std::string& plugin);

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
   *   \return SiconosVector
   */
  inline const SiconosVector getFVector() const
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

  // === Jacobianfx ===
  /** Return the name of the jacobianfx plug-in
   *  \return a string
   */
  inline const std::string getJacobianfxPlugin() const
  {
    if (!isJacobianfxPlugin())
      XMLException::selfThrow("FirstOrderNonLinearDSXML - getJacobianfxPlugin : jacobianfx is not calculated from a plugin since a jacobianfx matrix is given.");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianfxNode, "matrixPlugin");
  }

  /** return jacobianfx matrix
   *   \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianfxMatrix() const
  {
    if (isJacobianfxPlugin())
      XMLException::selfThrow("FirstOrderNonLinearDSXML - getJacobianfxMatrix : jacobianfx matrix is not given since jacobianfx is calculated using a plug-in");
    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianfxNode);
  }

  /** to save the jacobianfx matrix
   *   \param a SiconosMatrix
   */
  void setJacobianfxMatrix(const SiconosMatrix&v);

  /** to save the Jacobianfx plugin of the LagrangianDSXML
   *   \param a string (name of the plug-in)
   */
  void setJacobianfxPlugin(const std::string& plugin);

  /** Returns the xMemoryXML* of the DynamicalSystemXML
   *   \return SiconosMemoryXML*
   */
  inline SP::SiconosMemoryXML getXMemoryXML() const
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
    return (x0Node);
  }

  /** returns true if xNode is defined
   *  \return a bool
   */
  inline bool hasx() const
  {
    return (xNode);
  }

  /** returns true if xMemoryNode is defined
   *  \return a bool
   */
  inline bool hasXMemory()const
  {
    return (xMemoryNode);
  }

  /** returns true if MNode is defined
   *  \return a bool
   */
  inline bool hasM() const
  {
    return (MNode);
  }

  /** returns true if fNode is defined
   *  \return a bool
   */
  inline bool hasF() const
  {
    return (fNode);
  }

  /** returns true if jacobianfxNode is defined
   *  \return a bool
   */
  inline bool hasJacobianfx() const
  {
    return (jacobianfxNode);
  }

  /** Return true if M is calculated from a plugin
   *  \return a bool
   */
  inline bool isMPlugin() const
  {
    return xmlHasProp(MNode, (xmlChar *) MATRIXPLUGIN.c_str());
  }

  /** Return true if f is calculated from a plugin
   *  \return a bool
   */
  inline bool isFPlugin() const
  {
    return xmlHasProp(fNode, (xmlChar *) VECTORPLUGIN.c_str());
  }

  /** Return true if jacobianfx is calculated from a plugin
   *  \return a bool
   */
  inline bool isJacobianfxPlugin() const
  {
    return xmlHasProp(jacobianfxNode, (xmlChar *) MATRIXPLUGIN.c_str());
  }
};
DEFINE_SPTR(FirstOrderNonLinearDSXML);
#endif
