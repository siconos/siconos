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

/*! \file LagrangianDSXML.h

*/


#ifndef __LAGRANGIANNLDSXML__
#define __LAGRANGIANNLDSXML__

#include "DynamicalSystemXML.h"
#include "SiconosMemory.h"
#include "SimpleMatrix.h"
#include <vector>

class DynamicalSystemXML ;

const std::string LNLDS_Q = "q";
const std::string LNLDS_Q0 = "q0";
const std::string LNLDS_QMEMORY = "qMemory";
const std::string LNLDS_VELOCITY = "Velocity";
const std::string LNLDS_VELOCITY0 = "Velocity0";
const std::string LNLDS_VELOCITYMEMORY = "VelocityMemory";
const std::string LNLDS_QNLINERTIA = "NNL";
const std::string LNLDS_FINT = "FInt";
const std::string LNLDS_FEXT = "FExt";
const std::string LNLDS_Mass = "Mass";

/** XML management for LagrangianDS
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date 05/11/2004
 *
 */
class LagrangianDSXML : public DynamicalSystemXML
{
protected:

  xmlNodePtr qNode;
  xmlNodePtr q0Node;
  xmlNodePtr qMemoryNode;
  xmlNodePtr velocityNode;
  xmlNodePtr velocity0Node;
  xmlNodePtr velocityMemoryNode;
  xmlNodePtr MassNode;
  xmlNodePtr NNLNode;
  xmlNodePtr FIntNode;
  xmlNodePtr FExtNode;

  std::vector<xmlNodePtr> jacobianFIntNode;
  std::vector<xmlNodePtr> jacobianNNLNode;

  SiconosMemoryXMLSPtr qMemoryXML;
  SiconosMemoryXMLSPtr velocityMemoryXML;

  // === Constructors - Destructor ===
  LagrangianDSXML();

public:

  /** Build a LagrangianDSXML object from a DOM tree describing a DS
   *   \param xmlNodePtr LagrangianDSNode : the LagrangianDS DOM tree
   *   \param bool isBVP : if NonSmoothDynamicalSystem is BVP LagrangianDS have boundary condition
   */
  LagrangianDSXML(xmlNodePtr, bool);

  /** Destructor */
  virtual ~LagrangianDSXML();

  // Functions for members loading/setting

  // === q ===
  /** Return  q vector of the LagrangianDSXML
   *   \return SimpleVector : q vector of the LagrangianDSXML
   */
  inline const SimpleVector getQ() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(qNode);
  }

  /** allows to save the q of the LagrangianDSXML
   *   \param The q SiconosVector to save
   */
  inline void setQ(const SiconosVector &v)
  {
    if (!hasQ())
      qNode = SiconosDOMTreeTools::createVectorNode(rootNode, LNLDS_Q, v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(qNode, v);
  }

  // === q0 ===
  /** Return q0 vector of the LagrangianDSXML
   *   \return SimpleVector : q0 vector of the LagrangianDSXML
   */
  inline const SimpleVector getQ0() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(q0Node);
  }

  /** allows to save the q0 of the LagrangianDSXML
   *   \param The q0 SiconosVector to save
   */
  inline void  setQ0(const SiconosVector&v)
  {
    if (q0Node == NULL)
      q0Node = SiconosDOMTreeTools::createVectorNode(rootNode, LNLDS_Q0, v);
    else
      SiconosDOMTreeTools::setSiconosVectorNodeValue(q0Node, v);
  }

  // === qMemory ===
  /** Returns the qMemoryXML* of the DSXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXMLSPtr getQMemoryXML() const
  {
    return qMemoryXML;
  }

  /** allows to save the qMemory of the LagrangianDSXML
   *   \param SiconosMemory* : SiconosMemory to save
   */
  void setQMemory(const SiconosMemory& smem);

  /** Return the velocity of the LagrangianDSXML
   *   \return SimpleVector :  velocity vector of the LagrangianDSXML
   */

  // === Velocity ===
  inline  const SimpleVector getVelocity() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(velocityNode);
  }

  /** allows to save the velocity of the LagrangianDSXML
   *   \param The velocity SiconosVector to save
   */
  inline void setVelocity(const SiconosVector & v)
  {
    if (!hasVelocity())
      velocityNode = SiconosDOMTreeTools::createVectorNode(rootNode, LNLDS_VELOCITY, v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(velocityNode, v);
  }

  // === Velocity0 ===
  /** Return the initial velocity of the LagrangianDSXML
   *   \return SimpleVector : The velocity0 SiconosVector of the LagrangianDSXML
   */
  inline const SimpleVector getVelocity0() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(velocity0Node);
  }

  /** allows to save the velocity0 of the LagrangianDSXML
   *   \param The celocity0 SiconosVector to save
   */
  inline void setVelocity0(const SiconosVector&v)
  {
    if (velocity0Node == NULL)
    {
      velocity0Node = SiconosDOMTreeTools::createVectorNode(rootNode, LNLDS_VELOCITY0, v);
    }
    else
      SiconosDOMTreeTools::setSiconosVectorNodeValue(velocity0Node, v);
  }

  // === VelocityMemory ===
  /** Returns the velocityMemoryXML* of the DSXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXMLSPtr getVelocityMemoryXML() const
  {
    return velocityMemoryXML;
  }

  /** allows to save the velocityMemory of the LagrangianDSXML
   *   \param const SiconosMemory : SiconosMemory to save
   */
  void setVelocityMemory(const SiconosMemory& smem);

  // === NNL ===
  /** Return the NNL Plugin name of the LagrangianDSXML
   *   \return The NNL Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const std::string getNNLPlugin() const
  {
    if (!isNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getNNLPlugin : NNL is not calculated from a plugin ; NNL vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(NNLNode, VECTORPLUGIN);
  }

  /** Return the NNL vector of the LagrangianDSXML
   *   \return SimpleVector : NNL vector of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const SimpleVector getNNLVector() const
  {
    if (isNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getNNLVector : NNL vector is not given ; NNL is calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosVectorValue(NNLNode);
  }

  /** allows to save the NNL plugin of the LagrangianDSXML
   *   \param string : the string which contains the name and the location of the plugin
   */
  inline void setNNLPlugin(const std::string& plugin)
  {
    if (NNLNode == NULL)
    {
      NNLNode = SiconosDOMTreeTools::createSingleNode(rootNode, LNLDS_QNLINERTIA);
      xmlNewProp(NNLNode, (xmlChar*)(VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(NNLNode, VECTORPLUGIN, plugin);
  }

  /** allows to save the NNL vector of the LagrangianDSXML
   *   \return The NNL SiconosVector to save
   */
  inline void setNNLVector(const SiconosVector&v)
  {
    if (NNLNode == NULL)
    {
      NNLNode = SiconosDOMTreeTools::createVectorNode(rootNode, LNLDS_QNLINERTIA, v);
    }
    else
      SiconosDOMTreeTools::setSiconosVectorNodeValue(NNLNode, v);
  }

  // === FInt ===
  /** Return the FInt Plugin name of the LagrangianDSXML
   *   \return The FInt Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const std::string getFIntPlugin() const
  {
    if (!isFIntPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFIntPlugin : FInt is not calculated from a plugin ; FInt vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(FIntNode, VECTORPLUGIN);
  }

  /** Return the internal forces vector of the LagrangianDSXML
   *   \return SimpleVector : FInt SiconosVector of the LagrangianDSXML
   *   \exception XMLException
   */
  inline const SimpleVector getFIntVector() const
  {
    if (isFIntPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFIntVector : FInt vector is not given ; FInt is calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosVectorValue(FIntNode);
  }

  /** allows to save the FInt vector of the LagrangianDSXML
   *   \param The FInt SiconosVector to save
   */
  inline void setFIntVector(const SiconosVector&v)
  {
    if (!hasFInt())
    {
      FIntNode = SiconosDOMTreeTools::createVectorNode(rootNode, LNLDS_FINT, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(FIntNode, v);
  }

  /** allows to save the FExt plugin of the LagrangianDSXML
   *   \param string : the string which contains the name and the location of the plugin
   */
  inline void setFIntPlugin(const std::string& plugin)
  {
    if (FIntNode == NULL)
    {
      FIntNode = SiconosDOMTreeTools::createSingleNode(rootNode, LNLDS_FINT);
      xmlNewProp(FIntNode, (xmlChar*)(VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(FIntNode, VECTORPLUGIN, plugin);
  }

  // === FExt ===
  /** Return the FExt Plugin name of the LagrangianDSXML
   *   \return The FExt Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const std::string getFExtPlugin() const
  {
    if (!isFExtPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFExtPlugin : FExt is not calculated from a plugin ; FExt vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(FExtNode, VECTORPLUGIN);
  }

  /** Return the external forces vector of the LagrangianDSXML
   *   \return SimpleVector : FExt vector of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const SimpleVector getFExtVector() const
  {
    if (isFExtPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFExtVector : FExt matrix is not given ; FExt is calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosVectorValue(FExtNode);
  }

  /** allows to save the FInt vector of the LagrangianDSXML
   *   \param The FInt SiconosVector to save
   */
  inline void setFExtVector(const SiconosVector&v)
  {
    //SiconosDOMTreeTools::setSiconosVectorValue(FExtNode, v);
    SiconosDOMTreeTools::setSiconosVectorNodeValue(FExtNode, v);
  }

  /** allows to save the FExt plugin of the LagrangianDSXML
   *   \param string : the string which contains the name and the location of the plugin
   */
  inline void setFExtPlugin(const std::string& plugin)
  {
    if (FExtNode == NULL)
    {
      FExtNode = SiconosDOMTreeTools::createSingleNode(rootNode, LNLDS_FEXT);
      xmlNewProp(FExtNode, (xmlChar*)(VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else
    {
      SiconosDOMTreeTools::setStringAttributeValue(FExtNode, VECTORPLUGIN, plugin);
    }
  }

  /** Return the JacobianFInt Plug-in name
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return a string
   */
  inline const std::string getJacobianFIntPlugin(unsigned int i) const
  {
    if (!isJacobianFIntPlugin(i))
      XMLException::selfThrow("LagrangianDSXML - getJacobianFIntPlugin(i) : JacobianFInt is not calculated from a plugin ; JacobianFInt matrix is given. i=" + i);
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianFIntNode[i], MATRIXPLUGIN);
  }

  /** Return the JacobianFInt matrix
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return a SimpleMatrix
   */
  inline const SimpleMatrix getJacobianFIntMatrix(unsigned int i) const
  {
    if (isJacobianFIntPlugin(i))
      XMLException::selfThrow("LagrangianDSXML - getJacobianFIntMatrix(i) : JacobianFInt matrix is not given ; JacobianFInt is calculated from a plugin. i=" + i);

    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianFIntNode[i]);
  }

  /** To save the jacobianFInt plug-in
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param string : the string which contains the name and the location of the plugin
   */
  inline void setJacobianFIntPlugin(unsigned int i, const std::string& plugin)
  {
    std::string name = "Jacobian";
    if (i)
      name += "QFInt";
    else
      name += "VelocityFInt";
    if (jacobianFIntNode[i] == NULL)
    {
      jacobianFIntNode[i] = SiconosDOMTreeTools::createSingleNode(rootNode, name);
      xmlNewProp(jacobianFIntNode[i], (xmlChar*)(MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(jacobianFIntNode[i], name, plugin);
  }

  /**  To set the JacobianFInt matrix
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param The new value of JacobianFInt (SiconosMatrix)
   */
  inline void setJacobianFIntMatrix(unsigned int i, const SiconosMatrix& m)
  {
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianFIntNode[i], m);
  }

  /** Get the Jacobian NNL Plug-in name.
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return a string.
   */
  inline const std::string getJacobianNNLPlugin(unsigned int i) const
  {
    if (!isJacobianNNLPlugin(i))
      XMLException::selfThrow("LagrangianDSXML - getJacobianNNLPlugin : JacobianNNL is not calculated from a plugin ; JacobianNNL matrix is given. i=" + i);
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianNNLNode[i], MATRIXPLUGIN);
  }

  /** Get the Jacobian NNL matrix.
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return A SimpleMatrix.
   */
  inline const SimpleMatrix getJacobianNNLMatrix(unsigned int i) const
  {
    if (isJacobianNNLPlugin(i))
      XMLException::selfThrow("LagrangianDSXML - getJacobianNNLMatrix : JacobianNNL matrix is not given ; JacobianNNL is calculated from a plugin. i=" + i);

    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianNNLNode[i]);
  }

  /** To set the jacobian NNL plug-in name.
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param a string which contains the name and the location of the plug-in.
   */
  inline void setJacobianNNLPlugin(unsigned int i, const std::string& plugin)
  {
    std::string name = "Jacobian";
    if (i)
      name += "QNNL";
    else
      name += "VelocityNNL";
    if (jacobianNNLNode[i] == NULL)
    {
      jacobianNNLNode[i] = SiconosDOMTreeTools::createSingleNode(rootNode, name);
      xmlNewProp(jacobianNNLNode[i], (xmlChar*)(MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(jacobianNNLNode[i], name, plugin);
  }

  /** To set the Jacobian NNL matrix.
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param a SiconosMatrix.
   */
  inline void setJacobianNNLMatrix(unsigned int i, const SiconosMatrix& m)
  {
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianNNLNode[i], m);
  }

  /** Return the Mass Plug-in name of the LagrangianDSXML
   *   \return The Mass Plug-in name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const std::string getMassPlugin() const
  {
    if (!isMassPlugin())
      XMLException::selfThrow("LagrangianDSXML - getMPlugin : Mass is not calculated from a plugin ; Mass matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(MassNode, MATRIXPLUGIN);
  }

  /** Return the Mass matrix of the LagrangianDSXML
   *   \return The Mass SimpleMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const SimpleMatrix getMassMatrix() const
  {
    if (isMassPlugin())
      XMLException::selfThrow("LagrangianDSXML - getMMatrix : Mass matrix is not given ; Mass is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(MassNode);
  }

  /** To save the Mass plugin of the LagrangianDSXML
   *   \param string : the string which contains the name and the location of the plugin
   */
  inline void setMassPlugin(const std::string& plugin)
  {
    if (MassNode == NULL)
    {
      MassNode = SiconosDOMTreeTools::createSingleNode(rootNode, LNLDS_Mass);
      xmlNewProp(MassNode, (xmlChar*)(MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else
      SiconosDOMTreeTools::setStringAttributeValue(MassNode, MATRIXPLUGIN, plugin);
  }

  /** allows to save the M matrix of the LagrangianDSXML
   *   \return The M SiconosMatrix to save
   */
  inline void setMassMatrix(const SiconosMatrix& m)
  {
    if (MassNode == NULL)
    {
      MassNode = SiconosDOMTreeTools::createMatrixNode(rootNode, LNLDS_Mass, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(MassNode, m);
  }

  /** Return true if M is calculated from a plugin
   *  \return a bool
   *   \return True if M is calculated from plugin
   */
  inline const bool isMassPlugin() const
  {
    return xmlHasProp(MassNode, (xmlChar *) MATRIXPLUGIN.c_str());
  }

  /** Return true if M only given by a Matrix
   *  \return a bool
   */
  inline const bool isMMatrix() const
  {
    return !(xmlHasProp(MassNode, (xmlChar *) MATRIXPLUGIN.c_str()));
  }

  /** Return true if QLNInertia is calculated from a plugin
   *  \return a bool
   */
  inline const bool isNNLPlugin() const
  {
    return xmlHasProp(NNLNode, (xmlChar *) VECTORPLUGIN.c_str());
  }

  /** Return true if FInt is calculated from a plugin
   *  \return a bool
   */
  inline const bool isFIntPlugin() const
  {
    return xmlHasProp(FIntNode, (xmlChar *) VECTORPLUGIN.c_str());
  }

  /** Return true if FExt is calculated from a plugin
   *  \return a bool
   */
  inline const bool isFExtPlugin() const
  {
    return xmlHasProp(FExtNode, (xmlChar *) VECTORPLUGIN.c_str());
  }

  /** Return true if JacobianFInt is calculated from a plugin
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return a bool
   */
  inline const bool isJacobianFIntPlugin(unsigned int i) const
  {
    return xmlHasProp(jacobianFIntNode[i], (xmlChar *) MATRIXPLUGIN.c_str());
  }

  /** Return true if Jacobian NNL is calculated from a plugin
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return a bool
   */
  inline const bool isJacobianNNLPlugin(unsigned int i) const
  {
    return xmlHasProp(jacobianNNLNode[i], (xmlChar *) MATRIXPLUGIN.c_str());
  }

  /** determines if Mass is defined in the DOM tree
   *  \return a bool
   */
  inline const bool hasQ0() const
  {
    return (q0Node != NULL);
  }

  /** determines if Mass is defined in the DOM tree
   *  \return a bool
   */
  inline const bool hasVelocity0() const
  {
    return (velocity0Node != NULL);
  }

  /** determines if Mass is defined in the DOM tree
   *  \return a bool
   */
  inline const bool hasMass() const
  {
    return (MassNode != NULL);
  }

  /** determines if FInt is defined in the DOM tree
   *  \return a bool
   */
  inline const bool hasFInt() const
  {
    return (FIntNode != NULL);
  }

  /** determines if FExt is defined in the DOM tree
   *  \return a bool
   */
  inline const bool hasFExt() const
  {
    return (FExtNode != NULL);
  }

  /** determines if jacobianFInt is defined in the DOM tree
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return a bool
   */
  inline const bool hasJacobianFInt(unsigned int i) const
  {
    return (jacobianFIntNode[i] != NULL);
  }


  /** determines if jacobianNNL is defined in the DOM tree
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return a bool
   */
  inline const bool hasJacobianNNL(unsigned int i) const
  {
    return (jacobianNNLNode[i] != NULL);
  }

  /** determines if NNL is defined in the DOM tree
   *  \return a bool
   */
  inline const bool hasNNL() const
  {
    return (NNLNode != NULL);
  }

  /** returns true if qMemoryNode is defined
   *  \return a bool
   */
  inline const bool hasQMemory() const
  {
    return (qMemoryNode != NULL);
  }

  /** returns true if velocityMemoryNode is defined
   *  \return a bool
   */
  inline const bool hasVelocityMemory() const
  {
    return (velocityMemoryNode != NULL);
  }

  /** returns true if qNode is defined
   *  \return a bool
   */
  inline const bool hasQ() const
  {
    return (qNode != NULL);
  }

  /** returns true if velocityNode is defined
   *  \return a bool
   */
  inline const bool hasVelocity() const
  {
    return (velocityNode != NULL);
  }
};

#endif
