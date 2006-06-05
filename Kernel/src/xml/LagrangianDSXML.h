/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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

/** \class LagrangianDSXML
 *   \brief This class manages Lagrangian NLDS data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.2.0.
 *   \date 05/11/2004
 *
 *
 *
 * LagrangianDSXML allows to manage data of a LagrangianDS DOM tree.
 */

#ifndef __LAGRANGIANNLDSXML__
#define __LAGRANGIANNLDSXML__

#include "DynamicalSystemXML.h"

class DynamicalSystemXML ;

const std::string LNLDS_Q = "q";
const std::string LNLDS_Q0 = "q0";
const std::string LNLDS_QMEMORY = "qMemory";

const std::string LNLDS_VELOCITY = "Velocity";
const std::string LNLDS_VELOCITY0 = "Velocity0";
const std::string LNLDS_VELOCITYMEMORY = "VelocityMemory";

const std::string LNLDS_QNLINERTIA = "NNL";
const std::string LNLDS_FINT = "Fint";
const std::string LNLDS_FEXT = "Fext";

const std::string LNLDS_JACOBIANQFINT = "JacobianQFint";
const std::string LNLDS_JACOBIANVELOCITYFINT = "JacobianVelocityFint";
const std::string LNLDS_JACOBIANQQNLINERTIA = "JacobianQNNL";
const std::string LNLDS_JACOBIANVELOCITYQNLINERTIA = "JacobianVelocityNNL";

const std::string LNLDS_M = "M";
const std::string LNLDS_NDOF = "ndof";
const std::string LNLDS_MATRIXPLUGIN = "matrixPlugin";
const std::string LNLDS_VECTORPLUGIN = "vectorPlugin";

#include "check.h"

class LagrangianDSXML : public DynamicalSystemXML
{
protected:

  xmlNodePtr ndofNode;
  xmlNodePtr qNode;
  xmlNodePtr q0Node;
  xmlNodePtr qMemoryNode;

  xmlNodePtr velocityNode;
  xmlNodePtr velocity0Node;
  xmlNodePtr velocityMemoryNode;

  xmlNodePtr MNode;
  xmlNodePtr NNLNode;
  xmlNodePtr FintNode;
  xmlNodePtr FextNode;

  xmlNodePtr jacobianQFintNode;
  xmlNodePtr jacobianVelocityFintNode;
  xmlNodePtr jacobianQNNLNode;
  xmlNodePtr jacobianVelocityNNLNode;

  SiconosMemoryXML * qMemoryXML;
  SiconosMemoryXML * velocityMemoryXML;

public:

  // === Constructors - Destructor ===
  LagrangianDSXML();

  virtual ~LagrangianDSXML();

  /** \fn LagrangianDSXML(xmlNodePtr LagrangianDSNode, const bool isBVP)
   *   \brief Build a LagrangianDSXML object from a DOM tree describing a DS
   *   \param xmlNodePtr LagrangianDSNode : the LagrangianDS DOM tree
   *   \param bool isBVP : if NonSmoothDynamicalSystem is BVP LagrangianDS have boundary condition
   */
  LagrangianDSXML(xmlNodePtr, const bool&);

  // Functions for members loading/setting

  // === q ===
  /** \fn const SimpleVector getQ() const
   *   \brief Return  q vector of the LagrangianDSXML
   *   \return SimpleVector : q vector of the LagrangianDSXML
   */
  inline const SimpleVector getQ() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(qNode);
  }

  /** \fn void setQ(const SiconosVector &v)
   *   \brief allows to save the q of the LagrangianDSXML
   *   \param The q SiconosVector to save
   */
  inline void setQ(const SiconosVector &v)
  {
    if (!hasQ())
      qNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LNLDS_Q, v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(qNode, v);
  }

  // === q0 ===
  /** \fn const SimpleVector getQ()const
   *   \brief Return q0 vector of the LagrangianDSXML
   *   \return SimpleVector : q0 vector of the LagrangianDSXML
   */
  inline const SimpleVector getQ0() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(q0Node);
  }

  /** \fn void  setQ0(const SiconosVector &v)
   *   \brief allows to save the q0 of the LagrangianDSXML
   *   \param The q0 SiconosVector to save
   */
  inline void  setQ0(const SiconosVector&v)
  {
    if (q0Node == NULL)
      q0Node = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LNLDS_Q0, v);
    else
      SiconosDOMTreeTools::setSiconosVectorNodeValue(q0Node, v);
  }

  // === qMemory ===
  /** \fn SiconosMemoryXML* getQMemoryXML() const
   *   \brief Returns the qMemoryXML* of the DSXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXML* getQMemoryXML() const
  {
    return qMemoryXML;
  }

  /** \fn void setQMemory( const SiconosMemory& smem)
   *   \brief allows to save the qMemory of the LagrangianDSXML
   *   \param SiconosMemory* : SiconosMemory to save
   */
  void setQMemory(const SiconosMemory& smem);

  /** \fn  const SimpleVector getVelocity() const
   *   \brief Return the velocity of the LagrangianDSXML
   *   \return SimpleVector :  velocity vector of the LagrangianDSXML
   */

  // === Velocity ===
  inline  const SimpleVector getVelocity() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(velocityNode);
  }

  /** \fn void setVelocity( const SiconosVector & v)
   *   \brief allows to save the velocity of the LagrangianDSXML
   *   \param The velocity SiconosVector to save
   */
  inline void setVelocity(const SiconosVector & v)
  {
    if (!hasVelocity())
      velocityNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LNLDS_VELOCITY, v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(velocityNode, v);
  }

  // === Velocity0 ===
  /** \fn const SimpleVector getVelocity0() const
   *   \brief Return the initial velocity of the LagrangianDSXML
   *   \return SimpleVector : The velocity0 SiconosVector of the LagrangianDSXML
   */
  inline const SimpleVector getVelocity0() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(velocity0Node);
  }

  /** \fn void setVelocity0(const SiconosVector&v)
   *   \brief allows to save the velocity0 of the LagrangianDSXML
   *   \param The celocity0 SiconosVector to save
   */
  inline void setVelocity0(const SiconosVector&v)
  {
    if (velocity0Node == NULL)
    {
      velocity0Node = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LNLDS_VELOCITY0, v);
    }
    else
      SiconosDOMTreeTools::setSiconosVectorNodeValue(velocity0Node, v);
  }

  // === VelocityMemory ===
  /** \fn SiconosMemoryXML* getVelocityMemoryXML() const
   *   \brief Returns the velocityMemoryXML* of the DSXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXML* getVelocityMemoryXML() const
  {
    return velocityMemoryXML;
  }

  /** \fn void setVelocityMemory(const SiconosMemory& smem)
   *   \brief allows to save the velocityMemory of the LagrangianDSXML
   *   \param const SiconosMemory : SiconosMemory to save
   */
  void setVelocityMemory(const SiconosMemory& smem);

  // === NNL ===
  /** \fn inline const string getNNLPlugin() const
   *   \brief Return the NNL Plugin name of the LagrangianDSXML
   *   \return The NNL Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const std::string getNNLPlugin() const
  {
    if (!isNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getNNLPlugin : NNL is not calculated from a plugin ; NNL vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(NNLNode, LNLDS_VECTORPLUGIN);
  }

  /** \fn const SimpleVector getNNLVector() const
   *   \brief Return the NNL vector of the LagrangianDSXML
   *   \return SimpleVector : NNL vector of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const SimpleVector getNNLVector() const
  {
    if (isNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getNNLVector : NNL vector is not given ; NNL is calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosVectorValue(NNLNode);
  }

  /** \fn void setNNLPlugin(const string& plugin)
   *   \brief allows to save the NNL plugin of the LagrangianDSXML
   *   \param string : the string which contains the name and the location of the plugin
   */
  inline void setNNLPlugin(const std::string& plugin)
  {
    if (NNLNode == NULL)
    {
      NNLNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_QNLINERTIA);
      xmlNewProp(NNLNode, (xmlChar*)(LNLDS_VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(NNLNode, LNLDS_VECTORPLUGIN, plugin);
  }

  /** \fn void setNNLVector(const SiconosVector&v)
   *   \brief allows to save the NNL vector of the LagrangianDSXML
   *   \return The NNL SiconosVector to save
   */
  inline void setNNLVector(const SiconosVector&v)
  {
    if (NNLNode == NULL)
    {
      NNLNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LNLDS_QNLINERTIA, v);
    }
    else
      SiconosDOMTreeTools::setSiconosVectorNodeValue(NNLNode, v);
  }

  // === FInt ===
  /** \fn inline const string getFintPlugin() const
   *   \brief Return the Fint Plugin name of the LagrangianDSXML
   *   \return The Fint Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const std::string getFintPlugin() const
  {
    if (!isFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFintPlugin : Fint is not calculated from a plugin ; Fint vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(FintNode, LNLDS_VECTORPLUGIN);
  }

  /** \fn const SimpleVector getFintVector() const
   *   \brief Return the internal forces vector of the LagrangianDSXML
   *   \return SimpleVector : Fint SiconosVector of the LagrangianDSXML
   *   \exception XMLException
   */
  inline const SimpleVector getFintVector() const
  {
    if (isFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFintVector : Fint vector is not given ; Fint is calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosVectorValue(FintNode);
  }

  /** \fn void setFintVector(const SiconosVector&v)
   *   \brief allows to save the Fint vector of the LagrangianDSXML
   *   \param The Fint SiconosVector to save
   */
  inline void setFintVector(const SiconosVector&v)
  {
    if (!hasFint())
    {
      FintNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LNLDS_FINT, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(FintNode, v);
  }

  /** \fn void setFextPlugin(const string& plugin)
   *   \brief allows to save the Fext plugin of the LagrangianDSXML
   *   \param string : the string which contains the name and the location of the plugin
   */
  inline void setFintPlugin(const std::string& plugin)
  {
    if (FintNode == NULL)
    {
      FintNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_FINT);
      xmlNewProp(FintNode, (xmlChar*)(LNLDS_VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(FintNode, LNLDS_VECTORPLUGIN, plugin);
  }

  // === Fext ===
  /** \fn inline const string getFextPlugin() const
   *   \brief Return the Fext Plugin name of the LagrangianDSXML
   *   \return The Fext Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const std::string getFextPlugin() const
  {
    if (!isFextPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFextPlugin : Fext is not calculated from a plugin ; Fext vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(FextNode, LNLDS_VECTORPLUGIN);
  }

  /** \fn const SimpleVector getFextVector() const
   *   \brief Return the external forces vector of the LagrangianDSXML
   *   \return SimpleVector : Fext vector of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const SimpleVector getFextVector() const
  {
    if (isFextPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFextVector : Fext matrix is not given ; Fext is calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosVectorValue(FextNode);
  }

  /** \fn void setFextVector(const SiconosVector&v)
   *   \brief allows to save the Fint vector of the LagrangianDSXML
   *   \param The Fint SiconosVector to save
   */
  inline void setFextVector(const SiconosVector&v)
  {
    //SiconosDOMTreeTools::setSiconosVectorValue(FextNode, v);
    SiconosDOMTreeTools::setSiconosVectorNodeValue(FextNode, v);
  }

  /** \fn void setFextPlugin(const string& plugin)
   *   \brief allows to save the Fext plugin of the LagrangianDSXML
   *   \param string : the string which contains the name and the location of the plugin
   */
  inline void setFextPlugin(const std::string& plugin)
  {
    if (FextNode == NULL)
    {
      FextNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_FEXT);
      xmlNewProp(FextNode, (xmlChar*)(LNLDS_VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else
    {
      SiconosDOMTreeTools::setStringAttributeValue(FextNode, LNLDS_VECTORPLUGIN, plugin);
    }
  }

  /** \fn inline const string getJacobianQFintPlugin() const
   *   \brief Return the JacobianQFint Plugin name of the LagrangianDSXML
   *   \return The JacobianQFint Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const std::string getJacobianQFintPlugin() const
  {
    if (!isJacobianQFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianQFintPlugin : JacobianQFint is not calculated from a plugin ; JacobianQFint matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianQFintNode, LNLDS_MATRIXPLUGIN);
  }

  /** \fn const SimpleMatrix getJacobianQFintMatrix() const
   *   \brief Return the JacobianQFint matrix of the LagrangianDSXML
   *   \return The JacobianQFint SimpleMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const SimpleMatrix getJacobianQFintMatrix() const
  {
    if (isJacobianQFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianQFintMatrix : JacobianQFint matrix is not given ; JacobianQFint is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianQFintNode);
  }

  /** \fn void setJacobianQFintPlugin(const string& plugin)
   *   \brief allows to save the jacobianQFint plugin of the LagrangianDSXML
   *   \param string : the string which contains the name and the location of the plugin
   */
  inline void setJacobianQFintPlugin(const std::string& plugin)
  {
    if (jacobianQFintNode == NULL)
    {
      jacobianQFintNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_JACOBIANQFINT);
      xmlNewProp(jacobianQFintNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(jacobianQFintNode, LNLDS_JACOBIANQFINT, plugin);
  }

  /** \fn void setJacobianQFintMatrix(const SiconosMatrix& m)
   *   \brief allows to save the JacobianQFint matrix of the LagrangianDSXML
   *   \return The JacobianQFint SiconosMatrix to save
   */
  inline void setJacobianQFintMatrix(const SiconosMatrix& m)
  {
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianQFintNode, m);
  }

  /** \fn inline const string getJacobianVelocityFintPlugin() const
   *   \brief Return the JacobianVelocityFint Plugin name of the LagrangianDSXML
   *   \return The JacobianVelocityFint Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const std::string getJacobianVelocityFintPlugin() const
  {
    if (!isJacobianVelocityFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityFintPlugin : JacobianVelocityFint is not calculated from a plugin ; JacobianVelocityFint matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianVelocityFintNode, LNLDS_MATRIXPLUGIN);
  }

  /** \fn const SimpleMatrix getJacobianVelocityFintMatrix() const
   *   \brief Return the JacobianVelocityFint matrix of the LagrangianDSXML
   *   \return The JacobianVelocityFint SimpleMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const SimpleMatrix getJacobianVelocityFintMatrix() const
  {
    if (isJacobianVelocityFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityFintMatrix : JacobianVelocityFint matrix is not given ; JacobianVelocityFint is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianVelocityFintNode);
  }

  /** \fn void setJacobianVelocityFintPlugin(const string& plugin)
   *   \brief allows to save the jacobianVelocityFint plugin of the LagrangianDSXML
   *   \param string : the string which contains the name and the location of the plugin
   */
  inline void setJacobianVelocityFintPlugin(const std::string& plugin)
  {
    if (jacobianVelocityFintNode == NULL)
    {
      jacobianVelocityFintNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_JACOBIANVELOCITYFINT);
      xmlNewProp(jacobianVelocityFintNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(jacobianVelocityFintNode, LNLDS_JACOBIANVELOCITYFINT, plugin);
  }

  /** \fn void setJacobianVelocityFintMatrix(const SiconosMatrix& m)
   *   \brief allows to save the JacobianVelocityFint matrix of the LagrangianDSXML
   *   \return The JacobianVelocityFint SiconosMatrix to save
   */
  inline void setJacobianVelocityFintMatrix(const SiconosMatrix& m)
  {
    //SiconosDOMTreeTools::setSiconosMatrixValue(jacobianVelocityFintNode, m);
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianVelocityFintNode, m);
  }

  /** \fn inline const string getJacobianQQPlugin() const
   *   \brief Return the JacobianQQ Plugin name of the LagrangianDSXML
   *   \return The JacobianQQ Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const std::string getJacobianQNNLPlugin() const
  {
    if (!isJacobianQNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianQNNLPlugin : JacobianQNNL is not calculated from a plugin ; JacobianQNNL matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianQNNLNode, LNLDS_MATRIXPLUGIN);
  }

  /** \fn const SimpleMatrix getJacobianQQMatrix() const
   *   \brief Return the JacobianQQ matrix of the LagrangianDSXML
   *   \return The JacobianQQ SimpleMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const SimpleMatrix getJacobianQNNLMatrix() const
  {
    if (isJacobianQNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianQNNLMatrix : JacobianQNNL matrix is not given ; JacobianQNNL is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianQNNLNode);
  }

  /** \fn void setJacobianQNNLPlugin(const string& plugin)
   *   \brief allows to save the jacobianQNNL plugin of the LagrangianDSXML
   *   \param string : the string which contains the name and the location of the plugin
   */
  inline void setJacobianQNNLPlugin(const std::string& plugin)
  {
    if (jacobianQNNLNode == NULL)
    {
      jacobianQNNLNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_JACOBIANQQNLINERTIA);
      xmlNewProp(jacobianQNNLNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(jacobianQNNLNode, LNLDS_JACOBIANQQNLINERTIA, plugin);
  }

  /** \fn void setJacobianQQMatrix(const SiconosMatrix& m)
   *   \brief allows to save the JacobianQQ matrix of the LagrangianDSXML
   *   \return The JacobianQQ SiconosMatrix to save
   */
  inline void setJacobianQNNLMatrix(const SiconosMatrix& m)
  {
    //SiconosDOMTreeTools::setSiconosMatrixValue(jacobianQNNLNode, m);
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianQNNLNode, m);
  }

  /** \fn inline const string getJacobianVelocityNNLPlugin() const
   *   \brief Return the JacobianVelocityNNL Plugin name of the LagrangianDSXML
   *   \return The JacobianVelocityNNL Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const std::string getJacobianVelocityNNLPlugin() const
  {
    if (!isJacobianVelocityNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityNNLPlugin : JacobianVelocityNNL is not calculated from a plugin ; JacobianVelocityNNL matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianVelocityNNLNode, LNLDS_MATRIXPLUGIN);
  }

  /** \fn const SimpleMatrix getJacobianVelocityNNLMatrix() const
   *   \brief Return the JacobianVelocityNNL matrix of the LagrangianDSXML
   *   \return The JacobianVelocityNNL SimpleMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const SimpleMatrix getJacobianVelocityNNLMatrix() const
  {
    if (isJacobianVelocityNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityNNLMatrix : JacobianVelocityNNL matrix is not given ; JacobianVelocityNNL is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianVelocityNNLNode);
  }

  /** \fn void setJacobianVelocityNNLPlugin(const string& plugin)
   *   \brief allows to save the jacobianVelocityNNLPlugin plugin of the LagrangianDSXML
   *   \param string : the string which contains the name and the location of the plugin
   */
  inline void setJacobianVelocityNNLPlugin(const std::string& plugin)
  {
    if (jacobianVelocityNNLNode == NULL)
    {
      jacobianVelocityNNLNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_JACOBIANVELOCITYQNLINERTIA);
      xmlNewProp(jacobianVelocityNNLNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(jacobianVelocityNNLNode, LNLDS_JACOBIANVELOCITYQNLINERTIA, plugin);
  }

  /** \fn void setJacobianVelocityNNLMatrix(const SiconosMatrix& m)
   *   \brief allows to save the JacobianVelocityNNL matrix of the LagrangianDSXML
   *   \return The JacobianVelocityNNL SiconosMatrix to save
   */
  inline void setJacobianVelocityNNLMatrix(const SiconosMatrix& m)
  {
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianVelocityNNLNode, m);
  }

  /** \fn inline const string getMPlugin() const
   *   \brief Return the M Plugin name of the LagrangianDSXML
   *   \return The M Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const std::string getMPlugin() const
  {
    if (!isMPlugin())
      XMLException::selfThrow("LagrangianDSXML - getMPlugin : M is not calculated from a plugin ; M matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(MNode, LNLDS_MATRIXPLUGIN);
  }

  /** \fn const SimpleMatrix getMMatrix() const
   *   \brief Return the M matrix of the LagrangianDSXML
   *   \return The M SimpleMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline const SimpleMatrix getMMatrix() const
  {
    if (isMPlugin())
      XMLException::selfThrow("LagrangianDSXML - getMMatrix : M matrix is not given ; M is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(MNode);
  }

  /** \fn void setMPlugin(const string& plugin)
   *   \brief allows to save the jacobianVelocityNNLPlugin plugin of the LagrangianDSXML
   *   \param string : the string which contains the name and the location of the plugin
   */
  inline void setMPlugin(const std::string& plugin)
  {
    if (MNode == NULL)
    {
      MNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_M);
      xmlNewProp(MNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else
      SiconosDOMTreeTools::setStringAttributeValue(MNode, LNLDS_MATRIXPLUGIN, plugin);
  }

  /** \fn void setMMatrix(const SiconosMatrix& m)
   *   \brief allows to save the M matrix of the LagrangianDSXML
   *   \return The M SiconosMatrix to save
   */
  inline void setMMatrix(const SiconosMatrix& m)
  {
    if (MNode == NULL)
    {
      MNode = SiconosDOMTreeTools::createMatrixNode(rootDynamicalSystemXMLNode, LNLDS_M, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(MNode, m);
  }


  /** \fn const unsigned int getNdof() const
   *   \brief Return the ndof for the LagrangianDSXML
   *   \return The ndof integer for the LagrangianDSXML
   */
  inline const unsigned int getNdof() const
  {
    return  SiconosDOMTreeTools::getContentValue<unsigned int>(ndofNode);
  }

  /** \fn void setNdof(const unsigned int& i)
   *   \brief to save the ndof for the LagrangianDSXML
   *   \return The ndof integer to save
   */
  inline void setNdof(const unsigned int& i)
  {
    if (ndofNode == NULL)
      ndofNode = SiconosDOMTreeTools::createIntegerNode(rootDynamicalSystemXMLNode, LNLDS_NDOF, i);
    else SiconosDOMTreeTools::setIntegerContentValue(ndofNode, i);
  }


  /** \fn const bool isMPlugin() const
   *   \brief Return true if M is calculated from a plugin
   *   \return True if M is calculated from plugin
   */
  inline const bool isMPlugin() const
  {
    return xmlHasProp((xmlNodePtr)MNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str());
  }

  /** \fn const bool isMMatrix() const
   *   \brief Return true if M only given by a Matrix
   *   \return True if M is given by a Matrix
   */
  inline const bool isMMatrix() const
  {
    return !(xmlHasProp((xmlNodePtr)MNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str()));
  }

  /** \fn const bool isQLNInertiaPlugin() const
   *   \brief Return true if QLNInertia is calculated from a plugin
   *   \return True if QLNInertia is calculated from plugin
   */
  inline const bool isNNLPlugin() const
  {
    return xmlHasProp((xmlNodePtr)NNLNode, (xmlChar *) LNLDS_VECTORPLUGIN.c_str());
  }

  /** \fn const bool isFintPlugin() const
   *   \brief Return true if Fint is calculated from a plugin
   *   \return True if Fint is calculated from plugin
   */
  inline const bool isFintPlugin() const
  {
    return xmlHasProp((xmlNodePtr)FintNode, (xmlChar *) LNLDS_VECTORPLUGIN.c_str());
  }

  /** \fn const bool isFextPlugin() const
   *   \brief Return true if Fext is calculated from a plugin
   *   \return True if Fext is calculated from plugin
   */
  inline const bool isFextPlugin() const
  {
    return xmlHasProp((xmlNodePtr)FextNode, (xmlChar *) LNLDS_VECTORPLUGIN.c_str());
  }

  /** \fn const bool isJacobianQFintPlugin() const
   *   \brief Return true if JacobianQFint is calculated from a plugin
   *   \return True if JacobianQFint is calculated from plugin
   */
  inline const bool isJacobianQFintPlugin() const
  {
    return xmlHasProp((xmlNodePtr)jacobianQFintNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str());
  }

  /** \fn const bool isJacobianVelocityFintPlugin() const
   *   \brief Return true if JacobianVelocityFint is calculated from a plugin
   *   \return True if JacobianVelocityFint is calculated from plugin
   */
  inline const bool isJacobianVelocityFintPlugin() const
  {
    return xmlHasProp((xmlNodePtr)jacobianVelocityFintNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str());
  }

  /** \fn const bool isJacobianQQPlugin() const
   *   \brief Return true if JacobianQQ is calculated from a plugin
   *   \return True if JacobianQQ is calculated from plugin
   */
  inline const bool isJacobianQNNLPlugin() const
  {
    return xmlHasProp((xmlNodePtr)jacobianQNNLNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str());
  }

  /** \fn const bool isJacobianVelocityNNLPlugin() const
   *   \brief Return true if JacobianVelocityNNL is calculated from a plugin
   *   \return True if JacobianVelocityNNL is calculated from plugin
   */
  inline const bool isJacobianVelocityNNLPlugin() const
  {
    return xmlHasProp((xmlNodePtr)jacobianVelocityNNLNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str());
  }


  /** \fn const bool hasMass() const
   *  \brief determines if Mass is defined in the DOM tree
   *  \return bool : true if Mass is defined, false otherwise
   */
  inline const bool hasMass() const
  {
    return (MNode != NULL);
  }

  /** \fn const bool hasFint() const
   *  \brief determines if Fint is defined in the DOM tree
   *  \return const bool : true if Fint is defined, false otherwise
   */
  inline const bool hasFint() const
  {
    return (FintNode != NULL);
  }

  /** \fn const bool hasFext() const
   *  \brief determines if Fext is defined in the DOM tree
   *  \return const bool : true if Fext is defined, false otherwise
   */
  inline const bool hasFext() const
  {
    return (FextNode != NULL);
  }

  /** \fn const bool hasJacobianQFint() const
   *  \brief determines if jacobianQFint is defined in the DOM tree
   *  \return const bool : true if jacobianQFint is defined, false otherwise
   */
  inline const bool hasJacobianQFint() const
  {
    return (jacobianQFintNode != NULL);
  }

  /** \fn const bool hasJacobianVelocityFint() const
   *  \brief determines if jacobianVelocityFint is defined in the DOM tree
   *  \return const bool : true if jacobianVelocityFint is defined, false otherwise
   */
  inline const bool hasJacobianVelocityFint() const
  {
    return (jacobianVelocityFintNode != NULL);
  }

  /** \fn const bool hasJacobianQNNL() const
   *  \brief determines if jacobianQNNL is defined in the DOM tree
   *  \return const bool : true if jacobianQNNL is defined, false otherwise
   */
  inline const bool hasJacobianQNNL() const
  {
    return (jacobianQNNLNode != NULL);
  }

  /** \fn const bool hasJacobianVelocityNNL() const
   *  \brief determines if jacobianVelocityNNL is defined in the DOM tree
   *  \return const bool : true if jacobianVelocityNNL is defined, false otherwise
   */
  inline const bool hasJacobianVelocityNNL() const
  {
    return (jacobianVelocityNNLNode != NULL);
  }

  /** \fn const bool hasNNL() const
   *  \brief determines if NNL is defined in the DOM tree
   *  \return const bool : true if NNL is defined, false otherwise
   */
  inline const bool hasNNL() const
  {
    return (NNLNode != NULL);
  }

  /** \fn const bool hasQMemory() const
   *  \brief returns true if qMemoryNode is defined
   *  \return true if qMemoryNode is defined
   */
  inline const bool hasQMemory() const
  {
    return (qMemoryNode != NULL);
  }

  /** \fn const bool hasVelocityMemory() const
   *  \brief returns true if velocityMemoryNode is defined
   *  \return true if velocityMemoryNode is defined
   */
  inline const bool hasVelocityMemory() const
  {
    return (velocityMemoryNode != NULL);
  }

  /** \fn const bool hasQ() const
   *  \brief returns true if qNode is defined
   *  \return true if qNode is defined
   */
  inline const bool hasQ() const
  {
    return (qNode != NULL);
  }

  /** \fn const bool hasVelocity() const
   *  \brief returns true if velocityNode is defined
   *  \return true if velocityNode is defined
   */
  inline const bool hasVelocity() const
  {
    return (velocityNode != NULL);
  }


  /** \fn void updateDynamicalSystemXML(xmlNodePtr, DynamicalSystem*, BoundaryCondition*)
   *   \brief makes the operations to add a DynamicalSystem to the NonSmoothDynamicalSystemXML
   *   \param xmlNodePtr : the root node of this DynamicalSystem
   *   \param DynamicalSystem* : the DynamicalSystem of this DynamicalSystemXML
   *   \param BoundaryCondition* : the BoundaryCondition of the DS if the NonSmoothDynamicalSystem is BVP (optional)
   */
  void updateDynamicalSystemXML(xmlNodePtr, DynamicalSystem*, BoundaryCondition* bc = NULL);


};

#endif
