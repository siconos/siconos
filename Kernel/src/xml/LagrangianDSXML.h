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

/** \class LagrangianDSXML
 *   \brief This class manages Lagrangian NLDS data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.1.4.
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
public:

  // === Constructors - Destructor ===
  LagrangianDSXML();

  virtual ~LagrangianDSXML();

  /** \fn LagrangianDSXML(xmlNode * LagrangianDSNode, bool isBVP)
   *   \brief Build a LagrangianDSXML object from a DOM tree describing a DS
   *   \param xmlNode * LagrangianDSNode : the LagrangianDS DOM tree
   *   \param bool isBVP : if NonSmoothDynamicalSystem is BVP LagrangianDS have boundary condition
   */
  LagrangianDSXML(xmlNode * LagrangianDSNode, bool isBVP);

  // Functions for members loading/setting

  // === q ===
  /** \fn SimpleVector getQ()
   *   \brief Return  q vector of the LagrangianDSXML
   *   \return SimpleVector : q vector of the LagrangianDSXML
   */
  inline SimpleVector getQ()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(qNode);
  }

  /** \fn void setQ(SiconosVector *v)
   *   \brief allows to save the q of the LagrangianDSXML
   *   \param The q SiconosVector to save
   */
  inline void setQ(SiconosVector *v)
  {
    if (!hasQ())
      qNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LNLDS_Q, *v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(qNode, *v);
  }

  // === q0 ===
  /** \fn SimpleVector getQ()
   *   \brief Return q0 vector of the LagrangianDSXML
   *   \return SimpleVector : q0 vector of the LagrangianDSXML
   */
  inline SimpleVector getQ0()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(q0Node);
  }

  /** \fn void  setQ0(SiconosVector *v)
   *   \brief allows to save the q0 of the LagrangianDSXML
   *   \param The q0 SiconosVector to save
   */
  inline void  setQ0(SiconosVector *v)
  {
    if (q0Node == NULL)
      q0Node = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LNLDS_Q0, *v);
    else
      SiconosDOMTreeTools::setSiconosVectorNodeValue(q0Node, *v);
  }

  // === qMemory ===
  /** \fn SiconosMemoryXML* getQMemoryXML()
   *   \brief Returns the qMemoryXML* of the DSXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXML* getQMemoryXML()
  {
    return qMemoryXML;
  }

  /** \fn void setQMemory(SiconosMemory* smem)
   *   \brief allows to save the qMemory of the LagrangianDSXML
   *   \param SiconosMemory* : SiconosMemory to save
   */
  void setQMemory(SiconosMemory* smem);

  /** \fn SimpleVector getVelocity()
   *   \brief Return the velocity of the LagrangianDSXML
   *   \return SimpleVector :  velocity vector of the LagrangianDSXML
   */

  // === Velocity ===
  inline SimpleVector getVelocity()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(velocityNode);
  }

  /** \fn void setVelocity(SiconosVector *v)
   *   \brief allows to save the velocity of the LagrangianDSXML
   *   \param The velocity SiconosVector to save
   */
  inline void setVelocity(SiconosVector *v)
  {
    if (!hasVelocity())
      velocityNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LNLDS_VELOCITY, *v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(velocityNode, *v);
  }

  // === Velocity0 ===
  /** \fn SimpleVector getVelocity0()
   *   \brief Return the initial velocity of the LagrangianDSXML
   *   \return SimpleVector : The velocity0 SiconosVector of the LagrangianDSXML
   */
  inline SimpleVector getVelocity0()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(velocity0Node);
  }

  /** \fn void setVelocity0(SiconosVector *v)
   *   \brief allows to save the velocity0 of the LagrangianDSXML
   *   \param The celocity0 SiconosVector to save
   */
  inline void setVelocity0(SiconosVector *v)
  {
    if (velocity0Node == NULL)
    {
      velocity0Node = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LNLDS_VELOCITY0, *v);
    }
    else
    {
      //SiconosDOMTreeTools::setSiconosVectorValue(velocity0Node, *v);
      SiconosDOMTreeTools::setSiconosVectorNodeValue(velocity0Node, *v);
    }
  }

  // === VelocityMemory ===
  /** \fn SiconosMemoryXML* getVelocityMemoryXML()
   *   \brief Returns the velocityMemoryXML* of the DSXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXML* getVelocityMemoryXML()
  {
    return velocityMemoryXML;
  }

  /** \fn void setVelocityMemory(SiconosMemory* smem)
   *   \brief allows to save the velocityMemory of the LagrangianDSXML
   *   \param SiconosMemory* : SiconosMemory to save
   */
  void setVelocityMemory(SiconosMemory* smem);

  // === NNL ===
  /** \fn inline string getNNLPlugin()
   *   \brief Return the NNL Plugin name of the LagrangianDSXML
   *   \return The NNL Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline std::string getNNLPlugin()
  {
    if (!isNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getNNLPlugin : NNL is not calculated from a plugin ; NNL vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(NNLNode, LNLDS_VECTORPLUGIN);
  }

  /** \fn SimpleVector getNNLVector()
   *   \brief Return the NNL vector of the LagrangianDSXML
   *   \return SimpleVector : NNL vector of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SimpleVector getNNLVector()
  {
    if (isNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getNNLVector : NNL vector is not given ; NNL is calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosVectorValue(NNLNode);
  }

  /** \fn void setNNLPlugin(string plugin)
   *   \brief allows to save the NNL plugin of the LagrangianDSXML
   *   \param string : ths string which contains the name and the location of the plugin
   */
  inline void setNNLPlugin(std::string plugin)
  {
    if (NNLNode == NULL)
    {
      NNLNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_QNLINERTIA);
      xmlNewProp(NNLNode, (xmlChar*)(LNLDS_VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(NNLNode, LNLDS_VECTORPLUGIN, plugin);
  }

  /** \fn void setNNLVector(SiconosVector *v)
   *   \brief allows to save the NNL vector of the LagrangianDSXML
   *   \return The NNL SiconosVector to save
   */
  inline void setNNLVector(SiconosVector *v)
  {
    if (NNLNode == NULL)
    {
      NNLNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LNLDS_QNLINERTIA, *v);
    }
    else
      SiconosDOMTreeTools::setSiconosVectorNodeValue(NNLNode, *v);
  }

  // === FInt ===
  /** \fn inline string getFintPlugin()
   *   \brief Return the Fint Plugin name of the LagrangianDSXML
   *   \return The Fint Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline std::string getFintPlugin()
  {
    if (!isFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFintPlugin : Fint is not calculated from a plugin ; Fint vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(FintNode, LNLDS_VECTORPLUGIN);
  }

  /** \fn SimpleVector getFintVector()
   *   \brief Return the internal forces vector of the LagrangianDSXML
   *   \return SimpleVector : Fint SiconosVector of the LagrangianDSXML
   *   \exception XMLException
   */
  inline SimpleVector getFintVector()
  {
    if (isFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFintVector : Fint vector is not given ; Fint is calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosVectorValue(FintNode);
  }

  /** \fn void setFintVector(SiconosVector *v)
   *   \brief allows to save the Fint vector of the LagrangianDSXML
   *   \param The Fint SiconosVector to save
   */
  inline void setFintVector(SiconosVector *v)
  {
    if (hasFint())
    {
      FintNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, LNLDS_FINT, *v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(FintNode, *v);
  }

  /** \fn void setFextPlugin(string plugin)
   *   \brief allows to save the Fext plugin of the LagrangianDSXML
   *   \param string : ths string which contains the name and the location of the plugin
   */
  inline void setFintPlugin(std::string plugin)
  {
    if (FintNode == NULL)
    {
      FintNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_FINT);
      xmlNewProp(FintNode, (xmlChar*)(LNLDS_VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(FintNode, LNLDS_VECTORPLUGIN, plugin);
  }

  // === Fext ===
  /** \fn inline string getFextPlugin()
   *   \brief Return the Fext Plugin name of the LagrangianDSXML
   *   \return The Fext Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline std::string getFextPlugin()
  {
    if (!isFextPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFextPlugin : Fext is not calculated from a plugin ; Fext vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(FextNode, LNLDS_VECTORPLUGIN);
  }

  /** \fn SimpleVector getFextVector()
   *   \brief Return the external forces vector of the LagrangianDSXML
   *   \return SimpleVector : Fext vector of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SimpleVector getFextVector()
  {
    if (isFextPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFextVector : Fext matrix is not given ; Fext is calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosVectorValue(FextNode);
  }

  /** \fn void setFextVector(SiconosVector *v)
   *   \brief allows to save the Fint vector of the LagrangianDSXML
   *   \param The Fint SiconosVector to save
   */
  inline void setFextVector(SiconosVector *v)
  {
    //SiconosDOMTreeTools::setSiconosVectorValue(FextNode, *v);
    SiconosDOMTreeTools::setSiconosVectorNodeValue(FextNode, *v);
  }

  /** \fn void setFextPlugin(string plugin)
   *   \brief allows to save the Fext plugin of the LagrangianDSXML
   *   \param string : ths string which contains the name and the location of the plugin
   */
  inline void setFextPlugin(std::string plugin)
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

  /** \fn inline string getJacobianQFintPlugin()
   *   \brief Return the JacobianQFint Plugin name of the LagrangianDSXML
   *   \return The JacobianQFint Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline std::string getJacobianQFintPlugin()
  {
    if (!isJacobianQFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianQFintPlugin : JacobianQFint is not calculated from a plugin ; JacobianQFint matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianQFintNode, LNLDS_MATRIXPLUGIN);
  }

  /** \fn SimpleMatrix getJacobianQFintMatrix()
   *   \brief Return the JacobianQFint matrix of the LagrangianDSXML
   *   \return The JacobianQFint SimpleMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SimpleMatrix getJacobianQFintMatrix()
  {
    if (isJacobianQFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianQFintMatrix : JacobianQFint matrix is not given ; JacobianQFint is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianQFintNode);
  }

  /** \fn void setJacobianQFintPlugin(string plugin)
   *   \brief allows to save the jacobianQFint plugin of the LagrangianDSXML
   *   \param string : ths string which contains the name and the location of the plugin
   */
  inline void setJacobianQFintPlugin(std::string plugin)
  {
    if (jacobianQFintNode == NULL)
    {
      jacobianQFintNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_JACOBIANQFINT);
      xmlNewProp(jacobianQFintNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(jacobianQFintNode, LNLDS_JACOBIANQFINT, plugin);
  }

  /** \fn void setJacobianQFintMatrix(SiconosMatrix *m)
   *   \brief allows to save the JacobianQFint matrix of the LagrangianDSXML
   *   \return The JacobianQFint SiconosMatrix to save
   */
  inline void setJacobianQFintMatrix(SiconosMatrix *m)
  {
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianQFintNode, *m);
  }

  /** \fn inline string getJacobianVelocityFintPlugin()
   *   \brief Return the JacobianVelocityFint Plugin name of the LagrangianDSXML
   *   \return The JacobianVelocityFint Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline std::string getJacobianVelocityFintPlugin()
  {
    if (!isJacobianVelocityFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityFintPlugin : JacobianVelocityFint is not calculated from a plugin ; JacobianVelocityFint matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianVelocityFintNode, LNLDS_MATRIXPLUGIN);
  }

  /** \fn SimpleMatrix getJacobianVelocityFintMatrix()
   *   \brief Return the JacobianVelocityFint matrix of the LagrangianDSXML
   *   \return The JacobianVelocityFint SimpleMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SimpleMatrix getJacobianVelocityFintMatrix()
  {
    if (isJacobianVelocityFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityFintMatrix : JacobianVelocityFint matrix is not given ; JacobianVelocityFint is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianVelocityFintNode);
  }

  /** \fn void setJacobianVelocityFintPlugin(string plugin)
   *   \brief allows to save the jacobianVelocityFint plugin of the LagrangianDSXML
   *   \param string : ths string which contains the name and the location of the plugin
   */
  inline void setJacobianVelocityFintPlugin(std::string plugin)
  {
    if (jacobianVelocityFintNode == NULL)
    {
      jacobianVelocityFintNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_JACOBIANVELOCITYFINT);
      xmlNewProp(jacobianVelocityFintNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(jacobianVelocityFintNode, LNLDS_JACOBIANVELOCITYFINT, plugin);
  }

  /** \fn void setJacobianVelocityFintMatrix(SiconosMatrix *m)
   *   \brief allows to save the JacobianVelocityFint matrix of the LagrangianDSXML
   *   \return The JacobianVelocityFint SiconosMatrix to save
   */
  inline void setJacobianVelocityFintMatrix(SiconosMatrix *m)
  {
    //SiconosDOMTreeTools::setSiconosMatrixValue(jacobianVelocityFintNode, *m);
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianVelocityFintNode, *m);
  }

  /** \fn inline string getJacobianQQPlugin()
   *   \brief Return the JacobianQQ Plugin name of the LagrangianDSXML
   *   \return The JacobianQQ Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline std::string getJacobianQNNLPlugin()
  {
    if (!isJacobianQNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianQNNLPlugin : JacobianQNNL is not calculated from a plugin ; JacobianQNNL matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianQNNLNode, LNLDS_MATRIXPLUGIN);
  }

  /** \fn SimpleMatrix getJacobianQQMatrix()
   *   \brief Return the JacobianQQ matrix of the LagrangianDSXML
   *   \return The JacobianQQ SimpleMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SimpleMatrix getJacobianQNNLMatrix()
  {
    if (isJacobianQNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianQNNLMatrix : JacobianQNNL matrix is not given ; JacobianQNNL is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianQNNLNode);
  }

  /** \fn void setJacobianQNNLPlugin(string plugin)
   *   \brief allows to save the jacobianQNNL plugin of the LagrangianDSXML
   *   \param string : ths string which contains the name and the location of the plugin
   */
  inline void setJacobianQNNLPlugin(std::string plugin)
  {
    if (jacobianQNNLNode == NULL)
    {
      jacobianQNNLNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_JACOBIANQQNLINERTIA);
      xmlNewProp(jacobianQNNLNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(jacobianQNNLNode, LNLDS_JACOBIANQQNLINERTIA, plugin);
  }

  /** \fn void setJacobianQQMatrix(SiconosMatrix *m)
   *   \brief allows to save the JacobianQQ matrix of the LagrangianDSXML
   *   \return The JacobianQQ SiconosMatrix to save
   */
  inline void setJacobianQNNLMatrix(SiconosMatrix *m)
  {
    //SiconosDOMTreeTools::setSiconosMatrixValue(jacobianQNNLNode, *m);
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianQNNLNode, *m);
  }

  /** \fn inline string getJacobianVelocityNNLPlugin()
   *   \brief Return the JacobianVelocityNNL Plugin name of the LagrangianDSXML
   *   \return The JacobianVelocityNNL Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline std::string getJacobianVelocityNNLPlugin()
  {
    if (!isJacobianVelocityNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityNNLPlugin : JacobianVelocityNNL is not calculated from a plugin ; JacobianVelocityNNL matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianVelocityNNLNode, LNLDS_MATRIXPLUGIN);
  }

  /** \fn SimpleMatrix getJacobianVelocityNNLMatrix()
   *   \brief Return the JacobianVelocityNNL matrix of the LagrangianDSXML
   *   \return The JacobianVelocityNNL SimpleMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SimpleMatrix getJacobianVelocityNNLMatrix()
  {
    if (isJacobianVelocityNNLPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityNNLMatrix : JacobianVelocityNNL matrix is not given ; JacobianVelocityNNL is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianVelocityNNLNode);
  }

  /** \fn void setJacobianVelocityNNLPlugin(string plugin)
   *   \brief allows to save the jacobianVelocityNNLPlugin plugin of the LagrangianDSXML
   *   \param string : ths string which contains the name and the location of the plugin
   */
  inline void setJacobianVelocityNNLPlugin(std::string plugin)
  {
    if (jacobianVelocityNNLNode == NULL)
    {
      jacobianVelocityNNLNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_JACOBIANVELOCITYQNLINERTIA);
      xmlNewProp(jacobianVelocityNNLNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(jacobianVelocityNNLNode, LNLDS_JACOBIANVELOCITYQNLINERTIA, plugin);
  }

  /** \fn void setJacobianVelocityNNLMatrix(SiconosMatrix *m)
   *   \brief allows to save the JacobianVelocityNNL matrix of the LagrangianDSXML
   *   \return The JacobianVelocityNNL SiconosMatrix to save
   */
  inline void setJacobianVelocityNNLMatrix(SiconosMatrix *m)
  {
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianVelocityNNLNode, *m);
  }

  /** \fn inline string getMPlugin()
   *   \brief Return the M Plugin name of the LagrangianDSXML
   *   \return The M Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline std::string getMPlugin()
  {
    if (!isMPlugin())
      XMLException::selfThrow("LagrangianDSXML - getMPlugin : M is not calculated from a plugin ; M matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(MNode, LNLDS_MATRIXPLUGIN);
  }

  /** \fn SimpleMatrix getMMatrix()
   *   \brief Return the M matrix of the LagrangianDSXML
   *   \return The M SimpleMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SimpleMatrix getMMatrix()
  {
    if (isMPlugin())
      XMLException::selfThrow("LagrangianDSXML - getMMatrix : M matrix is not given ; M is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(MNode);
  }

  /** \fn void setMPlugin(string plugin)
   *   \brief allows to save the jacobianVelocityNNLPlugin plugin of the LagrangianDSXML
   *   \param string : ths string which contains the name and the location of the plugin
   */
  inline void setMPlugin(std::string plugin)
  {
    if (MNode == NULL)
    {
      MNode = SiconosDOMTreeTools::createSingleNode(rootDynamicalSystemXMLNode, LNLDS_M);
      xmlNewProp(MNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else
      SiconosDOMTreeTools::setStringAttributeValue(MNode, LNLDS_MATRIXPLUGIN, plugin);
  }

  /** \fn void setMMatrix(SiconosMatrix *m)
   *   \brief allows to save the M matrix of the LagrangianDSXML
   *   \return The M SiconosMatrix to save
   */
  inline void setMMatrix(SiconosMatrix *m)
  {
    if (MNode == NULL)
    {
      MNode = SiconosDOMTreeTools::createMatrixNode(rootDynamicalSystemXMLNode, LNLDS_M, *m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(MNode, *m);
  }


  /** \fn int getNdof()
   *   \brief Return the ndof for the LagrangianDSXML
   *   \return The ndof integer for the LagrangianDSXML
   */
  inline int getNdof()
  {
    return  SiconosDOMTreeTools::getIntegerContentValue(ndofNode);
  }

  /** \fn void setNdof(int i)
   *   \brief allows to save the ndof for the LagrangianDSXML
   *   \return The ndof integer to save
   */
  inline void setNdof(int i)
  {
    if (ndofNode == NULL)
    {
      ndofNode = SiconosDOMTreeTools::createIntegerNode(rootDynamicalSystemXMLNode, LNLDS_NDOF, i);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(ndofNode, i);
  }


  /** \fn bool isMPlugin()
   *   \brief Return true if M is calculated from a plugin
   *   \return True if M is calculated from plugin
   */
  inline bool isMPlugin()
  {
    return xmlHasProp((xmlNodePtr)MNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str());
  }

  /** \fn bool isMMatrix()
   *   \brief Return true if M only given by a Matrix
   *   \return True if M is given by a Matrix
   */
  inline bool isMMatrix()
  {
    return !(xmlHasProp((xmlNodePtr)MNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str()));
  }

  /** \fn bool isQLNInertiaPlugin()
   *   \brief Return true if QLNInertia is calculated from a plugin
   *   \return True if QLNInertia is calculated from plugin
   */
  inline bool isNNLPlugin()
  {
    return xmlHasProp((xmlNodePtr)NNLNode, (xmlChar *) LNLDS_VECTORPLUGIN.c_str());
  }

  /** \fn bool isFintPlugin()
   *   \brief Return true if Fint is calculated from a plugin
   *   \return True if Fint is calculated from plugin
   */
  inline bool isFintPlugin()
  {
    return xmlHasProp((xmlNodePtr)FintNode, (xmlChar *) LNLDS_VECTORPLUGIN.c_str());
  }

  /** \fn bool isFextPlugin()
   *   \brief Return true if Fext is calculated from a plugin
   *   \return True if Fext is calculated from plugin
   */
  inline bool isFextPlugin()
  {
    return xmlHasProp((xmlNodePtr)FextNode, (xmlChar *) LNLDS_VECTORPLUGIN.c_str());
  }

  /** \fn bool isJacobianQFintPlugin()
   *   \brief Return true if JacobianQFint is calculated from a plugin
   *   \return True if JacobianQFint is calculated from plugin
   */
  inline bool isJacobianQFintPlugin()
  {
    return xmlHasProp((xmlNodePtr)jacobianQFintNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str());
  }

  /** \fn bool isJacobianVelocityFintPlugin()
   *   \brief Return true if JacobianVelocityFint is calculated from a plugin
   *   \return True if JacobianVelocityFint is calculated from plugin
   */
  inline bool isJacobianVelocityFintPlugin()
  {
    return xmlHasProp((xmlNodePtr)jacobianVelocityFintNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str());
  }

  /** \fn bool isJacobianQQPlugin()
   *   \brief Return true if JacobianQQ is calculated from a plugin
   *   \return True if JacobianQQ is calculated from plugin
   */
  inline bool isJacobianQNNLPlugin()
  {
    return xmlHasProp((xmlNodePtr)jacobianQNNLNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str());
  }

  /** \fn bool isJacobianVelocityNNLPlugin()
   *   \brief Return true if JacobianVelocityNNL is calculated from a plugin
   *   \return True if JacobianVelocityNNL is calculated from plugin
   */
  inline bool isJacobianVelocityNNLPlugin()
  {
    return xmlHasProp((xmlNodePtr)jacobianVelocityNNLNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str());
  }


  /** \fn bool hasMass()
   *  \brief determines if Mass is defined in the DOM tree
   *  \return bool : true if Mass is defined, false otherwise
   */
  inline bool hasMass()
  {
    return (MNode != NULL);
  }

  /** \fn bool hasFint()
   *  \brief determines if Fint is defined in the DOM tree
   *  \return bool : true if Fint is defined, false otherwise
   */
  inline bool hasFint()
  {
    return (FintNode != NULL);
  }

  /** \fn bool hasFext()
   *  \brief determines if Fext is defined in the DOM tree
   *  \return bool : true if Fext is defined, false otherwise
   */
  inline bool hasFext()
  {
    return (FextNode != NULL);
  }

  /** \fn bool hasJacobianQFint()
   *  \brief determines if jacobianQFint is defined in the DOM tree
   *  \return bool : true if jacobianQFint is defined, false otherwise
   */
  inline bool hasJacobianQFint()
  {
    return (jacobianQFintNode != NULL);
  }

  /** \fn bool hasJacobianVelocityFint()
   *  \brief determines if jacobianVelocityFint is defined in the DOM tree
   *  \return bool : true if jacobianVelocityFint is defined, false otherwise
   */
  inline bool hasJacobianVelocityFint()
  {
    return (jacobianVelocityFintNode != NULL);
  }

  /** \fn bool hasJacobianQNNL()
   *  \brief determines if jacobianQNNL is defined in the DOM tree
   *  \return bool : true if jacobianQNNL is defined, false otherwise
   */
  inline bool hasJacobianQNNL()
  {
    return (jacobianQNNLNode != NULL);
  }

  /** \fn bool hasJacobianVelocityNNL()
   *  \brief determines if jacobianVelocityNNL is defined in the DOM tree
   *  \return bool : true if jacobianVelocityNNL is defined, false otherwise
   */
  inline bool hasJacobianVelocityNNL()
  {
    return (jacobianVelocityNNLNode != NULL);
  }

  /** \fn bool hasNNL()
   *  \brief determines if NNL is defined in the DOM tree
   *  \return bool : true if NNL is defined, false otherwise
   */
  inline bool hasNNL()
  {
    return (NNLNode != NULL);
  }

  /** \fn bool hasQMemory()
   *  \brief returns true if qMemoryNode is defined
   *  \return true if qMemoryNode is defined
   */
  inline bool hasQMemory()
  {
    return (qMemoryNode != NULL);
  }

  /** \fn bool hasVelocityMemory()
   *  \brief returns true if velocityMemoryNode is defined
   *  \return true if velocityMemoryNode is defined
   */
  inline bool hasVelocityMemory()
  {
    return (velocityMemoryNode != NULL);
  }

  /** \fn bool hasQ()
   *  \brief returns true if qNode is defined
   *  \return true if qNode is defined
   */
  inline bool hasQ()
  {
    return (qNode != NULL);
  }

  /** \fn bool hasVelocity()
   *  \brief returns true if velocityNode is defined
   *  \return true if velocityNode is defined
   */
  inline bool hasVelocity()
  {
    return (velocityNode != NULL);
  }


  /** \fn void updateDynamicalSystemXML(xmlNode*, DynamicalSystem*, BoundaryCondition*)
   *   \brief makes the operations to add a DynamicalSystem to the NonSmoothDynamicalSystemXML
   *   \param xmlNode* : the root node of this DynamicalSystem
   *   \param DynamicalSystem* : the DynamicalSystem of this DynamicalSystemXML
   *   \param BoundaryCondition* : the BoundaryCondition of the DS if the NonSmoothDynamicalSystem is BVP (optional)
   */
  void updateDynamicalSystemXML(xmlNode*, DynamicalSystem*, BoundaryCondition* bc = NULL);


protected:

  SiconosMemoryXML * qMemoryXML;
  SiconosMemoryXML * velocityMemoryXML;

  //Nodes

  xmlNode * qNode;
  xmlNode * q0Node;
  xmlNode * qMemoryNode;

  xmlNode * velocityNode;
  xmlNode * velocity0Node;
  xmlNode * velocityMemoryNode;

  xmlNode * NNLNode;
  xmlNode * FintNode;
  xmlNode * FextNode;

  xmlNode * jacobianQFintNode;
  xmlNode * jacobianVelocityFintNode;
  xmlNode * jacobianQNNLNode;
  xmlNode * jacobianVelocityNNLNode;


  xmlNode * MNode;
  xmlNode * ndofNode;

};

#endif
