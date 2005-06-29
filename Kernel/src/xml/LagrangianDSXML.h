
/** \class LagrangianDSXML
 *   \brief This class manages Lagrangian NLDS data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 05/11/2004
 *
 *
 *
 * LagrangianDSXML allows to manage data of a LagrangianDS DOM tree.
 */

#ifndef __LAGRANGIANNLDSXML__
#define __LAGRANGIANNLDSXML__

#include "DSXML.h"

class DSXML;

const std::string LNLDS_Q = "q";
const std::string LNLDS_Q0 = "q0";
const std::string LNLDS_QMEMORY = "qMemory";

const std::string LNLDS_VELOCITY = "Velocity";
const std::string LNLDS_VELOCITY0 = "Velocity0";
const std::string LNLDS_VELOCITYMEMORY = "VelocityMemory";

const std::string LNLDS_QNLINERTIA = "QNLInertia";
const std::string LNLDS_FINT = "Fint";
const std::string LNLDS_FEXT = "Fext";

const std::string LNLDS_JACOBIANQFINT = "JacobianQFint";
const std::string LNLDS_JACOBIANVELOCITYFINT = "JacobianVelocityFint";
const std::string LNLDS_JACOBIANQQNLINERTIA = "JacobianQQNLInertia";
const std::string LNLDS_JACOBIANVELOCITYQNLINERTIA = "JacobianVelocityQNLInertia";

const std::string LNLDS_M = "M";
const std::string LNLDS_NDOF = "ndof";
const std::string LNLDS_MATRIXPLUGIN = "matrixPlugin";
const std::string LNLDS_VECTORPLUGIN = "vectorPlugin";
//#include "XMLTagsName.h"

#include "check.h"

class LagrangianDSXML : public DSXML
{
public:
  LagrangianDSXML();

  virtual ~LagrangianDSXML();

  /** \fn LagrangianDSXML(xmlNode * LagrangianDSNode, bool isBVP)
   *   \brief Build a LagrangianDSXML object from a DOM tree describing a DS
   *   \param xmlNode * LagrangianDSNode : the LagrangianDS DOM tree
   *   \param bool isBVP : if NSDS is BVP LagrangianDS have boundary condition
   */
  LagrangianDSXML(xmlNode * LagrangianDSNode, bool isBVP);

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
    if (hasQ() == false)
    {
      qNode = SiconosDOMTreeTools::createVectorNode(rootDSXMLNode, LNLDS_Q, *v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(qNode, *v);
  }

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
    {
      q0Node = SiconosDOMTreeTools::createVectorNode(rootDSXMLNode, LNLDS_Q0, *v);
    }
    else
    {
      //SiconosDOMTreeTools::setSiconosVectorValue(q0Node, *v);
      SiconosDOMTreeTools::setSiconosVectorNodeValue(q0Node, *v);
    }
  }


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
  inline void setQMemory(SiconosMemory* smem)
  {
    if (hasQMemory() == false)
    {
      qMemoryXML = new SiconosMemoryXML(NULL, rootDSXMLNode, LNLDS_QMEMORY);
      qMemoryNode = qMemoryXML->getSiconosMemoryXMLNode();

      qMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      qMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
    else
    {
      qMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      qMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
  }

  /** \fn SimpleVector getVelocity()
   *   \brief Return the velocity of the LagrangianDSXML
   *   \return SimpleVector :  velocity vector of the LagrangianDSXML
   */
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
    if (hasVelocity() == false)
    {
      velocityNode = SiconosDOMTreeTools::createVectorNode(rootDSXMLNode, LNLDS_VELOCITY, *v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(velocityNode, *v);
    //SiconosDOMTreeTools::setSiconosVectorValue(velocityNode, v);
  }

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
      velocity0Node = SiconosDOMTreeTools::createVectorNode(rootDSXMLNode, LNLDS_VELOCITY0, *v);
    }
    else
    {
      //SiconosDOMTreeTools::setSiconosVectorValue(velocity0Node, *v);
      SiconosDOMTreeTools::setSiconosVectorNodeValue(velocity0Node, *v);
    }
  }

  //    /** \fn SiconosMemory getVelocityMemory()
  //    *   \brief Return the velocityMemory of the LagrangianDSXML
  //    *   \return SiconosMemory velocityMemory of the LagrangianDSXML
  //    */
  //    SiconosMemory getVelocityMemory();

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
  inline void setVelocityMemory(SiconosMemory* smem)
  {
    if (hasVelocityMemory() == false)
    {
      velocityMemoryXML = new SiconosMemoryXML(NULL, rootDSXMLNode, LNLDS_VELOCITYMEMORY);
      velocityMemoryNode = velocityMemoryXML->getSiconosMemoryXMLNode();

      velocityMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      velocityMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
    else
    {
      velocityMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      velocityMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
  }

  /** \fn inline string getQNLInertiaPlugin()
   *   \brief Return the QNLInertia Plugin name of the LagrangianDSXML
   *   \return The QNLInertia Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline std::string getQNLInertiaPlugin()
  {
    if (!isQNLInertiaPlugin())
      XMLException::selfThrow("LagrangianDSXML - getQNLInertiaPlugin : QNLInertia is not calculated from a plugin ; QNLInertia vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(QNLInertiaNode, LNLDS_VECTORPLUGIN);
  }

  /** \fn SimpleVector getQNLInertiaVector()
   *   \brief Return the QNLInertia vector of the LagrangianDSXML
   *   \return SimpleVector : QNLInertia vector of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SimpleVector getQNLInertiaVector()
  {
    if (isQNLInertiaPlugin())
      XMLException::selfThrow("LagrangianDSXML - getQNLInertiaVector : QNLInertia vector is not given ; QNLInertia is calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosVectorValue(QNLInertiaNode);
  }

  /** \fn void setQNLInertiaPlugin(string plugin)
   *   \brief allows to save the QNLInertia plugin of the LagrangianDSXML
   *   \param string : ths string which contains the name and the location of the plugin
   */
  inline void setQNLInertiaPlugin(std::string plugin)
  {
    if (QNLInertiaNode == NULL)
    {
      QNLInertiaNode = SiconosDOMTreeTools::createSingleNode(rootDSXMLNode, LNLDS_QNLINERTIA);
      xmlNewProp(QNLInertiaNode, (xmlChar*)(LNLDS_VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(QNLInertiaNode, LNLDS_VECTORPLUGIN, plugin);
  }

  /** \fn void setQNLInertiaVector(SiconosVector *v)
   *   \brief allows to save the QNLInertia vector of the LagrangianDSXML
   *   \return The QNLInertia SiconosVector to save
   */
  inline void setQNLInertiaVector(SiconosVector *v)
  {
    if (QNLInertiaNode == NULL)
    {
      QNLInertiaNode = SiconosDOMTreeTools::createVectorNode(rootDSXMLNode, LNLDS_QNLINERTIA, *v);
    }
    else
    {
      //SiconosDOMTreeTools::setSiconosVectorValue(QNLInertiaNode, *v);
      SiconosDOMTreeTools::setSiconosVectorNodeValue(QNLInertiaNode, *v);
    }
  }


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
      FintNode = SiconosDOMTreeTools::createVectorNode(rootDSXMLNode, LNLDS_FINT, *v);
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
      FintNode = SiconosDOMTreeTools::createSingleNode(rootDSXMLNode, LNLDS_FINT);
      xmlNewProp(FintNode, (xmlChar*)(LNLDS_VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(FintNode, LNLDS_VECTORPLUGIN, plugin);
  }

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
      FextNode = SiconosDOMTreeTools::createSingleNode(rootDSXMLNode, LNLDS_FEXT);
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

  /** \fn SiconosMatrix getJacobianQFintMatrix()
   *   \brief Return the JacobianQFint matrix of the LagrangianDSXML
   *   \return The JacobianQFint SiconosMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SiconosMatrix getJacobianQFintMatrix()
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
      jacobianQFintNode = SiconosDOMTreeTools::createSingleNode(rootDSXMLNode, LNLDS_JACOBIANQFINT);
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

  /** \fn SiconosMatrix getJacobianVelocityFintMatrix()
   *   \brief Return the JacobianVelocityFint matrix of the LagrangianDSXML
   *   \return The JacobianVelocityFint SiconosMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SiconosMatrix getJacobianVelocityFintMatrix()
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
      jacobianVelocityFintNode = SiconosDOMTreeTools::createSingleNode(rootDSXMLNode, LNLDS_JACOBIANVELOCITYFINT);
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
  inline std::string getJacobianQQNLInertiaPlugin()
  {
    if (!isJacobianQQNLInertiaPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianQQNLInertiaPlugin : JacobianQQNLInertia is not calculated from a plugin ; JacobianQQNLInertia matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianQQNLInertiaNode, LNLDS_MATRIXPLUGIN);
  }

  /** \fn SiconosMatrix getJacobianQQMatrix()
   *   \brief Return the JacobianQQ matrix of the LagrangianDSXML
   *   \return The JacobianQQ SiconosMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SiconosMatrix getJacobianQQNLInertiaMatrix()
  {
    if (isJacobianQQNLInertiaPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianQQNLInertiaMatrix : JacobianQQNLInertia matrix is not given ; JacobianQQNLInertia is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianQQNLInertiaNode);
  }

  /** \fn void setJacobianQQNLInertiaPlugin(string plugin)
   *   \brief allows to save the jacobianQQNLInertia plugin of the LagrangianDSXML
   *   \param string : ths string which contains the name and the location of the plugin
   */
  inline void setJacobianQQNLInertiaPlugin(std::string plugin)
  {
    if (jacobianQQNLInertiaNode == NULL)
    {
      jacobianQQNLInertiaNode = SiconosDOMTreeTools::createSingleNode(rootDSXMLNode, LNLDS_JACOBIANQQNLINERTIA);
      xmlNewProp(jacobianQQNLInertiaNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(jacobianQQNLInertiaNode, LNLDS_JACOBIANQQNLINERTIA, plugin);
  }

  /** \fn void setJacobianQQMatrix(SiconosMatrix *m)
   *   \brief allows to save the JacobianQQ matrix of the LagrangianDSXML
   *   \return The JacobianQQ SiconosMatrix to save
   */
  inline void setJacobianQQNLInertiaMatrix(SiconosMatrix *m)
  {
    //SiconosDOMTreeTools::setSiconosMatrixValue(jacobianQQNLInertiaNode, *m);
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianQQNLInertiaNode, *m);
  }

  /** \fn inline string getJacobianVelocityQNLInertiaPlugin()
   *   \brief Return the JacobianVelocityQNLInertia Plugin name of the LagrangianDSXML
   *   \return The JacobianVelocityQNLInertia Plugin name of the LagrangianDSXML
   *  \exception XMLException
   */
  inline std::string getJacobianVelocityQNLInertiaPlugin()
  {
    if (!isJacobianVelocityQNLInertiaPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityQNLInertiaPlugin : JacobianVelocityQNLInertia is not calculated from a plugin ; JacobianVelocityQNLInertia matrix is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianVelocityQNLInertiaNode, LNLDS_MATRIXPLUGIN);
  }

  /** \fn SiconosMatrix getJacobianVelocityQNLInertiaMatrix()
   *   \brief Return the JacobianVelocityQNLInertia matrix of the LagrangianDSXML
   *   \return The JacobianVelocityQNLInertia SiconosMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SiconosMatrix getJacobianVelocityQNLInertiaMatrix()
  {
    if (isJacobianVelocityQNLInertiaPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityQNLInertiaMatrix : JacobianVelocityQNLInertia matrix is not given ; JacobianVelocityQNLInertia is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianVelocityQNLInertiaNode);
  }

  /** \fn void setJacobianVelocityQNLInertiaPlugin(string plugin)
   *   \brief allows to save the jacobianVelocityQNLInertiaPlugin plugin of the LagrangianDSXML
   *   \param string : ths string which contains the name and the location of the plugin
   */
  inline void setJacobianVelocityQNLInertiaPlugin(std::string plugin)
  {
    if (jacobianVelocityQNLInertiaNode == NULL)
    {
      jacobianVelocityQNLInertiaNode = SiconosDOMTreeTools::createSingleNode(rootDSXMLNode, LNLDS_JACOBIANVELOCITYQNLINERTIA);
      xmlNewProp(jacobianVelocityQNLInertiaNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(jacobianVelocityQNLInertiaNode, LNLDS_JACOBIANVELOCITYQNLINERTIA, plugin);
  }

  /** \fn void setJacobianVelocityQNLInertiaMatrix(SiconosMatrix *m)
   *   \brief allows to save the JacobianVelocityQNLInertia matrix of the LagrangianDSXML
   *   \return The JacobianVelocityQNLInertia SiconosMatrix to save
   */
  inline void setJacobianVelocityQNLInertiaMatrix(SiconosMatrix *m)
  {
    SiconosDOMTreeTools::setSiconosMatrixNodeValue(jacobianVelocityQNLInertiaNode, *m);
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

  /** \fn SiconosMatrix getMMatrix()
   *   \brief Return the M matrix of the LagrangianDSXML
   *   \return The M SiconosMatrix of the LagrangianDSXML
   *  \exception XMLException
   */
  inline SiconosMatrix getMMatrix()
  {
    if (isMPlugin())
      XMLException::selfThrow("LagrangianDSXML - getMMatrix : M matrix is not given ; M is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(MNode);
  }

  /** \fn void setMPlugin(string plugin)
   *   \brief allows to save the jacobianVelocityQNLInertiaPlugin plugin of the LagrangianDSXML
   *   \param string : ths string which contains the name and the location of the plugin
   */
  inline void setMPlugin(std::string plugin)
  {
    if (MNode == NULL)
    {
      MNode = SiconosDOMTreeTools::createSingleNode(rootDSXMLNode, LNLDS_M);
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
      MNode = SiconosDOMTreeTools::createMatrixNode(rootDSXMLNode, LNLDS_M, *m);
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
      ndofNode = SiconosDOMTreeTools::createIntegerNode(rootDSXMLNode, LNLDS_NDOF, i);
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
  inline bool isQNLInertiaPlugin()
  {
    return xmlHasProp((xmlNodePtr)QNLInertiaNode, (xmlChar *) LNLDS_VECTORPLUGIN.c_str());
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
  inline bool isJacobianQQNLInertiaPlugin()
  {
    return xmlHasProp((xmlNodePtr)jacobianQQNLInertiaNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str());
  }

  /** \fn bool isJacobianVelocityQNLInertiaPlugin()
   *   \brief Return true if JacobianVelocityQNLInertia is calculated from a plugin
   *   \return True if JacobianVelocityQNLInertia is calculated from plugin
   */
  inline bool isJacobianVelocityQNLInertiaPlugin()
  {
    return xmlHasProp((xmlNodePtr)jacobianVelocityQNLInertiaNode, (xmlChar *) LNLDS_MATRIXPLUGIN.c_str());
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

  /** \fn bool hasJacobianQQNLInertia()
   *  \brief determines if jacobianQQNLInertia is defined in the DOM tree
   *  \return bool : true if jacobianQQNLInertia is defined, false otherwise
   */
  inline bool hasJacobianQQNLInertia()
  {
    return (jacobianQQNLInertiaNode != NULL);
  }

  /** \fn bool hasJacobianVelocityQNLInertia()
   *  \brief determines if jacobianVelocityQNLInertia is defined in the DOM tree
   *  \return bool : true if jacobianVelocityQNLInertia is defined, false otherwise
   */
  inline bool hasJacobianVelocityQNLInertia()
  {
    return (jacobianVelocityQNLInertiaNode != NULL);
  }

  /** \fn bool hasQNLInertia()
   *  \brief determines if QNLInertia is defined in the DOM tree
   *  \return bool : true if QNLInertia is defined, false otherwise
   */
  inline bool hasQNLInertia()
  {
    return (QNLInertiaNode != NULL);
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
   *   \brief makes the operations to add a DynamicalSystem to the NSDSXML
   *   \param xmlNode* : the root node of this DynamicalSystem
   *   \param DynamicalSystem* : the DynamicalSystem of this DSXML
   *   \param BoundaryCondition* : the BoundaryCondition of the DS if the NSDS is BVP (optional)
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

  xmlNode * QNLInertiaNode;
  xmlNode * FintNode;
  xmlNode * FextNode;

  xmlNode * jacobianQFintNode;
  xmlNode * jacobianVelocityFintNode;
  xmlNode * jacobianQQNLInertiaNode;
  xmlNode * jacobianVelocityQNLInertiaNode;


  xmlNode * MNode;
  xmlNode * ndofNode;


  //Methods

  /** \fn loadLagrangianDSProperties()
   *   \brief load the different properties of a LagrangianDS
   *   \exception XMLException : if a property of the LagrangianDS lacks in the DOM tree
   */
  void loadLagrangianDSProperties();

};

#endif
