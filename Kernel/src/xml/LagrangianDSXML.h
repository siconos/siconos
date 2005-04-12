
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


#include <vector>
#include <string>
#include <libxml/tree.h>

#include "DSXML.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include "SiconosMemory.h"
#include "XMLException.h"
#include "SiconosDOMTreeTools.h"


//using namespace std;


class DSXML;

const string LNLDS_Q = "q";
const string LNLDS_Q0 = "q0";
const string LNLDS_QMEMORY = "qMemory";

const string LNLDS_VELOCITY = "Velocity";
const string LNLDS_VELOCITY0 = "Velocity0";
const string LNLDS_VELOCITYMEMORY = "VelocityMemory";

const string LNLDS_QNLINERTIA = "QNLInertia";
const string LNLDS_FINT = "Fint";
const string LNLDS_FEXT = "Fext";

const string LNLDS_JACOBIANQFINT = "JacobianQFint";
const string LNLDS_JACOBIANVELOCITYFINT = "JacobianVelocityFint";
const string LNLDS_JACOBIANQQNLINERTIA = "JacobianQQNLInertia";
const string LNLDS_JACOBIANVELOCITYQNLINERTIA = "JacobianVelocityQNLInertia";

const string LNLDS_M = "M";
const string LNLDS_NDOF = "ndof";
const string LNLDS_MATRIXPLUGIN = "matrixPlugin";
const string LNLDS_VECTORPLUGIN = "vectorPlugin";
//#include "XMLTagsName.h"


class LagrangianDSXML : public DSXML
{
public:
  LagrangianDSXML();
  ~LagrangianDSXML();

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
  inline /*SiconosVector*/ SimpleVector getQ()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->qNode);
  }

  /** \fn void setQ(SiconosVector *v)
  *   \brief allows to save the q of the LagrangianDSXML
  *   \param The q SiconosVector to save
  */
  inline void setQ(SiconosVector *v)
  {
    if (this->hasQ() == false)
    {
      this->qNode = SiconosDOMTreeTools::createVectorNode(this->rootDSXMLNode, LNLDS_Q, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->qNode, v);
  }

  /** \fn SimpleVector getQ0()
  *   \brief Return q0 vector of the LagrangianDSXML
  *   \return SimpleVector : q0 vector of the LagrangianDSXML
  */
  inline SimpleVector getQ0()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->q0Node);
  }

  /** \fn void  setQ0(SiconosVector *v)
  *   \brief allows to save the q0 of the LagrangianDSXML
  *   \param The q0 SiconosVector to save
  */
  inline void  setQ0(SiconosVector *v)
  {
    if (this->q0Node == NULL)
    {
      this->q0Node = SiconosDOMTreeTools::createVectorNode(this->rootDSXMLNode, LNLDS_Q0, v);
    }
    else
    {
      //SiconosDOMTreeTools::setSiconosVectorValue(this->q0Node, *v);
      SiconosDOMTreeTools::setSiconosVectorValue(this->q0Node, v);
    }
  }


  /** \fn SiconosMemoryXML* getQMemoryXML()
  *   \brief Returns the qMemoryXML* of the DSXML
  *   \return SiconosMemoryXML*
  */
  inline SiconosMemoryXML* getQMemoryXML()
  {
    return this->qMemoryXML;
  }

  /** \fn void setQMemory(SiconosMemory* smem)
  *   \brief allows to save the qMemory of the LagrangianDSXML
  *   \param SiconosMemory* : SiconosMemory to save
  */
  inline void setQMemory(SiconosMemory* smem)
  {
    if (this->hasQMemory() == false)
    {
      this->qMemoryXML = new SiconosMemoryXML(NULL, this->rootDSXMLNode, LNLDS_QMEMORY);
      this->qMemoryNode = this->qMemoryXML->getSiconosMemoryXMLNode();

      this->qMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->qMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
    else
    {
      this->qMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->qMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
  }

  /** \fn SimpleVector getVelocity()
  *   \brief Return the velocity of the LagrangianDSXML
  *   \return SimpleVector :  velocity vector of the LagrangianDSXML
  */
  inline SimpleVector getVelocity()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->velocityNode);
  }

  /** \fn void setVelocity(SiconosVector *v)
  *   \brief allows to save the velocity of the LagrangianDSXML
  *   \param The velocity SiconosVector to save
  */
  inline void setVelocity(SiconosVector *v)
  {
    if (this->hasVelocity() == false)
    {
      this->velocityNode = SiconosDOMTreeTools::createVectorNode(this->rootDSXMLNode, LNLDS_VELOCITY, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->velocityNode, v);
    //SiconosDOMTreeTools::setSiconosVectorValue(this->velocityNode, v);
  }

  /** \fn SimpleVector getVelocity0()
  *   \brief Return the initial velocity of the LagrangianDSXML
  *   \return SimpleVector : The velocity0 SiconosVector of the LagrangianDSXML
  */
  inline SimpleVector getVelocity0()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->velocity0Node);
  }

  /** \fn void setVelocity0(SiconosVector *v)
  *   \brief allows to save the velocity0 of the LagrangianDSXML
  *   \param The celocity0 SiconosVector to save
  */
  inline void setVelocity0(SiconosVector *v)
  {
    if (this->velocity0Node == NULL)
    {
      this->velocity0Node = SiconosDOMTreeTools::createVectorNode(this->rootDSXMLNode, LNLDS_VELOCITY0, v);
    }
    else
    {
      //SiconosDOMTreeTools::setSiconosVectorValue(this->velocity0Node, *v);
      SiconosDOMTreeTools::setSiconosVectorValue(this->velocity0Node, v);
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
    return this->velocityMemoryXML;
  }

  /** \fn void setVelocityMemory(SiconosMemory* smem)
  *   \brief allows to save the velocityMemory of the LagrangianDSXML
  *   \param SiconosMemory* : SiconosMemory to save
  */
  inline void setVelocityMemory(SiconosMemory* smem)
  {
    if (this->hasVelocityMemory() == false)
    {
      this->velocityMemoryXML = new SiconosMemoryXML(NULL, this->rootDSXMLNode, LNLDS_VELOCITYMEMORY);
      this->velocityMemoryNode = this->velocityMemoryXML->getSiconosMemoryXMLNode();

      this->velocityMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->velocityMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
    else
    {
      this->velocityMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->velocityMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
  }

  /** \fn inline string getQNLInertiaPlugin()
  *   \brief Return the QNLInertia Plugin name of the LagrangianDSXML
  *   \return The QNLInertia Plugin name of the LagrangianDSXML
  *  \exception XMLException
  */
  inline string getQNLInertiaPlugin()
  {
    if (this->isQNLInertiaPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->QNLInertiaNode, LNLDS_VECTORPLUGIN);
    else
      XMLException::selfThrow("LagrangianDSXML - getQNLInertiaPlugin : QNLInertia is not calculated from a plugin ; QNLInertia vector is given");
  };

  /** \fn SimpleVector getQNLInertiaVector()
  *   \brief Return the QNLInertia vector of the LagrangianDSXML
  *   \return SimpleVector : QNLInertia vector of the LagrangianDSXML
  *  \exception XMLException
  */
  inline SimpleVector getQNLInertiaVector()
  {
    if (this->isQNLInertiaPlugin())
      XMLException::selfThrow("LagrangianDSXML - getQNLInertiaVector : QNLInertia vector is not given ; QNLInertia is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(this->QNLInertiaNode);
  }

  /** \fn void setQNLInertiaPlugin(string plugin)
  *   \brief allows to save the QNLInertia plugin of the LagrangianDSXML
  *   \param string : ths string which contains the name and the location of the plugin
  */
  inline void setQNLInertiaPlugin(string plugin)
  {
    if (this->QNLInertiaNode == NULL)
    {
      this->QNLInertiaNode = SiconosDOMTreeTools::createSingleNode(this->rootDSXMLNode, LNLDS_QNLINERTIA);
      xmlNewProp(this->QNLInertiaNode, (xmlChar*)(LNLDS_VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->QNLInertiaNode, LNLDS_VECTORPLUGIN, plugin);
  }

  /** \fn void setQNLInertiaVector(SiconosVector *v)
  *   \brief allows to save the QNLInertia vector of the LagrangianDSXML
  *   \return The QNLInertia SiconosVector to save
  */
  inline void setQNLInertiaVector(SiconosVector *v)
  {
    if (this->QNLInertiaNode == NULL)
    {
      this->QNLInertiaNode = SiconosDOMTreeTools::createVectorNode(this->rootDSXMLNode, LNLDS_QNLINERTIA, v);
    }
    else
    {
      //SiconosDOMTreeTools::setSiconosVectorValue(this->QNLInertiaNode, *v);
      SiconosDOMTreeTools::setSiconosVectorValue(this->QNLInertiaNode, v);
    }
  }


  /** \fn inline string getFintPlugin()
  *   \brief Return the Fint Plugin name of the LagrangianDSXML
  *   \return The Fint Plugin name of the LagrangianDSXML
  *  \exception XMLException
  */
  inline string getFintPlugin()
  {
    if (this->isFintPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->FintNode, LNLDS_VECTORPLUGIN);
    XMLException::selfThrow("LagrangianDSXML - getFintPlugin : Fint is not calculated from a plugin ; Fint vector is given");
  }

  /** \fn SimpleVector getFintVector()
  *   \brief Return the internal forces vector of the LagrangianDSXML
  *   \return SimpleVector : Fint SiconosVector of the LagrangianDSXML
  *   \exception XMLException
  */
  inline SimpleVector getFintVector()
  {
    if (this->isFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFintVector : Fint vector is not given ; Fint is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(this->FintNode);
  }

  /** \fn void setFintVector(SiconosVector *v)
  *   \brief allows to save the Fint vector of the LagrangianDSXML
  *   \param The Fint SiconosVector to save
  */
  inline void setFintVector(SiconosVector *v)
  {
    if (this->hasFint())
    {
      this->FintNode = SiconosDOMTreeTools::createVectorNode(this->rootDSXMLNode, LNLDS_FINT, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->FintNode, v);
  }

  /** \fn void setFextPlugin(string plugin)
  *   \brief allows to save the Fext plugin of the LagrangianDSXML
  *   \param string : ths string which contains the name and the location of the plugin
  */
  inline void setFintPlugin(string plugin)
  {
    if (this->FintNode == NULL)
    {
      this->FintNode = SiconosDOMTreeTools::createSingleNode(this->rootDSXMLNode, LNLDS_FINT);
      xmlNewProp(this->FintNode, (xmlChar*)(LNLDS_VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->FintNode, LNLDS_VECTORPLUGIN, plugin);
  }

  /** \fn inline string getFextPlugin()
  *   \brief Return the Fext Plugin name of the LagrangianDSXML
  *   \return The Fext Plugin name of the LagrangianDSXML
  *  \exception XMLException
  */
  inline string getFextPlugin()
  {
    if (this->isFextPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->FextNode, LNLDS_VECTORPLUGIN);
    XMLException::selfThrow("LagrangianDSXML - getFextPlugin : Fext is not calculated from a plugin ; Fext vector is given");
  }

  /** \fn SimpleVector getFextVector()
  *   \brief Return the external forces vector of the LagrangianDSXML
  *   \return SimpleVector : Fext vector of the LagrangianDSXML
  *  \exception XMLException
  */
  inline SimpleVector getFextVector()
  {
    if (this->isFextPlugin())
      XMLException::selfThrow("LagrangianDSXML - getFextVector : Fext matrix is not given ; Fext is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(this->FextNode);
  }

  /** \fn void setFextVector(SiconosVector *v)
  *   \brief allows to save the Fint vector of the LagrangianDSXML
  *   \param The Fint SiconosVector to save
  */
  inline void setFextVector(SiconosVector *v)
  {
    //SiconosDOMTreeTools::setSiconosVectorValue(this->FextNode, *v);
    SiconosDOMTreeTools::setSiconosVectorValue(this->FextNode, v);
  }

  /** \fn void setFextPlugin(string plugin)
  *   \brief allows to save the Fext plugin of the LagrangianDSXML
  *   \param string : ths string which contains the name and the location of the plugin
  */
  inline void setFextPlugin(string plugin)
  {
    if (this->FextNode == NULL)
    {
      this->FextNode = SiconosDOMTreeTools::createSingleNode(this->rootDSXMLNode, LNLDS_FEXT);
      xmlNewProp(this->FextNode, (xmlChar*)(LNLDS_VECTORPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else
    {
      SiconosDOMTreeTools::setStringAttributeValue(this->FextNode, LNLDS_VECTORPLUGIN, plugin);
    }
  }

  /** \fn inline string getJacobianQFintPlugin()
  *   \brief Return the JacobianQFint Plugin name of the LagrangianDSXML
  *   \return The JacobianQFint Plugin name of the LagrangianDSXML
  *  \exception XMLException
  */
  inline string getJacobianQFintPlugin()
  {
    if (this->isJacobianQFintPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->jacobianQFintNode, LNLDS_MATRIXPLUGIN);
    XMLException::selfThrow("LagrangianDSXML - getJacobianQFintPlugin : JacobianQFint is not calculated from a plugin ; JacobianQFint matrix is given");
  }

  /** \fn SiconosMatrix getJacobianQFintMatrix()
  *   \brief Return the JacobianQFint matrix of the LagrangianDSXML
  *   \return The JacobianQFint SiconosMatrix of the LagrangianDSXML
  *  \exception XMLException
  */
  inline SiconosMatrix getJacobianQFintMatrix()
  {
    if (this->isJacobianQFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianQFintMatrix : JacobianQFint matrix is not given ; JacobianQFint is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->jacobianQFintNode);
  }

  /** \fn void setJacobianQFintPlugin(string plugin)
  *   \brief allows to save the jacobianQFint plugin of the LagrangianDSXML
  *   \param string : ths string which contains the name and the location of the plugin
  */
  inline void setJacobianQFintPlugin(string plugin)
  {
    if (this->jacobianQFintNode == NULL)
    {
      this->jacobianQFintNode = SiconosDOMTreeTools::createSingleNode(this->rootDSXMLNode, LNLDS_JACOBIANQFINT);
      xmlNewProp(this->jacobianQFintNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->jacobianQFintNode, LNLDS_JACOBIANQFINT, plugin);
  }

  /** \fn void setJacobianQFintMatrix(SiconosMatrix *m)
  *   \brief allows to save the JacobianQFint matrix of the LagrangianDSXML
  *   \return The JacobianQFint SiconosMatrix to save
  */
  inline void setJacobianQFintMatrix(SiconosMatrix *m)
  {
    SiconosDOMTreeTools::setSiconosMatrixValue(this->jacobianQFintNode, m);
  }

  /** \fn inline string getJacobianVelocityFintPlugin()
  *   \brief Return the JacobianVelocityFint Plugin name of the LagrangianDSXML
  *   \return The JacobianVelocityFint Plugin name of the LagrangianDSXML
  *  \exception XMLException
  */
  inline string getJacobianVelocityFintPlugin()
  {
    if (this->isJacobianVelocityFintPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->jacobianVelocityFintNode, LNLDS_MATRIXPLUGIN);
    XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityFintPlugin : JacobianVelocityFint is not calculated from a plugin ; JacobianVelocityFint matrix is given");
  }

  /** \fn SiconosMatrix getJacobianVelocityFintMatrix()
  *   \brief Return the JacobianVelocityFint matrix of the LagrangianDSXML
  *   \return The JacobianVelocityFint SiconosMatrix of the LagrangianDSXML
  *  \exception XMLException
  */
  inline SiconosMatrix getJacobianVelocityFintMatrix()
  {
    if (this->isJacobianVelocityFintPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityFintMatrix : JacobianVelocityFint matrix is not given ; JacobianVelocityFint is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->jacobianVelocityFintNode);
  }

  /** \fn void setJacobianVelocityFintPlugin(string plugin)
  *   \brief allows to save the jacobianVelocityFint plugin of the LagrangianDSXML
  *   \param string : ths string which contains the name and the location of the plugin
  */
  inline void setJacobianVelocityFintPlugin(string plugin)
  {
    if (this->jacobianVelocityFintNode == NULL)
    {
      this->jacobianVelocityFintNode = SiconosDOMTreeTools::createSingleNode(this->rootDSXMLNode, LNLDS_JACOBIANVELOCITYFINT);
      xmlNewProp(this->jacobianVelocityFintNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->jacobianVelocityFintNode, LNLDS_JACOBIANVELOCITYFINT, plugin);
  }

  /** \fn void setJacobianVelocityFintMatrix(SiconosMatrix *m)
  *   \brief allows to save the JacobianVelocityFint matrix of the LagrangianDSXML
  *   \return The JacobianVelocityFint SiconosMatrix to save
  */
  inline void setJacobianVelocityFintMatrix(SiconosMatrix *m)
  {
    //SiconosDOMTreeTools::setSiconosMatrixValue(this->jacobianVelocityFintNode, *m);
    SiconosDOMTreeTools::setSiconosMatrixValue(this->jacobianVelocityFintNode, m);
  }

  /** \fn inline string getJacobianQQPlugin()
  *   \brief Return the JacobianQQ Plugin name of the LagrangianDSXML
  *   \return The JacobianQQ Plugin name of the LagrangianDSXML
  *  \exception XMLException
  */
  inline string getJacobianQQNLInertiaPlugin()
  {
    if (this->isJacobianQQNLInertiaPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->jacobianQQNLInertiaNode, LNLDS_MATRIXPLUGIN);
    XMLException::selfThrow("LagrangianDSXML - getJacobianQQNLInertiaPlugin : JacobianQQNLInertia is not calculated from a plugin ; JacobianQQNLInertia matrix is given");
  }

  /** \fn SiconosMatrix getJacobianQQMatrix()
  *   \brief Return the JacobianQQ matrix of the LagrangianDSXML
  *   \return The JacobianQQ SiconosMatrix of the LagrangianDSXML
  *  \exception XMLException
  */
  inline SiconosMatrix getJacobianQQNLInertiaMatrix()
  {
    if (this->isJacobianQQNLInertiaPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianQQNLInertiaMatrix : JacobianQQNLInertia matrix is not given ; JacobianQQNLInertia is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->jacobianQQNLInertiaNode);
  }

  /** \fn void setJacobianQQNLInertiaPlugin(string plugin)
  *   \brief allows to save the jacobianQQNLInertia plugin of the LagrangianDSXML
  *   \param string : ths string which contains the name and the location of the plugin
  */
  inline void setJacobianQQNLInertiaPlugin(string plugin)
  {
    if (this->jacobianQQNLInertiaNode == NULL)
    {
      this->jacobianQQNLInertiaNode = SiconosDOMTreeTools::createSingleNode(this->rootDSXMLNode, LNLDS_JACOBIANQQNLINERTIA);
      xmlNewProp(this->jacobianQQNLInertiaNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->jacobianQQNLInertiaNode, LNLDS_JACOBIANQQNLINERTIA, plugin);
  }

  /** \fn void setJacobianQQMatrix(SiconosMatrix *m)
  *   \brief allows to save the JacobianQQ matrix of the LagrangianDSXML
  *   \return The JacobianQQ SiconosMatrix to save
  */
  inline void setJacobianQQNLInertiaMatrix(SiconosMatrix *m)
  {
    //SiconosDOMTreeTools::setSiconosMatrixValue(this->jacobianQQNLInertiaNode, *m);
    SiconosDOMTreeTools::setSiconosMatrixValue(this->jacobianQQNLInertiaNode, m);
  }

  /** \fn inline string getJacobianVelocityQNLInertiaPlugin()
  *   \brief Return the JacobianVelocityQNLInertia Plugin name of the LagrangianDSXML
  *   \return The JacobianVelocityQNLInertia Plugin name of the LagrangianDSXML
  *  \exception XMLException
  */
  inline string getJacobianVelocityQNLInertiaPlugin()
  {
    if (this->isJacobianVelocityQNLInertiaPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->jacobianVelocityQNLInertiaNode, LNLDS_MATRIXPLUGIN);
    XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityQNLInertiaPlugin : JacobianVelocityQNLInertia is not calculated from a plugin ; JacobianVelocityQNLInertia matrix is given");
  }

  /** \fn SiconosMatrix getJacobianVelocityQNLInertiaMatrix()
  *   \brief Return the JacobianVelocityQNLInertia matrix of the LagrangianDSXML
  *   \return The JacobianVelocityQNLInertia SiconosMatrix of the LagrangianDSXML
  *  \exception XMLException
  */
  inline SiconosMatrix getJacobianVelocityQNLInertiaMatrix()
  {
    if (this->isJacobianVelocityQNLInertiaPlugin())
      XMLException::selfThrow("LagrangianDSXML - getJacobianVelocityQNLInertiaMatrix : JacobianVelocityQNLInertia matrix is not given ; JacobianVelocityQNLInertia is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->jacobianVelocityQNLInertiaNode);
  }

  /** \fn void setJacobianVelocityQNLInertiaPlugin(string plugin)
  *   \brief allows to save the jacobianVelocityQNLInertiaPlugin plugin of the LagrangianDSXML
  *   \param string : ths string which contains the name and the location of the plugin
  */
  inline void setJacobianVelocityQNLInertiaPlugin(string plugin)
  {
    if (this->jacobianVelocityQNLInertiaNode == NULL)
    {
      this->jacobianVelocityQNLInertiaNode = SiconosDOMTreeTools::createSingleNode(this->rootDSXMLNode, LNLDS_JACOBIANVELOCITYQNLINERTIA);
      xmlNewProp(this->jacobianVelocityQNLInertiaNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->jacobianVelocityQNLInertiaNode, LNLDS_JACOBIANVELOCITYQNLINERTIA, plugin);
  }

  /** \fn void setJacobianVelocityQNLInertiaMatrix(SiconosMatrix *m)
  *   \brief allows to save the JacobianVelocityQNLInertia matrix of the LagrangianDSXML
  *   \return The JacobianVelocityQNLInertia SiconosMatrix to save
  */
  inline void setJacobianVelocityQNLInertiaMatrix(SiconosMatrix *m)
  {
    SiconosDOMTreeTools::setSiconosMatrixValue(this->jacobianVelocityQNLInertiaNode, m);
  }

  /** \fn inline string getMPlugin()
  *   \brief Return the M Plugin name of the LagrangianDSXML
  *   \return The M Plugin name of the LagrangianDSXML
  *  \exception XMLException
  */
  inline string getMPlugin()
  {
    if (this->isMPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->MNode, LNLDS_MATRIXPLUGIN);
    XMLException::selfThrow("LagrangianDSXML - getMPlugin : M is not calculated from a plugin ; M matrix is given");
  }

  /** \fn SiconosMatrix getMMatrix()
  *   \brief Return the M matrix of the LagrangianDSXML
  *   \return The M SiconosMatrix of the LagrangianDSXML
  *  \exception XMLException
  */
  inline SiconosMatrix getMMatrix()
  {
    if (this->isMPlugin())
      XMLException::selfThrow("LagrangianDSXML - getMMatrix : M matrix is not given ; M is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->MNode);
  }

  /** \fn void setMPlugin(string plugin)
  *   \brief allows to save the jacobianVelocityQNLInertiaPlugin plugin of the LagrangianDSXML
  *   \param string : ths string which contains the name and the location of the plugin
  */
  inline void setMPlugin(string plugin)
  {
    if (this->MNode == NULL)
    {
      this->MNode = SiconosDOMTreeTools::createSingleNode(this->rootDSXMLNode, LNLDS_M);
      xmlNewProp(this->MNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str());
    }
    //      else if (xmlHasProp((xmlNodePtr)this->MNode, (xmlChar *)LNLDS_MATRIXPLUGIN.c_str()))
    //      {
    //        xmlSetProp((xmlNode *) this->MNode, (xmlChar *)LNLDS_MATRIXPLUGIN.c_str(), (xmlChar *) (plugin.c_str()));
    //      }
    else
    {
      //cout<<"MNode == "<< this->MNode <<endl;
      //cout<<"MNode name == "<< this->MNode->name <<endl;
      //cout<<"MNode content == "<< this->MNode->content <<endl;
      SiconosDOMTreeTools::setStringAttributeValue(this->MNode, LNLDS_MATRIXPLUGIN, plugin);
      //        xmlNewProp(this->MNode, (xmlChar*)(LNLDS_MATRIXPLUGIN.c_str()), (xmlChar*)plugin.c_str() );
    }
  }

  /** \fn void setMMatrix(SiconosMatrix *m)
  *   \brief allows to save the M matrix of the LagrangianDSXML
  *   \return The M SiconosMatrix to save
  */
  inline void setMMatrix(SiconosMatrix *m)
  {
    if (this->MNode == NULL)
    {
      this->MNode = SiconosDOMTreeTools::createMatrixNode(this->rootDSXMLNode, LNLDS_M, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->MNode, m);
  }


  /** \fn int getNdof()
  *   \brief Return the ndof for the LagrangianDSXML
  *   \return The ndof integer for the LagrangianDSXML
  */
  inline int getNdof()
  {
    return  SiconosDOMTreeTools::getIntegerContentValue(this->ndofNode);
  }

  /** \fn void setNdof(int i)
  *   \brief allows to save the ndof for the LagrangianDSXML
  *   \return The ndof integer to save
  */
  inline void setNdof(int i)
  {
    if (this->ndofNode == NULL)
    {
      this->ndofNode = SiconosDOMTreeTools::createIntegerNode(this->rootDSXMLNode, LNLDS_NDOF, i);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(this->ndofNode, i);
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
    return (this->MNode != NULL);
  }

  /** \fn bool hasFint()
   *  \brief determines if Fint is defined in the DOM tree
   *  \return bool : true if Fint is defined, false otherwise
   */
  inline bool hasFint()
  {
    return (this->FintNode != NULL);
  }

  /** \fn bool hasFext()
   *  \brief determines if Fext is defined in the DOM tree
   *  \return bool : true if Fext is defined, false otherwise
   */
  inline bool hasFext()
  {
    return (this->FextNode != NULL);
  }

  /** \fn bool hasJacobianQFint()
   *  \brief determines if jacobianQFint is defined in the DOM tree
   *  \return bool : true if jacobianQFint is defined, false otherwise
   */
  inline bool hasJacobianQFint()
  {
    return (this->jacobianQFintNode != NULL);
  }

  /** \fn bool hasJacobianVelocityFint()
   *  \brief determines if jacobianVelocityFint is defined in the DOM tree
   *  \return bool : true if jacobianVelocityFint is defined, false otherwise
   */
  inline bool hasJacobianVelocityFint()
  {
    return (this->jacobianVelocityFintNode != NULL);
  }

  /** \fn bool hasJacobianQQNLInertia()
   *  \brief determines if jacobianQQNLInertia is defined in the DOM tree
   *  \return bool : true if jacobianQQNLInertia is defined, false otherwise
   */
  inline bool hasJacobianQQNLInertia()
  {
    return (this->jacobianQQNLInertiaNode != NULL);
  }

  /** \fn bool hasJacobianVelocityQNLInertia()
   *  \brief determines if jacobianVelocityQNLInertia is defined in the DOM tree
   *  \return bool : true if jacobianVelocityQNLInertia is defined, false otherwise
   */
  inline bool hasJacobianVelocityQNLInertia()
  {
    return (this->jacobianVelocityQNLInertiaNode != NULL);
  }

  /** \fn bool hasQNLInertia()
   *  \brief determines if QNLInertia is defined in the DOM tree
   *  \return bool : true if QNLInertia is defined, false otherwise
   */
  inline bool hasQNLInertia()
  {
    return (this->QNLInertiaNode != NULL);
  }

  /** \fn bool hasQMemory()
   *  \brief returns true if qMemoryNode is defined
   *  \return true if qMemoryNode is defined
   */
  inline bool hasQMemory()
  {
    return (this->qMemoryNode != NULL);
  }

  /** \fn bool hasVelocityMemory()
   *  \brief returns true if velocityMemoryNode is defined
   *  \return true if velocityMemoryNode is defined
   */
  inline bool hasVelocityMemory()
  {
    return (this->velocityMemoryNode != NULL);
  }

  /** \fn bool hasQ()
   *  \brief returns true if qNode is defined
   *  \return true if qNode is defined
   */
  inline bool hasQ()
  {
    return (this->qNode != NULL);
  }

  /** \fn bool hasVelocity()
   *  \brief returns true if velocityNode is defined
   *  \return true if velocityNode is defined
   */
  inline bool hasVelocity()
  {
    return (this->velocityNode != NULL);
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
