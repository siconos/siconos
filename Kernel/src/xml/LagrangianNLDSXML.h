
/** \class LagrangianNLDSXML
*   \brief This class manages Lagrangian NLDS data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/11/2004
*
*
*
* LagrangianNLDSXML allows to manage data of a LagrangianNLDS DOM tree.
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


using namespace std;


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


class LagrangianNLDSXML : public DSXML
{
public:
  LagrangianNLDSXML();
  ~LagrangianNLDSXML();

  /** \fn LagrangianNLDSXML(xmlNode * LagrangianNLDSNode, bool isBVP)
  *   \brief Build a LagrangianNLDSXML object from a DOM tree describing a DS
  *   \param xmlNode * LagrangianNLDSNode : the LagrangianNLDS DOM tree
  *   \param bool isBVP : if NSDS is BVP LagrangianNLDS have boundary condition
  */
  LagrangianNLDSXML(xmlNode * LagrangianNLDSNode, bool isBVP);

  /** \fn SimpleVector getQ()
  *   \brief Return  q vector of the LagrangianNLDSXML
  *   \return SimpleVector : q vector of the LagrangianNLDSXML
  */
  inline /*SiconosVector*/ SimpleVector getQ()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->qNode);
  }

  /** \fn void setQ(SiconosVector *v)
  *   \brief allows to save the q of the LagrangianNLDSXML
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
  *   \brief Return q0 vector of the LagrangianNLDSXML
  *   \return SimpleVector : q0 vector of the LagrangianNLDSXML
  */
  inline SimpleVector getQ0()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->q0Node);
  }

  /** \fn void  setQ0(SiconosVector *v)
  *   \brief allows to save the q0 of the LagrangianNLDSXML
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
  *   \brief allows to save the qMemory of the LagrangianNLDSXML
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
  *   \brief Return the velocity of the LagrangianNLDSXML
  *   \return SimpleVector :  velocity vector of the LagrangianNLDSXML
  */
  inline SimpleVector getVelocity()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->velocityNode);
  }

  /** \fn void setVelocity(SiconosVector *v)
  *   \brief allows to save the velocity of the LagrangianNLDSXML
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
  *   \brief Return the initial velocity of the LagrangianNLDSXML
  *   \return SimpleVector : The velocity0 SiconosVector of the LagrangianNLDSXML
  */
  inline SimpleVector getVelocity0()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->velocity0Node);
  }

  /** \fn void setVelocity0(SiconosVector *v)
  *   \brief allows to save the velocity0 of the LagrangianNLDSXML
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
  //    *   \brief Return the velocityMemory of the LagrangianNLDSXML
  //    *   \return SiconosMemory velocityMemory of the LagrangianNLDSXML
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
  *   \brief allows to save the velocityMemory of the LagrangianNLDSXML
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
  *   \brief Return the QNLInertia Plugin name of the LagrangianNLDSXML
  *   \return The QNLInertia Plugin name of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline string getQNLInertiaPlugin()
  {
    if (this->isQNLInertiaPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->QNLInertiaNode, LNLDS_VECTORPLUGIN);
    else
      XMLException::selfThrow("LagrangianNLDSXML - getQNLInertiaPlugin : QNLInertia is not calculated from a plugin ; QNLInertia vector is given");
  };

  /** \fn SimpleVector getQNLInertiaVector()
  *   \brief Return the QNLInertia vector of the LagrangianNLDSXML
  *   \return SimpleVector : QNLInertia vector of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline SimpleVector getQNLInertiaVector()
  {
    if (this->isQNLInertiaPlugin())
      XMLException::selfThrow("LagrangianNLDSXML - getQNLInertiaVector : QNLInertia vector is not given ; QNLInertia is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(this->QNLInertiaNode);
  }

  /** \fn void setQNLInertiaPlugin(string plugin)
  *   \brief allows to save the QNLInertia plugin of the LagrangianNLDSXML
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
  *   \brief allows to save the QNLInertia vector of the LagrangianNLDSXML
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
  *   \brief Return the Fint Plugin name of the LagrangianNLDSXML
  *   \return The Fint Plugin name of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline string getFintPlugin()
  {
    if (this->isFintPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->FintNode, LNLDS_VECTORPLUGIN);
    XMLException::selfThrow("LagrangianNLDSXML - getFintPlugin : Fint is not calculated from a plugin ; Fint vector is given");
  }

  /** \fn SimpleVector getFintVector()
  *   \brief Return the internal forces vector of the LagrangianNLDSXML
  *   \return SimpleVector : Fint SiconosVector of the LagrangianNLDSXML
  *   \exception XMLException
  */
  inline SimpleVector getFintVector()
  {
    if (this->isFintPlugin())
      XMLException::selfThrow("LagrangianNLDSXML - getFintVector : Fint vector is not given ; Fint is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(this->FintNode);
  }

  /** \fn void setFintVector(SiconosVector *v)
  *   \brief allows to save the Fint vector of the LagrangianNLDSXML
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
  *   \brief allows to save the Fext plugin of the LagrangianNLDSXML
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
  *   \brief Return the Fext Plugin name of the LagrangianNLDSXML
  *   \return The Fext Plugin name of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline string getFextPlugin()
  {
    if (this->isFextPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->FextNode, LNLDS_VECTORPLUGIN);
    XMLException::selfThrow("LagrangianNLDSXML - getFextPlugin : Fext is not calculated from a plugin ; Fext vector is given");
  }

  /** \fn SimpleVector getFextVector()
  *   \brief Return the external forces vector of the LagrangianNLDSXML
  *   \return SimpleVector : Fext vector of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline SimpleVector getFextVector()
  {
    if (this->isFextPlugin())
      XMLException::selfThrow("LagrangianNLDSXML - getFextVector : Fext matrix is not given ; Fext is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(this->FextNode);
  }

  /** \fn void setFextVector(SiconosVector *v)
  *   \brief allows to save the Fint vector of the LagrangianNLDSXML
  *   \param The Fint SiconosVector to save
  */
  inline void setFextVector(SiconosVector *v)
  {
    //SiconosDOMTreeTools::setSiconosVectorValue(this->FextNode, *v);
    SiconosDOMTreeTools::setSiconosVectorValue(this->FextNode, v);
  }

  /** \fn void setFextPlugin(string plugin)
  *   \brief allows to save the Fext plugin of the LagrangianNLDSXML
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
  *   \brief Return the JacobianQFint Plugin name of the LagrangianNLDSXML
  *   \return The JacobianQFint Plugin name of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline string getJacobianQFintPlugin()
  {
    if (this->isJacobianQFintPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->jacobianQFintNode, LNLDS_MATRIXPLUGIN);
    XMLException::selfThrow("LagrangianNLDSXML - getJacobianQFintPlugin : JacobianQFint is not calculated from a plugin ; JacobianQFint matrix is given");
  }

  /** \fn SiconosMatrix getJacobianQFintMatrix()
  *   \brief Return the JacobianQFint matrix of the LagrangianNLDSXML
  *   \return The JacobianQFint SiconosMatrix of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline SiconosMatrix getJacobianQFintMatrix()
  {
    if (this->isJacobianQFintPlugin())
      XMLException::selfThrow("LagrangianNLDSXML - getJacobianQFintMatrix : JacobianQFint matrix is not given ; JacobianQFint is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->jacobianQFintNode);
  }

  /** \fn void setJacobianQFintPlugin(string plugin)
  *   \brief allows to save the jacobianQFint plugin of the LagrangianNLDSXML
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
  *   \brief allows to save the JacobianQFint matrix of the LagrangianNLDSXML
  *   \return The JacobianQFint SiconosMatrix to save
  */
  inline void setJacobianQFintMatrix(SiconosMatrix *m)
  {
    SiconosDOMTreeTools::setSiconosMatrixValue(this->jacobianQFintNode, m);
  }

  /** \fn inline string getJacobianVelocityFintPlugin()
  *   \brief Return the JacobianVelocityFint Plugin name of the LagrangianNLDSXML
  *   \return The JacobianVelocityFint Plugin name of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline string getJacobianVelocityFintPlugin()
  {
    if (this->isJacobianVelocityFintPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->jacobianVelocityFintNode, LNLDS_MATRIXPLUGIN);
    XMLException::selfThrow("LagrangianNLDSXML - getJacobianVelocityFintPlugin : JacobianVelocityFint is not calculated from a plugin ; JacobianVelocityFint matrix is given");
  }

  /** \fn SiconosMatrix getJacobianVelocityFintMatrix()
  *   \brief Return the JacobianVelocityFint matrix of the LagrangianNLDSXML
  *   \return The JacobianVelocityFint SiconosMatrix of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline SiconosMatrix getJacobianVelocityFintMatrix()
  {
    if (this->isJacobianVelocityFintPlugin())
      XMLException::selfThrow("LagrangianNLDSXML - getJacobianVelocityFintMatrix : JacobianVelocityFint matrix is not given ; JacobianVelocityFint is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->jacobianVelocityFintNode);
  }

  /** \fn void setJacobianVelocityFintPlugin(string plugin)
  *   \brief allows to save the jacobianVelocityFint plugin of the LagrangianNLDSXML
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
  *   \brief allows to save the JacobianVelocityFint matrix of the LagrangianNLDSXML
  *   \return The JacobianVelocityFint SiconosMatrix to save
  */
  inline void setJacobianVelocityFintMatrix(SiconosMatrix *m)
  {
    //SiconosDOMTreeTools::setSiconosMatrixValue(this->jacobianVelocityFintNode, *m);
    SiconosDOMTreeTools::setSiconosMatrixValue(this->jacobianVelocityFintNode, m);
  }

  /** \fn inline string getJacobianQQPlugin()
  *   \brief Return the JacobianQQ Plugin name of the LagrangianNLDSXML
  *   \return The JacobianQQ Plugin name of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline string getJacobianQQNLInertiaPlugin()
  {
    if (this->isJacobianQQNLInertiaPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->jacobianQQNLInertiaNode, LNLDS_MATRIXPLUGIN);
    XMLException::selfThrow("LagrangianNLDSXML - getJacobianQQNLInertiaPlugin : JacobianQQNLInertia is not calculated from a plugin ; JacobianQQNLInertia matrix is given");
  }

  /** \fn SiconosMatrix getJacobianQQMatrix()
  *   \brief Return the JacobianQQ matrix of the LagrangianNLDSXML
  *   \return The JacobianQQ SiconosMatrix of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline SiconosMatrix getJacobianQQNLInertiaMatrix()
  {
    if (this->isJacobianQQNLInertiaPlugin())
      XMLException::selfThrow("LagrangianNLDSXML - getJacobianQQNLInertiaMatrix : JacobianQQNLInertia matrix is not given ; JacobianQQNLInertia is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->jacobianQQNLInertiaNode);
  }

  /** \fn void setJacobianQQNLInertiaPlugin(string plugin)
  *   \brief allows to save the jacobianQQNLInertia plugin of the LagrangianNLDSXML
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
  *   \brief allows to save the JacobianQQ matrix of the LagrangianNLDSXML
  *   \return The JacobianQQ SiconosMatrix to save
  */
  inline void setJacobianQQNLInertiaMatrix(SiconosMatrix *m)
  {
    //SiconosDOMTreeTools::setSiconosMatrixValue(this->jacobianQQNLInertiaNode, *m);
    SiconosDOMTreeTools::setSiconosMatrixValue(this->jacobianQQNLInertiaNode, m);
  }

  /** \fn inline string getJacobianVelocityQNLInertiaPlugin()
  *   \brief Return the JacobianVelocityQNLInertia Plugin name of the LagrangianNLDSXML
  *   \return The JacobianVelocityQNLInertia Plugin name of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline string getJacobianVelocityQNLInertiaPlugin()
  {
    if (this->isJacobianVelocityQNLInertiaPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->jacobianVelocityQNLInertiaNode, LNLDS_MATRIXPLUGIN);
    XMLException::selfThrow("LagrangianNLDSXML - getJacobianVelocityQNLInertiaPlugin : JacobianVelocityQNLInertia is not calculated from a plugin ; JacobianVelocityQNLInertia matrix is given");
  }

  /** \fn SiconosMatrix getJacobianVelocityQNLInertiaMatrix()
  *   \brief Return the JacobianVelocityQNLInertia matrix of the LagrangianNLDSXML
  *   \return The JacobianVelocityQNLInertia SiconosMatrix of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline SiconosMatrix getJacobianVelocityQNLInertiaMatrix()
  {
    if (this->isJacobianVelocityQNLInertiaPlugin())
      XMLException::selfThrow("LagrangianNLDSXML - getJacobianVelocityQNLInertiaMatrix : JacobianVelocityQNLInertia matrix is not given ; JacobianVelocityQNLInertia is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->jacobianVelocityQNLInertiaNode);
  }

  /** \fn void setJacobianVelocityQNLInertiaPlugin(string plugin)
  *   \brief allows to save the jacobianVelocityQNLInertiaPlugin plugin of the LagrangianNLDSXML
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
  *   \brief allows to save the JacobianVelocityQNLInertia matrix of the LagrangianNLDSXML
  *   \return The JacobianVelocityQNLInertia SiconosMatrix to save
  */
  inline void setJacobianVelocityQNLInertiaMatrix(SiconosMatrix *m)
  {
    SiconosDOMTreeTools::setSiconosMatrixValue(this->jacobianVelocityQNLInertiaNode, m);
  }

  /** \fn inline string getMPlugin()
  *   \brief Return the M Plugin name of the LagrangianNLDSXML
  *   \return The M Plugin name of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline string getMPlugin()
  {
    if (this->isMPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->MNode, LNLDS_MATRIXPLUGIN);
    XMLException::selfThrow("LagrangianNLDSXML - getMPlugin : M is not calculated from a plugin ; M matrix is given");
  }

  /** \fn SiconosMatrix getMMatrix()
  *   \brief Return the M matrix of the LagrangianNLDSXML
  *   \return The M SiconosMatrix of the LagrangianNLDSXML
  *  \exception XMLException
  */
  inline SiconosMatrix getMMatrix()
  {
    if (this->isMPlugin())
      XMLException::selfThrow("LagrangianNLDSXML - getMMatrix : M matrix is not given ; M is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->MNode);
  }

  /** \fn void setMPlugin(string plugin)
  *   \brief allows to save the jacobianVelocityQNLInertiaPlugin plugin of the LagrangianNLDSXML
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
  *   \brief allows to save the M matrix of the LagrangianNLDSXML
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
  *   \brief Return the ndof for the LagrangianNLDSXML
  *   \return The ndof integer for the LagrangianNLDSXML
  */
  inline int getNdof()
  {
    return  SiconosDOMTreeTools::getIntegerContentValue(this->ndofNode);
  }

  /** \fn void setNdof(int i)
  *   \brief allows to save the ndof for the LagrangianNLDSXML
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

  /** \fn loadLagrangianNLDSProperties()
  *   \brief load the different properties of a LagrangianNLDS
  *   \exception XMLException : if a property of the LagrangianNLDS lacks in the DOM tree
  */
  void loadLagrangianNLDSProperties();

};

#endif
//$Log: LagrangianNLDSXML.h,v $
//Revision 1.42  2005/03/08 12:41:38  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.41  2005/02/24 15:50:21  jbarbier
//- LCP prepared to changes needed for several interactions
//
//- new function for the SiconosMatrices to copy a block matrix into another matrix
//
//- tag names of BoundaryConditionXML, DSInputOutputXML, DSXML, InteractionXML, LagrangianLinearRXML, LagrangianNLDSXML put in XMLTagNames.h
//
//Revision 1.40  2005/01/18 10:35:17  jbarbier
//- attribute "r" no longer used for Moreau integrator
//
//- modificatoin in the tests for Moreau integrator
//
//- file XMLTagsName.h for further use to regroup all xml tags name...
//
//Revision 1.39  2004/09/27 08:24:26  charlety
//
//_ Modifications in doxygen comments.
//
//Revision 1.38  2004/09/22 10:54:44  jbarbier
//- light modification according to the attribute mass of the lagrangian dynamical
//systems. The lagrangianNLDS take always an function from a plugin to compute the
//mass, whereas the lagrangianTIDS needs only a matrix.
//
//- xml input files have been modified in consequence
//
//Revision 1.37  2004/09/10 11:26:27  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.36  2004/08/23 14:30:02  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.35  2004/08/20 15:26:45  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NSDS and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.34  2004/08/10 12:04:30  jbarbier
//- save of the plugin's name for fInt
//
//Revision 1.33  2004/08/09 15:00:55  jbarbier
//- changes in the cardinality of some attributes of the DynamicalSystem,
//OneStepIntegrator
//
//- modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//
//- corrections in the test xml files
//
//Revision 1.32  2004/08/06 11:27:53  jbarbier
//- new tests with the XML and the optional attributes
//
//- tests on the save of the XML data
//
//Revision 1.31  2004/08/04 11:03:23  jbarbier
//- about the SiconosMemory : when a SiconosMemory has a maxSize greater than the
//number of steps in memory required by an integrator, the oldest SiconosVector
//are deleted
//
//- the way to initialize the SiconosMemory by the integrator has been updated to
//match with these changes
//
//Revision 1.30  2004/07/30 14:37:15  jbarbier
//- saving methods for DynamicalSystemXML and LagrangianNLDSXML
//
//Revision 1.29  2004/07/29 14:25:42  jbarbier
