/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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

/** \class DynamicalSystemXML
 *   \brief This class manages DynamicalSystem data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.1.2.
 *   \date 04/04/2004
 *
 *
 * DynamicalSystemXML allows to manage data of a DynamicalSystem DOM tree.
 */
#ifndef __DynamicalSystemXML__
#define __DynamicalSystemXML__

#include "SiconosDOMTreeTools.h"
#include "DynamicalSystem.h"

#include "SiconosMemory.h"
#include "SiconosMemoryXML.h"
#include "BoundaryConditionXML.h"
#include "DSInputOutputXML.h"
#include "BoundaryCondition.h"

class DynamicalSystem;
class BoundaryCondition;
class BoundaryConditionXML;
class DSInputOutputXML;
class SiconosMemoryXML;

//Tags
const std::string DS_N = "n";
const std::string DS_X0 = "x0";
const std::string DS_X = "x";
const std::string DS_XDOT = "xDot";
const std::string DS_R = "R";
const std::string DS_U = "u";
const std::string DS_T = "T";
const std::string DS_XMEMORY = "xMemory";
const std::string DS_XDOTMEMORY = "xDotMemory";
const std::string DS_RMEMORY = "RMemory";
const std::string DS_STEPSINMEMORY = "StepsInMemory";
const std::string DS_VECTORFIELD = "vectorField";
const std::string DS_COMPUTEJACOBIANX = "computeJacobianX";
const std::string DS_MATRIXPLUGIN = "matrixPlugin";
const std::string DS_VECTORPLUGIN = "vectorPlugin";

class DynamicalSystemXML
{
public:

  DynamicalSystemXML();

  virtual ~DynamicalSystemXML();

  /** \fn DynamicalSystemXML(xmlNode * DynamicalSystemNode)
   *   \brief Build a DynamicalSystemXML object from a DOM tree describing a DynamicalSystem
   *   \param xmlNode * DynamicalSystemNode : the DynamicalSystem DOM tree
   *   \param bool isBVP : if NonSmoothDynamicalSystem is IBP DynamicalSystem have boundary condition
   */
  DynamicalSystemXML(xmlNode * DynamicalSystemNode, const bool&);

  /** \fn int getNumber()
   *   \brief Return the number of the DynamicalSystemXML
   *   \return The integer number of the DynamicalSystemXML
   */
  inline int getNumber()
  {
    return SiconosDOMTreeTools::getIntegerAttributeValue(this->rootDynamicalSystemXMLNode, NUMBER_ATTRIBUTE);
  }


  /** \fn int getType()
   *   \brief Return the type of the DynamicalSystemXML
   *   \return The string type of the DynamicalSystemXML
   */
  inline std::string getType()
  {
    //return SiconosDOMTreeTools::getStringAttributeValue(this->rootDynamicalSystemXMLNode, TYPE_ATTRIBUTE);
    std::string res((char*)this->rootDynamicalSystemXMLNode->name);
    return res;
  }


  /** \fn string getId()
   *   \brief Return the id of the DynamicalSystemXML
   *   \return The string id of the DynamicalSystemXML
   */
  inline std::string getId()
  {
    return  SiconosDOMTreeTools::getStringContentValue(this->idNode);
  }

  /** \fn void setId(string s)
   *   \brief allows to save the id of the DynamicalSystemXML
   *   \param Integer : The string s of the DynamicalSystemXML
   */
  inline void setId(std::string s)
  {
    if (this->hasId() == false)
    {
      this->idNode = SiconosDOMTreeTools::createStringNode(this->rootDynamicalSystemXMLNode, ID_ATTRIBUTE, s);
    }
    else SiconosDOMTreeTools::setStringContentValue(this->idNode, s);
  }

  /** \fn int getN()
   *   \brief Return the n of the DynamicalSystemXML
   *   \return The integer n of the DynamicalSystemXML
   */
  inline int getN()
  {
    return  SiconosDOMTreeTools::getIntegerContentValue(this->nNode);
  }

  /** \fn void setN(int n)
   *   \brief allows to save the n of the DynamicalSystemXML
   *   \param Integer : The integer n of the DynamicalSystemXML
   */
  inline void setN(int n)
  {
    if (this->hasN() == false)
    {
      this->nNode = SiconosDOMTreeTools::createIntegerNode(this->rootDynamicalSystemXMLNode, DS_N, n);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(this->nNode, n);
  }


  /** \fn SimpleVector getX0()
   *   \brief Returns the X0 vector of the DynamicalSystemXML
   *   \return SimpleVector : X0 vector of the DynamicalSystemXML
   */
  inline /*SiconosVector*/SimpleVector getX0()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->x0Node);
  }

  /** \fn void setX0(SiconosVector *v)
   *   \brief allows to set the X0 of the DynamicalSystemXML
   *   \param The X0 SiconosVector to save
   */
  inline void setX0(SiconosVector *v)
  {
    if (!hasX0())
      x0Node = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, DS_X0, *v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(x0Node, *v);
  }

  /** \fn SimpleVector getX()
   *   \brief Returns the X vector of the DynamicalSystemXML
   *   \return SimpleVector : X vector of the DynamicalSystemXML
   */
  inline /*SiconosVector*/SimpleVector getX()
  {
    OUT("inline /*SiconosVector*/SimpleVector getX()");
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->xNode);
  }

  /** \fn void setX(SiconosVector *v)
   *   \brief allows to save the X of the DynamicalSystemXML
   *   \return The X SiconosVector to save
   */
  inline void setX(SiconosVector *v)
  {
    if (this->hasX() == false)
    {
      this->xNode = SiconosDOMTreeTools::createVectorNode(this->rootDynamicalSystemXMLNode, DS_X, *v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(this->xNode, *v);
  }

  /** \fn SimpleVector getXDot()
   *   \brief Returns the XDot vector of the DynamicalSystemXML
   *   \return SimpleVector : XDot vector of the DynamicalSystemXML
   */
  inline /*SiconosVector*/SimpleVector getXDot()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->xDotNode);
  }

  /** \fn void setXDot(SiconosVector *v)
   *   \brief allows to save the XDot of the DynamicalSystemXML
   *   \param The XDot SiconosVector to save
   */
  inline void setXDot(SiconosVector *v)
  {
    if (this->hasXDot() == false)
    {
      this->xDotNode = SiconosDOMTreeTools::createVectorNode(this->rootDynamicalSystemXMLNode, DS_XDOT, *v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(this->xDotNode, *v);
  }


  /** \fn SiconosMemoryXML* getXMemoryXML()
   *   \brief Returns the xMemoryXML* of the DynamicalSystemXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXML* getXMemoryXML()
  {
    return this->xMemoryXML;
  }

  /** \fn void setXMemory(SiconosMemory* smem)
   *   \brief allows to save the XMemory of the DynamicalSystemXML
   *   \param SiconosMemory* smem : SiconosMemory to save
   */
  inline void setXMemory(SiconosMemory* smem)
  {
    if (this->hasXMemory() == false)
    {
      this->xMemoryXML = new SiconosMemoryXML(NULL, this->rootDynamicalSystemXMLNode, DS_XMEMORY);
      this->xMemoryNode = this->xMemoryXML->getSiconosMemoryXMLNode();

      this->xMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->xMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
    else
    {
      this->xMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->xMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
  }


  /** \fn SiconosMemoryXML* getXDotMemoryXML()
   *   \brief Returns the xDotMemoryXML* of the DynamicalSystemXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXML* getXDotMemoryXML()
  {
    return xDotMemoryXML;
  }

  /** \fn void setXDotMemory(SiconosMemory* smem)
   *   \brief allows to save the xDotMemory of the DynamicalSystemXML
   *   \param SiconosMemory* smem : SiconosMemory to save
   */
  inline void setXDotMemory(SiconosMemory* smem)
  {

    if (this->hasXDotMemory() == false)
    {
      this->xDotMemoryXML = new SiconosMemoryXML(NULL, this->rootDynamicalSystemXMLNode, DS_XDOTMEMORY);
      this->xDotMemoryNode = this->xDotMemoryXML->getSiconosMemoryXMLNode();

      this->xDotMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->xDotMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
    else
    {
      this->xDotMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->xDotMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
  }

  /** \fn void setR(SimpleVector *r)
   *   \brief allows to save the R of the DynamicalSystemXML
   *   \param SiconosVector R of the DynamicalSystemXML
   */
  inline void setR(SimpleVector *r)
  {
    if (this->hasR() == false)
    {
      this->rNode = SiconosDOMTreeTools::createVectorNode(this->rootDynamicalSystemXMLNode, DS_R, *r);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(this->rNode, *r);
  }

  /** \fn SiconosMemoryXML* getRMemoryXML()
   *   \brief Returns the rMemoryXML* of the DynamicalSystemXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXML* getRMemoryXML()
  {
    return this->rMemoryXML;
  }

  /** \fn void setRMemory(SiconosMemory* smem)
   *   \brief allows to save the rMemory of the DynamicalSystemXML
   *   \param SiconosMemory* smem : SiconosMemory to save
   */
  inline void setRMemory(SiconosMemory* smem)
  {
    if (this->hasRMemory() == false)
    {
      this->rMemoryXML = new SiconosMemoryXML(NULL, this->rootDynamicalSystemXMLNode, DS_RMEMORY);
      this->rMemoryNode = this->rMemoryXML->getSiconosMemoryXMLNode();

      this->rMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->rMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
    else
    {
      this->rMemoryXML->setSiconosMemorySize(smem->getMemorySize());
      this->rMemoryXML->setSiconosMemoryVector(smem->getVectorMemory());
    }
  }

  /** \fn int getStepsInMemory()
   *   \brief Returns the steps in memory for the DynamicalSystemXML
   *   \return The integer number of steps in memory for the DynamicalSystemXML
   */
  inline int getStepsInMemory()
  {
    return  SiconosDOMTreeTools::getIntegerContentValue(this->stepsInMemoryNode);
  }

  /** \fn inline void setStepsInMemory(int nb)
   *   \brief allows to save the steps in memory for the DynamicalSystemXML
   *   \param The integer number of steps in memory to save
   */
  inline void setStepsInMemory(int nb)
  {
    if (this->hasStepsInMemory() == false)
    {
      this->stepsInMemoryNode = SiconosDOMTreeTools::createIntegerNode(this->rootDynamicalSystemXMLNode, DS_STEPSINMEMORY, nb);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(this->stepsInMemoryNode, nb);
  }

  /** \fn BoundaryConditionXML * getBoundaryConditionXML()
   *   \brief Returns the BoundaryConditionXML pointer of the DynamicalSystemXML
   *   \return the BoundaryConditionXML pointer of the DynamicalSystemXML ; NULL if DynamicalSystemXML does not have
   */
  inline BoundaryConditionXML * getBoundaryConditionXML()
  {
    return this->boundaryConditionXML;
  }


  /** \fn int getVectorFieldPlugin()
   *   \brief Returns the plugin for the DynamicalSystemXML
   *   \return string which defines the plugin for the DynamicalSystemXML
   */
  inline std::string getVectorFieldPlugin()
  {
    return  SiconosDOMTreeTools::getStringAttributeValue(this->vectorFieldNode, PLUGIN_ATTRIBUTE);
  }

  /** \fn inline void setVectorFieldPlugin(string plugin)
   *   \brief allows to save the the vectorFieldPlugin for the DynamicalSystemXML
   *   \param The string corresponding to the plugin to save
   */
  inline void setVectorFieldPlugin(std::string plugin)
  {
    if (this->hasVectorFieldPlugin() == false)
    {
      this->vectorFieldNode = SiconosDOMTreeTools::createSingleNode(this->rootDynamicalSystemXMLNode, DS_VECTORFIELD);
      xmlNewProp(this->vectorFieldNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->vectorFieldNode, PLUGIN_ATTRIBUTE, plugin);
  }

  /** \fn int getComputeJacobianXPlugin()
   *   \brief Returns the plugin for the DynamicalSystemXML
   *   \return string which defines the plugin for the DynamicalSystemXML
   */
  inline std::string getComputeJacobianXPlugin()
  {
    return  SiconosDOMTreeTools::getStringAttributeValue(this->computeJacobianXNode, PLUGIN_ATTRIBUTE);
  }

  /** \fn inline void setComputeJacobianXPlugin(string plugin)
   *   \brief allows to save the the vectorFieldPlugin for the DynamicalSystemXML
   *   \param The string corresponding to the plugin to save
   */
  inline void setComputeJacobianXPlugin(std::string plugin)
  {
    if (this->hasComputeJacobianXPlugin() == false)
    {
      this->computeJacobianXNode = SiconosDOMTreeTools::createSingleNode(this->rootDynamicalSystemXMLNode, DS_COMPUTEJACOBIANX);
      xmlNewProp(this->computeJacobianXNode, (xmlChar*)(PLUGIN_ATTRIBUTE.c_str()), (xmlChar*)plugin.c_str());
    }
    else SiconosDOMTreeTools::setStringAttributeValue(this->computeJacobianXNode, PLUGIN_ATTRIBUTE, plugin);
  }

  /** \fn bool hasN()
   *  \brief returns true if nNode is defined
   *  \return true if nNode is defined
   */
  inline bool hasN()
  {
    return (this->nNode != NULL);
  }

  /** \fn bool hasId()
   *  \brief returns true if idNode is defined
   *  \return true if idNode is defined
   */
  inline bool hasId()
  {
    return (this->idNode != NULL);
  }

  /** \fn bool hasX()
   *  \brief returns true if xNode is defined
   *  \return true if xNode is defined
   */
  inline bool hasX()
  {
    return (this->xNode != NULL);
  }

  /** \fn bool hasXDot()
   *  \brief returns true if xDotNode is defined
   *  \return true if xDotNode is defined
   */
  inline bool hasXDot()
  {
    return (this->xDotNode != NULL);
  }

  /** \fn bool hasXMemory()
   *  \brief returns true if xMemoryNode is defined
   *  \return true if xMemoryNode is defined
   */
  inline bool hasXMemory()
  {
    return (this->xMemoryNode != NULL);
  }

  /** \fn bool hasXDotMemory()
   *  \brief returns true if xDotMemoryNode is defined
   *  \return true if xDotMemoryNode is defined
   */
  inline bool hasXDotMemory()
  {
    return (this->xDotMemoryNode != NULL);
  }

  /** \fn bool hasX0()
   *  \brief returns true if x0Node is defined
   *  \return true if x0Node is defined
   */
  inline bool hasX0()
  {
    return (this->x0Node != NULL);
  }

  /** \fn bool hasStepsInMemory()
   *  \brief returns true if stepsInMemoryNode is defined
   *  \return true if stepsInMemoryNode is defined
   */
  inline bool hasStepsInMemory()
  {
    return (this->stepsInMemoryNode != NULL);
  }

  /** \fn bool hasR()
   *  \brief returns true if R is defined
   *  \return true if R is defined
   */
  inline bool hasR()
  {
    return (this->rNode != NULL);
  }

  /** \fn bool hasRMemory()
   *  \brief returns true if RMemory is defined
   *  \return true if RMemory is defined
   */
  inline bool hasRMemory()
  {
    return (this->rMemoryNode != NULL);
  }

  /** \fn bool hasBoundaryCondition()
   *  \brief returns true if boundaryConditionNode is defined
   *  \return true if boundaryConditionNode is defined
   */
  inline bool hasBoundaryCondition()
  {
    return (this->boundaryConditionNode != NULL);
  }

  /** \fn bool hasComputeJacobianXPlugin()
   *  \brief returns true if computeJacobianXNode is defined
   *  \return true if computeJacobianXNode is defined
   */
  inline bool hasComputeJacobianXPlugin()
  {
    return (this->computeJacobianXNode != NULL);
  }

  /** \fn bool hasVectorFieldPlugin()
   *  \brief returns true if vectorFieldNode is defined
   *  \return true if vectorFieldNode is defined
   */
  inline bool hasVectorFieldPlugin()
  {
    return (this->vectorFieldNode != NULL);
  }

  /** \fn void updateDynamicalSystemXML(xmlNode*, DynamicalSystem*, BoundaryCondition*)
   *   \brief makes the operations to add a DynamicalSystem to the NonSmoothDynamicalSystemXML
   *   \param xmlNode* : the root node of this DynamicalSystem
   *   \param DynamicalSystem* : the DynamicalSystem of this DynamicalSystemXML
   *   \param BoundaryCondition* : the BoundaryCondition of the DynamicalSystem if the NonSmoothDynamicalSystem is BVP (optional)
   */
  void updateDynamicalSystemXML(xmlNode*, DynamicalSystem*, BoundaryCondition* bc = NULL);

  /** \fn void loadDynamicalSystem( DynamicalSystem* )
   *   \brief loads the depending data of the DynamicalSystem into the DynamicalSystemXML (the BoundaryCondition if exists)
   *   \param DynamicalSystem* : the DynamicalSystem of this DynamicalSystemXML
   */
  void loadDynamicalSystem(DynamicalSystem*);


  /** \fn DSInputOutputXML* getDSInputOutputXML(int number)
   *   \brief Return the DSInputOutputXML with id number
   *   \param number : int number : the number of the DSInputOutputXML to return
   *  \exception XMLException
   *   \return the DSInputOutputXML of number number, NULL if doesn't exist
   */
  DSInputOutputXML* getDSInputOutputXML(int number);

  /** \fn inline vector<int> getDSInputOutputNumbers();
   *   \brief Allows to know the defined DSInputOutputs
   *   \return vector DSInputOutputs integer numbers
   */
  inline std::vector<int> getDSInputOutputNumbers()
  {
    return definedDSInputOutputNumbers;
  }

  /** \fn inline vector<int> getDSInputOutputNumbers()
   *   \brief Allows to know the defined DSInputOutputs
   *   \return vector DSInputOutputs integer numbers
   */
  void setDSInputOutputXML(std::map<int, DSInputOutputXML*> m);

  /** \fn void loadDSInputOutputXML(xmlNode * )
   *   \brief Builds DSInputOutputXML objects from a DOM tree describing DSInputOutputs
   *   \param xmlNode* : the DSInputOutputs DOM tree
   *   \exception XMLException : if a number relating to an DSInputOutput declares in the NSDS is already used
   */
  //void loadDSInputOutputXML(xmlNode * rootdsioNode);

  // ===== U AND T MANAGEMENT =====

  // --- uSize ---

  /** \fn int getUSize()
   *   \brief get size of vector u
   */
  inline const unsigned int getUSize() const
  {
    if (!hasUSize())
      XMLException::selfThrow("DynamicalSystemXML - getUSize: this node does not exist");
    return  SiconosDOMTreeTools::getIntegerContentValue(uSizeNode);
  }

  // --- u ---

  /** \fn inline string getUPlugin()
   *   \brief Return the u Plugin name of the DynamicalSystemXML
   *   \return The u Plugin name of the DynamicalSystemXML
   *  \exception XMLException
   */
  inline const std::string getUPlugin() const
  {
    if (!isUPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getUPlugin : u is not calculated from a plugin ; u vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(uNode, DS_VECTORPLUGIN);
  }

  /** \fn inline SimpleVector getUVector()
   *   \brief Return u vector of the DynamicalSystemXML
   *   \return SimpleVector : u of DynamicalSystemXML
   *  \exception XMLException
   */
  inline const SimpleVector getUVector() const
  {
    if (isUPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getUVector : u vector is not given ; u is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(uNode);
  }

  /** \fn inline void setUVector(SiconosVector *v)
   *   \brief allows to save the u vector of the DynamicalSystemXML
   *   \param SiconosVector *u : SiconosVector U to save
   */
  inline void setUVector(const SiconosVector& v)
  {
    if (uNode != NULL)
      SiconosDOMTreeTools::setSiconosVectorNodeValue(uNode, v);
    else uNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, DS_U, v);
  }

  // --- T ---

  /** \fn SiconosMatrix getTMatrix()
   *   \brief get T Matrix
   *   \return a SiconosMatrix
   */
  inline const SiconosMatrix getTMatrix() const
  {
    if (isTPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getT: T is not given but calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosMatrixValue(TNode);
  }

  /** \fn inline string getTPlugin()
   *  \brief get Plugin name to compute T
   *  \exception XMLException
   */
  inline const std::string getTPlugin() const
  {
    if (!isTPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getTPlugin : T is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(TNode, DS_MATRIXPLUGIN);
  }

  /** \fn void setT(SiconosMatrix *m)
   *   \brief save T
   *   \param The SiconosMatrix to save
   */
  inline void setT(const SiconosMatrix &m)
  {
    if (TNode != NULL)
      SiconosDOMTreeTools::setSiconosMatrixNodeValue(TNode, m);
    else TNode = SiconosDOMTreeTools::createMatrixNode(rootDynamicalSystemXMLNode, DS_T, m);
  }
  /** \fn bool isUPlugin()
   *   \brief Return true if u is calculated from a plugin
   */
  inline bool isUPlugin() const
  {
    return xmlHasProp((xmlNodePtr)uNode, (xmlChar *) DS_VECTORPLUGIN.c_str());
  }

  /** \fn bool isTPlugin()
   *   \brief Return true if T is calculated from a plugin
   */
  inline bool isTPlugin() const
  {
    return xmlHasProp((xmlNodePtr)TNode, (xmlChar *) DS_MATRIXPLUGIN.c_str());
  }

  /** \fn bool hasXX()
   * \brief return true if XXnode exists */
  inline bool hasUSize() const
  {
    return (uSizeNode != NULL);
  }
  inline bool hasU() const
  {
    return (uNode != NULL);
  }
  inline bool hasT() const
  {
    return (TNode != NULL);
  }

protected:
  xmlNode * rootDynamicalSystemXMLNode;
  xmlNode * parentNode;

  //Object
  BoundaryConditionXML * boundaryConditionXML;  //Maybe not defined (if not BVP NonSmoothDynamicalSystem)


  /** \fn void loadBoundaryConditionXML(xmlNode * rootBoundaryConditionNode)
   *   \brief Build BoundaryConditionXML object from a DOM tree describing BoundaryCondition
   *   \param rootBoundaryConditionXMLNode : the BoundaryCondition DOM tree
   *   \exception XMLException : if the type of the BoundaryCondition given in the DOM tree does not exist
   */
  void loadBoundaryConditionXML(xmlNode * rootBoundaryConditionNode);

  SiconosMemoryXML * xMemoryXML;
  SiconosMemoryXML * xDotMemoryXML;
  SiconosMemoryXML * rMemoryXML;

  /* Map of DSInputOutputs */
  std::map<int, DSInputOutputXML*> dsInputOutputXMLMap;

  /* vector of DSInputOutput numbers*/
  std::vector<int> definedDSInputOutputNumbers;

private:
  //Nodes

  xmlNode * idNode;
  xmlNode * nNode;
  xmlNode * x0Node;
  xmlNode * xNode;
  xmlNode * xDotNode;
  xmlNode * xMemoryNode;
  xmlNode * xDotMemoryNode;
  xmlNode * stepsInMemoryNode;
  xmlNode * vectorFieldNode;
  xmlNode * computeJacobianXNode;
  xmlNode * boundaryConditionNode;
  xmlNode * dsInputOutputNode;
  xmlNode * rNode;
  xmlNode * rMemoryNode;
  xmlNode * uSizeNode;
  xmlNode * uNode;
  xmlNode * TNode;

};

#endif
