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

/** \class DynamicalSystemXML
 *   \brief This class manages DynamicalSystem data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.1.4.
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
const std::string DS_JACOBIANX = "jacobianX";
const std::string DS_MATRIXPLUGIN = "matrixPlugin";
const std::string DS_VECTORPLUGIN = "vectorPlugin";

class DynamicalSystemXML
{

protected:

  xmlNodePtr rootDynamicalSystemXMLNode;/**< root node  */
  xmlNodePtr parentNode; /**< node that owns root node (NonSmoothDynamicalSystemNode)  */
  xmlNodePtr nNode;  /**< nimber of degrees of freedom  */
  xmlNodePtr x0Node;/**< initial state */
  xmlNodePtr xNode;/**<  state (usefull is start from recovery xml-file*/
  xmlNodePtr stepsInMemoryNode; /**< size of memory */
  xmlNodePtr xMemoryNode;/**<  memory vector for x*/
  xmlNodePtr xDotMemoryNode;/**< memory vector for xDot (usefull?) */
  xmlNodePtr rMemoryNode;/**< memory vector for r */
  xmlNodePtr vectorFieldNode;/**< f(x,t) */
  xmlNodePtr jacobianXNode;/**< jacobian of f according to x */
  xmlNodePtr boundaryConditionNode;/**< boundary conditions */
  xmlNodePtr dsInputOutputNode;/**< ds input-output */
  xmlNodePtr uSizeNode;/**< size of control vector */
  xmlNodePtr uNode;/**< control term */
  xmlNodePtr TNode;/**< coefficient of control term */

  BoundaryConditionXML * boundaryConditionXML;/**< Boundary conditions obxml object */
  SiconosMemoryXML * xMemoryXML;/** <xml object for xMemory*/
  SiconosMemoryXML * xDotMemoryXML;/**<xml object for xDotMemory*/
  SiconosMemoryXML * rMemoryXML;/**<xml object for rMemory*/
  std::map<int, DSInputOutputXML*> dsInputOutputXMLMap;/**< list of dsinput-output, with an int as a key identifier*/
  std::vector<int> definedDSInputOutputNumbers;/**<  useless at the time */

  /** \fn void loadBoundaryConditionXML(xmlNodePtr rootBoundaryConditionNode)
   *   \brief Build BoundaryConditionXML object from a DOM tree describing BoundaryCondition
   *   \param rootBoundaryConditionXMLNode : the BoundaryCondition DOM tree
   *   \exception XMLException : if the type of the BoundaryCondition given in the DOM tree does not exist
   */
  void loadBoundaryConditionXML(xmlNodePtr rootBoundaryConditionNode);

public:

  /** \fn DynamicalSystemXML()
   *   \brief Default constructor */
  DynamicalSystemXML();

  /** \fn ~DynamicalSystemXML()
   *   \brief Destructor */
  virtual ~DynamicalSystemXML();

  /** \fn DynamicalSystemXML(xmlNodePtr DynamicalSystemNode, const bool&)
   *   \brief Build a DynamicalSystemXML object from the DynamicalSystem node of the xml DOMtree
   *   \param an xmlNodePtr DynamicalSystemNode
   *   \param bool isBVP : if true, NonSmoothDynamicalSystem is a Boundary value problem
   */
  DynamicalSystemXML(xmlNodePtr DynamicalSystemNode, const bool&);

  /** \fn const int getNumber() const
   *   \brief Return the number of the DynamicalSystem
   *   \return an integer
   */
  inline const int getNumber() const
  {
    return SiconosDOMTreeTools::getIntegerAttributeValue(rootDynamicalSystemXMLNode, NUMBER_ATTRIBUTE);
  }

  /** \fn const string getType() const
   *   \brief Return the type of the DynamicalSystem
   *   \return a string
   */
  inline const std::string getType() const
  {
    std::string res((char*)rootDynamicalSystemXMLNode->name);
    return res;
  }

  /** \fn const std::string getId() const
   *   \brief Return the id of the DynamicalSystem
   *   \return The string id of the DynamicalSystem
   */
  inline const std::string getId() const
  {
    return SiconosDOMTreeTools::getStringAttributeValue(rootDynamicalSystemXMLNode, "Id");
  }

  /** \fn bool hasId() const
   *  \brief return true if id is given
   *  \return a bool
   */
  inline bool hasId() const
  {
    return (xmlHasProp(rootDynamicalSystemXMLNode, (xmlChar*)"Id"));
  }

  /** \fn   void setId(const std::string& s);
   *   \brief to save the id of the DynamicalSystem
   *   \param a string
   */
  inline void setId(const std::string& s)
  {
    SiconosDOMTreeTools::setStringAttributeValue(rootDynamicalSystemXMLNode, "Id", s);
  }

  /** \fn const unsigned int getN() const
   *   \brief Return the number of degrees of freedom of the DynamicalSystem
   *   \return an unsigned integer
   */
  inline const unsigned int getN() const
  {
    return  SiconosDOMTreeTools::getIntegerContentValue(nNode);
  }

  /** \fn void setN(const unsigned int& n)
   *   \brief to save the number of degrees of freedom of the DynamicalSystem
   *   \param an unsigned int
   */
  inline void setN(const unsigned int& n)
  {
    if (!hasN())
      nNode = SiconosDOMTreeTools::createIntegerNode(rootDynamicalSystemXMLNode, DS_N, n);
    else SiconosDOMTreeTools::setIntegerContentValue(nNode, n);
  }

  /** \fn  const SimpleVector getX0() const
   *   \brief Returns initial state of the DynamicalSystem (x0)
   *   \return a SimpleVector
   */
  inline const SimpleVector getX0() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(x0Node);
  }

  /** \fn void setX0(const SiconosVector& v)
   *   \brief save x0 of the DynamicalSystem
   *   \param a SiconosVector
   */
  inline void setX0(const SiconosVector& v)
  {
    if (!hasX0())
      x0Node = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, DS_X0, v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(x0Node, v);
  }

  /** \fn const SimpleVector getX() const
   *   \brief Returns the x state-vector of the DynamicalSystem
   *   \return SimpleVector
   */
  inline const SimpleVector getX() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(xNode);
  }

  /** \fn void setX(const SiconosVector& v)
   *   \brief save x of the DynamicalSystem
   *   \param a SiconosVector
   */
  inline void setX(const SiconosVector& v)
  {
    if (!hasX())
      xNode = SiconosDOMTreeTools::createVectorNode(rootDynamicalSystemXMLNode, DS_X, v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(xNode, v);
  }

  /** \fn int getStepsInMemory()
   *   \brief Returns the steps in memory for the DynamicalSystemXML
   *   \return The integer number of steps in memory for the DynamicalSystemXML
   */
  inline const unsigned int getStepsInMemory() const
  {
    return  SiconosDOMTreeTools::getIntegerContentValue(stepsInMemoryNode);
  }

  /** \fn inline void setStepsInMemory(const unsigned int& nb)
   *   \brief to save the steps in memory for the DynamicalSystemXML
   *   \param an unsigned int
   */
  void setStepsInMemory(const unsigned int&);

  /** \fn SiconosMemoryXML* getXMemoryXML() const
   *   \brief Returns the xMemoryXML* of the DynamicalSystemXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXML* getXMemoryXML() const
  {
    return xMemoryXML;
  }

  /** \fn void setXMemory(const SiconosMemory& smem)
   *   \brief to save the XMemory of the DynamicalSystemXML
   *   \param SiconosMemory* smem : SiconosMemory to save
   */
  void setXMemory(const SiconosMemory&);

  /** \fn SiconosMemoryXML* getXDotMemoryXML() const
   *   \brief Returns the xDotMemoryXML* of the DynamicalSystemXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXML* getXDotMemoryXML() const
  {
    return xDotMemoryXML;
  }

  /** \fn void setXDotMemory(const SiconosMemory& )
   *   \brief to save the xDotMemory of the DynamicalSystemXML
   *   \param SiconosMemory* smem : SiconosMemory to save
   */
  void setXDotMemory(const SiconosMemory&);

  /** \fn SiconosMemoryXML* getRMemoryXML() const
   *   \brief Returns the rMemoryXML* of the DynamicalSystemXML
   *   \return SiconosMemoryXML*
   */
  inline SiconosMemoryXML* getRMemoryXML() const
  {
    return rMemoryXML;
  }

  /** \fn void setRMemory( const SiconosMemory& );
   *   \brief to save the rMemory of the DynamicalSystemXML
   *   \param SiconosMemory*
   */
  void setRMemory(const SiconosMemory&);

  /** \fn BoundaryConditionXML * getBoundaryConditionXML() const
   *   \brief Returns the BoundaryConditionXML pointer of the DynamicalSystemXML
   *   \return the BoundaryConditionXML pointer of the DynamicalSystemXML ; NULL if DynamicalSystemXML does not have
   */
  inline BoundaryConditionXML * getBoundaryConditionXML() const
  {
    return boundaryConditionXML;
  }

  // === VectorField ===
  /** \fn inline string getVectorFieldPlugin()
   *  \brief Return the name of the vectorField plug-in
   *  \return a string
   */
  inline const std::string getVectorFieldPlugin() const
  {
    if (!isVectorFieldPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getVectorFieldPlugin : vectorField is not calculated from a plugin since a vectorField vector is given.");
    return  SiconosDOMTreeTools::getStringAttributeValue(vectorFieldNode, "vectorPlugin");
  }

  /** \fn const SimpleVector getVectorFieldVector() const
   *   \brief return vectorField vector
   *   \return SimpleVector
   */
  inline const SimpleVector getVectorFieldVector() const
  {
    if (isVectorFieldPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getVectorFieldVector : vectorField vector is not given since vectorField is calculated using a plug-in");
    return  SiconosDOMTreeTools::getSiconosVectorValue(vectorFieldNode);
  }

  /** \fn void setVectorFieldVector(const SiconosVector& v)
   *   \brief to save the vectorField vector
   *   \param a SiconosVector
   */
  void setVectorFieldVector(const SiconosVector&v);

  /** \fn void setVectorFieldPlugin(const string& plugin)
   *   \brief to save the VectorField plugin
   *   \param a string (name of the plug-in)
   */
  void setVectorFieldPlugin(const std::string& plugin);

  // === JacobianX ===
  /** \fn inline string getJacobianXPlugin()
   *  \brief Return the name of the jacobianX plug-in
   *  \return a string
   */
  inline const std::string getJacobianXPlugin() const
  {
    if (!isJacobianXPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getJacobianXPlugin : jacobianX is not calculated from a plugin since a jacobianX matrix is given.");
    return  SiconosDOMTreeTools::getStringAttributeValue(jacobianXNode, "matrixPlugin");
  }

  /** \fn const SimpleMatrix getJacobianXMatrix() const
   *   \brief return jacobianX matrix
   *   \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianXMatrix() const
  {
    if (isJacobianXPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getJacobianXMatrix : jacobianX matrix is not given since jacobianX is calculated using a plug-in");
    return  SiconosDOMTreeTools::getSiconosMatrixValue(jacobianXNode);
  }

  /** \fn void setJacobianXMatrix(const SiconosMatrix& v)
   *   \brief to save the jacobianX matrix
   *   \param a SiconosMatrix
   */
  void setJacobianXMatrix(const SiconosMatrix&v);

  /** \fn void setJacobianXPlugin(const string& plugin)
   *   \brief to save the JacobianX plugin of the LagrangianDSXML
   *   \param a string (name of the plug-in)
   */
  void setJacobianXPlugin(const std::string& plugin);

  /** \fn int getUSize()
   *   \brief get size of vector u
   */
  inline const unsigned int getUSize() const
  {
    if (!hasUSize())
      XMLException::selfThrow("DynamicalSystemXML - getUSize: this node does not exist");
    return  SiconosDOMTreeTools::getIntegerContentValue(uSizeNode);
  }

  /** \fn void setUSize(const unsigned int& us)
   *   \brief to save the size of vector u
   *   \param an unsigned int
   */
  inline void setUSize(const unsigned int& us)
  {
    if (!hasUSize())
      uSizeNode = SiconosDOMTreeTools::createIntegerNode(rootDynamicalSystemXMLNode, "uSize", us);
    else SiconosDOMTreeTools::setIntegerContentValue(uSizeNode, us);
  }

  // --- u ---

  /** \fn inline const string getUPlugin() const
   *   \brief Return the u plug-in name
   *   \return a string
   */
  inline const std::string getUPlugin() const
  {
    if (!isUPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getUPlugin : u is not calculated from a plugin ; u vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(uNode, DS_VECTORPLUGIN);
  }

  /** \fn inline const SimpleVector getUVector() const
   *   \brief Return u vector
   *   \return a SimpleVector
   *  \exception XMLException
   */
  inline const SimpleVector getUVector() const
  {
    if (isUPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getUVector : u vector is not given ; u is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(uNode);
  }

  /** \fn void setUVector(const SiconosVector & v)
   *   \brief to save the u vector
   *   \param a SiconosVector
   */
  void setUVector(const SiconosVector&);

  /** \fn void setUPlugin(const string& plugin)
   *   \brief to save the u plugin
   *   \param a string (name of the plug-in)
   */
  void setUPlugin(const std::string& plugin);

  // --- T ---

  /** \fn const SimpleMatrix getTMatrix() const
   *   \brief get T Matrix
   *   \return a SimpleMatrix
   */
  inline const SimpleMatrix getTMatrix() const
  {
    if (isTPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getT: T is not given but calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosMatrixValue(TNode);
  }

  /** \fn inline const string getTPlugin() const
   *  \brief get Plugin name to compute T
   *  \exception XMLException
   */
  inline const std::string getTPlugin() const
  {
    if (!isTPlugin())
      XMLException::selfThrow("DynamicalSystemXML - getTPlugin : T is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(TNode, DS_MATRIXPLUGIN);
  }

  /** \fn void setT(const SiconosMatrix& m)
   *   \brief save T
   *   \param The SiconosMatrix to save
   */
  void setT(const SiconosMatrix &m);

  /** \fn void setTPlugin(const string& plugin)
   *   \brief to save the T plugin
   *   \param a string (name of the plug-in)
   */
  void setTPlugin(const std::string& plugin);

  /** \fn bool hasN() const
   *  \brief returns true if nNode is defined
   *  \return true if nNode is defined
   */
  inline bool hasN() const
  {
    return (nNode != NULL);
  }

  /** \fn bool hasX0() const
   *  \brief returns true if x0Node is defined
   *  \return true if x0Node is defined
   */
  inline bool hasX0() const
  {
    return (x0Node != NULL);
  }

  /** \fn bool hasX() const
   *  \brief returns true if xNode is defined
   *  \return true if xNode is defined
   */
  inline bool hasX() const
  {
    return (xNode != NULL);
  }

  /** \fn bool hasStepsInMemory() const
   *  \brief returns true if stepsInMemoryNode is defined
   *  \return true if stepsInMemoryNode is defined
   */
  inline bool hasStepsInMemory() const
  {
    return (stepsInMemoryNode != NULL);
  }

  /** \fn bool hasXMemory() const
   *  \brief returns true if xMemoryNode is defined
   *  \return true if xMemoryNode is defined
   */
  inline bool hasXMemory()const
  {
    return (xMemoryNode != NULL);
  }

  /** \fn bool hasXDotMemory()  const
   *  \brief returns true if xDotMemoryNode is defined
   *  \return true if xDotMemoryNode is defined
   */
  inline bool hasXDotMemory() const
  {
    return (xDotMemoryNode != NULL);
  }

  /** \fn bool hasRMemory() const
   *  \brief returns true if RMemory is defined
   *  \return true if RMemory is defined
   */
  inline bool hasRMemory() const
  {
    return (rMemoryNode != NULL);
  }

  /** \fn bool hasVectorField() const
   *  \brief returns true if vectorFieldNode is defined
   *  \return true if vectorFieldNode is defined
   */
  inline bool hasVectorField() const
  {
    return (vectorFieldNode != NULL);
  }

  /** \fn bool hasJacobianX() const
   *  \brief returns true if jacobianXNode is defined
   *  \return true if jacobianXNode is defined
   */
  inline bool hasJacobianX() const
  {
    return (jacobianXNode != NULL);
  }

  /** \fn bool hasBoundaryCondition() const
   *  \brief returns true if boundaryConditionNode is defined
   *  \return true if boundaryConditionNode is defined
   */
  inline bool hasBoundaryCondition() const
  {
    return (boundaryConditionNode != NULL);
  }

  /** \fn bool hasUSize() const
  *  \brief returns true if uSizeNode is defined
  *  \return true if jacobianXNode is defined
  */
  inline bool hasUSize() const
  {
    return (uSizeNode != NULL);
  }

  /** \fn bool hasU() const
  *  \brief returns true if uNode is defined
  *  \return true if jacobianXNode is defined
  */
  inline bool hasU() const
  {
    return (uNode != NULL);
  }

  /** \fn bool hasT() const
  *  \brief returns true if TNode is defined
  *  \return true if jacobianXNode is defined
  */
  inline bool hasT() const
  {
    return (TNode != NULL);
  }

  /** \fn bool isVectorFieldPlugin()
   *   \brief Return true if vectorField is calculated from a plugin
   */
  inline bool isVectorFieldPlugin() const
  {
    return xmlHasProp((xmlNodePtr)vectorFieldNode, (xmlChar *) DS_VECTORPLUGIN.c_str());
  }

  /** \fn bool isJacobianXPlugin()
   *   \brief Return true if jacobianX is calculated from a plugin
   */
  inline bool isJacobianXPlugin() const
  {
    return xmlHasProp((xmlNodePtr)jacobianXNode, (xmlChar *) DS_MATRIXPLUGIN.c_str());
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

  /** \fn void updateDynamicalSystemXML(xmlNodePtr, DynamicalSystem*, BoundaryCondition*)
   *   \brief prepare object(s) to add a DynamicalSystem to the NonSmoothDynamicalSystemXML
   *   \param xmlNodePtr : the root node of this DynamicalSystem
   *   \param DynamicalSystem* : the DynamicalSystem of this DynamicalSystemXML
   *   \param BoundaryCondition* : the BoundaryCondition of the DynamicalSystem if the NonSmoothDynamicalSystem is BVP (optional)
   */
  void updateDynamicalSystemXML(xmlNodePtr, DynamicalSystem*, BoundaryCondition* bc = NULL);

  /** \fn void loadDynamicalSystem( DynamicalSystem* )
   *   \brief loads the depending data of the DynamicalSystem into the DynamicalSystemXML (the BoundaryCondition if exists)
   *   \param DynamicalSystem* : the DynamicalSystem of this DynamicalSystemXML
   */
  void loadDynamicalSystem(DynamicalSystem*);

  /** \fn DSInputOutputXML* getDSInputOutputXML(const int& number);
   *   \brief Return the DSInputOutputXML with id number
   *   \param number : int number : the number of the DSInputOutputXML to return
   *  \exception XMLException
   *   \return the DSInputOutputXML of number number, NULL if doesn't exist
   */
  DSInputOutputXML* getDSInputOutputXML(const int&);

  /** \fn inline const vector<int> getDSInputOutputNumbers(); const
   *   \brief Allows to know the defined DSInputOutputs
   *   \return vector DSInputOutputs integer numbers
   */
  inline const std::vector<int> getDSInputOutputNumbers() const
  {
    return definedDSInputOutputNumbers;
  }

  /** \fn inline vector<int> getDSInputOutputNumbers()
   *   \brief Allows to know the defined DSInputOutputs
   *   \return vector DSInputOutputs integer numbers
   */
  void setDSInputOutputXML(std::map<int, DSInputOutputXML*> m);

  /** \fn void loadDSInputOutputXML(xmlNodePtr )
   *   \brief Builds DSInputOutputXML objects from a DOM tree describing DSInputOutputs
   *   \param xmlNodePtr : the DSInputOutputs DOM tree
   *   \exception XMLException : if a number relating to an DSInputOutput declares in the NSDS is already used
   */
  //void loadDSInputOutputXML(xmlNodePtr rootdsioNode);
};

#endif
