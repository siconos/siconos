
/** \class LinearDSXML
 *   \brief This class manages LinearSystem DS data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 05/11/2004
 *
 *
 *
 * LinearDSXML allows to manage data of a LinearDS DOM tree.
 */

#ifndef __LINEARSYSTEMDSXML__
#define __LINEARSYSTEMDSXML__

#include "DSXML.h"

const std::string LDS_A = "A";
const std::string LDS_B = "B";

const std::string LDS_U = "u";
const std::string LDS_F = "f";

const std::string LDS_MATRIXPLUGIN = "matrixPlugin";
const std::string LDS_VECTORPLUGIN = "vectorPlugin";


class LinearDSXML : public DSXML
{
public:
  LinearDSXML();

  /** \fn LinearDSXML(xmlNode * LagrangianDSNode, bool isBVP)
   *   \brief Build a LinearDSXML object from a DOM tree describing a DS
   *   \param xmlNode * linearSystemDSNode : the linearSystemDS DOM tree
   *   \param bool isBVP : if NSDS is BVP, linearSystemDS has boundary condition
   */
  LinearDSXML(xmlNode * linearSystemDSNode, const bool& isBVP);

  ~LinearDSXML();

  /** \fn SiconosMatrix getA()
   *   \brief Return the A of the LinearDSXML
   *   \return The A SiconosMatrix of the LinearDSXML
   */
  inline const SiconosMatrix getA() const
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(ANode);
  }

  /** \fn inline string getAPlugin()
   *   \brief Return the A Plugin name of the LinearDSXML
   *  \exception XMLException
   */
  inline const std::string getAPlugin() const
  {
    if (!isAPlugin())
      XMLException::selfThrow("LinearDSXML - getAPlugin : A is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(ANode, LDS_MATRIXPLUGIN);
  }

  /** \fn void setA(SiconosMatrix *m)
   *   \brief allows to save the A of the LinearDSXML
   *   \return The A SiconosMatrix to save
   */
  inline void setA(const SiconosMatrix& m)
  {
    if (ANode != NULL)
      SiconosDOMTreeTools::setSiconosMatrixNodeValue(ANode, m);
    else
      ANode = SiconosDOMTreeTools::createMatrixNode(rootDSXMLNode, LDS_A, m);
  }

  /** \fn SiconosMatrix getB()
   *   \brief Return the B of the LinearDSXML
   *   \return The B SiconosMatrix of the LinearDSXML
   */
  inline const SiconosMatrix getB() const
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(BNode);
  }

  /** \fn inline string getBPlugin()
   *   \brief Return the B Plugin name of the LinearDSXML
   *  \exception XMLException
   */
  inline const std::string getBPlugin() const
  {
    if (!isBPlugin())
      XMLException::selfThrow("LinearDSXML - getAPlugin : A is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(BNode, LDS_MATRIXPLUGIN);
  }

  /** \fn void setB(SiconosMatrix *m)
   *   \brief allows to save the B of the LinearDSXML
   *   \return The B SiconosMatrix to save
   */
  inline void setB(const SiconosMatrix &m)
  {
    if (BNode != NULL)
      SiconosDOMTreeTools::setSiconosMatrixNodeValue(BNode, m);
    else BNode = SiconosDOMTreeTools::createMatrixNode(rootDSXMLNode, LDS_B, m);
  }

  /** \fn int getUSize()
   *   \brief get size of vector u
   */
  inline const unsigned int getUSize() const
  {
    if (!hasUSize())
      XMLException::selfThrow("LinearDSXML - getUSize: this node does not exist");
    return  SiconosDOMTreeTools::getIntegerContentValue(uSizeNode);
  }

  /** \fn inline string getUPlugin()
   *   \brief Return the u Plugin name of the LinearDSXML
   *   \return The u Plugin name of the LinearDSXML
   *  \exception XMLException
   */
  inline const std::string getUPlugin() const
  {
    if (!isUPlugin())
      XMLException::selfThrow("LinearDSXML - getUPlugin : u is not calculated from a plugin ; u vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(uNode, LDS_VECTORPLUGIN);
  }

  /** \fn inline SimpleVector getUVector()
   *   \brief Return u vector of the LinearDSXML
   *   \return SimpleVector : u of LinearDSXML
   *  \exception XMLException
   */
  inline const SimpleVector getUVector() const
  {
    if (isUPlugin())
      XMLException::selfThrow("LinearDSXML - getUVector : u vector is not given ; u is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(uNode);
  }

  /** \fn inline void setUVector(SiconosVector *v)
   *   \brief allows to save the u vector of the LinearDSXML
   *   \param SiconosVector *u : SiconosVector U to save
   */
  inline void setUVector(const SiconosVector& v)
  {
    if (uNode != NULL)
      SiconosDOMTreeTools::setSiconosVectorNodeValue(uNode, v);
    else uNode = SiconosDOMTreeTools::createVectorNode(rootDSXMLNode, LDS_U, v);
  }

  /** \fn inline string getFPlugin()
   *   \brief Return the f Plugin name of the LinearDSXML
   *   \return The f Plugin name of the LinearDSXML
   *  \exception XMLException
   */
  inline const std::string getFPlugin() const
  {
    if (!isFPlugin())
      XMLException::selfThrow("LinearDSXML - getUPlugin : f is not calculated from a plugin ; f vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(fNode, LDS_VECTORPLUGIN);
  }

  /** \fn inline SimpleVector getFVector()
   *   \brief Return f vector of the LinearDSXML
   *   \return SimpleVector : value of f of LinearDSXML
   *  \exception XMLException
   */
  inline const SimpleVector getFVector() const
  {
    if (isFPlugin())
      XMLException::selfThrow("LinearDSXML - getFVector : f vector is not given ; f is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(fNode);
  }

  /** \fn inline void setFVector(SiconosVector *v)
   *   \brief allows to save the f vector of the LinearDSXML
   *   \return The f SimpleVector to save
   */
  inline void setFVector(const SiconosVector& v)
  {
    if (fNode != NULL)
      SiconosDOMTreeTools::setSiconosVectorNodeValue(fNode, v);
    else fNode = SiconosDOMTreeTools::createVectorNode(rootDSXMLNode, LDS_F, v);
  }

  /** \fn bool isAPlugin()
   *   \brief Return true if A is calculated from a plugin
   */
  inline bool isAPlugin() const
  {
    return xmlHasProp((xmlNodePtr)ANode, (xmlChar *) LDS_MATRIXPLUGIN.c_str());
  }

  /** \fn bool isFPlugin()
   *   \brief Return true if f is calculated from a plugin
   */
  inline bool isFPlugin() const
  {
    return xmlHasProp((xmlNodePtr)fNode, (xmlChar *) LDS_VECTORPLUGIN.c_str());
  }
  /** \fn bool isUPlugin()
   *   \brief Return true if u is calculated from a plugin
   */
  inline bool isUPlugin() const
  {
    return xmlHasProp((xmlNodePtr)uNode, (xmlChar *) LDS_VECTORPLUGIN.c_str());
  }

  /** \fn bool isBPlugin()
   *   \brief Return true if B is calculated from a plugin
   */
  inline bool isBPlugin() const
  {
    return xmlHasProp((xmlNodePtr)BNode, (xmlChar *) LDS_MATRIXPLUGIN.c_str());
  }

  /** \fn bool hasXX()
   * \brief return true if XXnode exists */
  inline bool hasA() const
  {
    return (ANode != NULL);
  }
  inline bool hasF() const
  {
    return (fNode != NULL);
  }
  inline bool hasUSize() const
  {
    return (uSizeNode != NULL);
  }
  inline bool hasU() const
  {
    return (uNode != NULL);
  }
  inline bool hasB() const
  {
    return (BNode != NULL);
  }


  /** \fn void updateDynamicalSystemXML(xmlNode*, DynamicalSystem*, BoundaryCondition*)
   *   \brief makes the operations to add a DynamicalSystem to the NSDSXML
   *   \param xmlNode* : the root node of this DynamicalSystem
   *   \param DynamicalSystem* : the DynamicalSystem of this DSXML
   *   \param BoundaryCondition* : the BoundaryCondition of the DS if the NSDS is BVP (optional)
   */
  void updateDynamicalSystemXML(xmlNode*, DynamicalSystem*, BoundaryCondition* bc = NULL);


private:

  //Nodes
  xmlNode * ANode;
  xmlNode * fNode;
  xmlNode * uSizeNode;
  xmlNode * uNode;
  xmlNode * BNode;
};

#endif
