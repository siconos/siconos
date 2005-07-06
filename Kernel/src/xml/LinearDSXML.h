
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
const std::string LDS_B = "b";
const std::string LDS_E = "E";
const std::string LDS_U = "u";

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
    if (isAPlugin())
      XMLException::selfThrow("LinearDSXML - getA: A is not given but calculated from a plugin");
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

  /** \fn SiconosMatrix getE()
   *   \brief Return the E of the LinearDSXML
   *   \return The E SiconosMatrix of the LinearDSXML
   */
  inline const SiconosMatrix getE() const
  {
    if (isEPlugin())
      XMLException::selfThrow("LinearDSXML - getE: E is not given but calculated from a plugin");
    return  SiconosDOMTreeTools::getSiconosMatrixValue(ENode);
  }

  /** \fn inline string getEPlugin()
   *   \brief Return the E Plugin name of the LinearDSXML
   *  \exception XMLException
   */
  inline const std::string getEPlugin() const
  {
    if (!isEPlugin())
      XMLException::selfThrow("LinearDSXML - getAPlugin : E is not loaded from a plugin");
    return  SiconosDOMTreeTools::getStringAttributeValue(ENode, LDS_MATRIXPLUGIN);
  }

  /** \fn void setE(SiconosMatrix *m)
   *   \brief allows to save the E of the LinearDSXML
   *   \return The E SiconosMatrix to save
   */
  inline void setE(const SiconosMatrix &m)
  {
    if (ENode != NULL)
      SiconosDOMTreeTools::setSiconosMatrixNodeValue(ENode, m);
    else ENode = SiconosDOMTreeTools::createMatrixNode(rootDSXMLNode, LDS_E, m);
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

  /** \fn inline string getBPlugin()
   *   \brief Return the b Plugin name of the LinearDSXML
   *   \return The b Plugin name of the LinearDSXML
   *  \exception XMLException
   */
  inline const std::string getBPlugin() const
  {
    if (!isBPlugin())
      XMLException::selfThrow("LinearDSXML - getUPlugin : b is not calculated from a plugin ; b vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(bNode, LDS_VECTORPLUGIN);
  }

  /** \fn inline SimpleVector getBVector()
   *   \brief Return b vector of the LinearDSXML
   *   \return SimpleVector : value of b of LinearDSXML
   *  \exception XMLException
   */
  inline const SimpleVector getBVector() const
  {
    if (isBPlugin())
      XMLException::selfThrow("LinearDSXML - getBVector : b vector is not given ; b is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(bNode);
  }

  /** \fn inline void setBVector(SiconosVector *v)
   *   \brief allows to save the b vector of the LinearDSXML
   *   \return The b SimpleVector to save
   */
  inline void setBVector(const SiconosVector& v)
  {
    if (bNode != NULL)
      SiconosDOMTreeTools::setSiconosVectorNodeValue(bNode, v);
    else bNode = SiconosDOMTreeTools::createVectorNode(rootDSXMLNode, LDS_B, v);
  }

  /** \fn bool isAPlugin()
   *   \brief Return true if A is calculated from a plugin
   */
  inline bool isAPlugin() const
  {
    return xmlHasProp((xmlNodePtr)ANode, (xmlChar *) LDS_MATRIXPLUGIN.c_str());
  }

  /** \fn bool isBPlugin()
   *   \brief Return true if b is calculated from a plugin
   */
  inline bool isBPlugin() const
  {
    return xmlHasProp((xmlNodePtr)bNode, (xmlChar *) LDS_VECTORPLUGIN.c_str());
  }
  /** \fn bool isUPlugin()
   *   \brief Return true if u is calculated from a plugin
   */
  inline bool isUPlugin() const
  {
    return xmlHasProp((xmlNodePtr)uNode, (xmlChar *) LDS_VECTORPLUGIN.c_str());
  }

  /** \fn bool isEPlugin()
   *   \brief Return true if E is calculated from a plugin
   */
  inline bool isEPlugin() const
  {
    return xmlHasProp((xmlNodePtr)ENode, (xmlChar *) LDS_MATRIXPLUGIN.c_str());
  }

  /** \fn bool hasXX()
   * \brief return true if XXnode exists */
  inline bool hasA() const
  {
    return (ANode != NULL);
  }
  inline bool hasB() const
  {
    return (bNode != NULL);
  }
  inline bool hasUSize() const
  {
    return (uSizeNode != NULL);
  }
  inline bool hasU() const
  {
    return (uNode != NULL);
  }
  inline bool hasE() const
  {
    return (ENode != NULL);
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
  xmlNode * bNode;
  xmlNode * uSizeNode;
  xmlNode * uNode;
  xmlNode * ENode;
};

#endif
