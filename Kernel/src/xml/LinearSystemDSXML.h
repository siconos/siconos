
/** \class LinearSystemDSXML
 *   \brief This class manages LinearSystem DS data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 05/11/2004
 *
 *
 *
 * LinearSystemDSXML allows to manage data of a LinearSystemDS DOM tree.
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


class LinearSystemDSXML : public DSXML
{
public:
  LinearSystemDSXML();

  /** \fn LinearSystemDSXML(xmlNode * LagrangianDSNode, bool isBVP)
   *   \brief Build a LinearSystemDSXML object from a DOM tree describing a DS
   *   \param xmlNode * linearSystemDSNode : the linearSystemDS DOM tree
   *   \param bool isBVP : if NSDS is BVP, linearSystemDS has boundary condition
   */
  LinearSystemDSXML(xmlNode * linearSystemDSNode, bool isBVP);

  ~LinearSystemDSXML();

  /** \fn SiconosMatrix getA()
   *   \brief Return the A of the LinearSystemDSXML
   *   \return The A SiconosMatrix of the LinearSystemDSXML
   */
  inline SiconosMatrix getA()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(ANode);
  }

  /** \fn SiconosMatrix getB()
   *   \brief Return the B of the LinearSystemDSXML
   *   \return The B SiconosMatrix of the LinearSystemDSXML
   */
  inline SiconosMatrix getB()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(BNode);
  }

  /** \fn void setA(SiconosMatrix *m)
   *   \brief allows to save the A of the LinearSystemDSXML
   *   \return The A SiconosMatrix to save
   */
  inline void setA(SiconosMatrix *m)
  {
    if (ANode != NULL)
    {
      //SiconosDOMTreeTools::setSiconosMatrixValue(ANode, *m);
      SiconosDOMTreeTools::setSiconosMatrixNodeValue(ANode, *m);
    }
    else
    {
      ANode = SiconosDOMTreeTools::createMatrixNode(rootDSXMLNode, LDS_A, *m);
    }
  }

  /** \fn void setB(SiconosMatrix *m)
   *   \brief allows to save the B of the LinearSystemDSXML
   *   \return The B SiconosMatrix to save
   */
  inline void setB(SiconosMatrix *m)
  {
    if (BNode != NULL)
    {
      //SiconosDOMTreeTools::setSiconosMatrixValue(BNode, *m);
      SiconosDOMTreeTools::setSiconosMatrixNodeValue(BNode, *m);
    }
    else BNode = SiconosDOMTreeTools::createMatrixNode(rootDSXMLNode, LDS_B, *m);
  }
  /////////////////////////////

  /** \fn inline string getUPlugin()
   *   \brief Return the u Plugin name of the LinearSystemDSXML
   *   \return The u Plugin name of the LinearSystemDSXML
   *  \exception XMLException
   */
  inline std::string getUPlugin()
  {
    if (!isUPlugin())
      XMLException::selfThrow("LinearSystemDSXML - getUPlugin : u is not calculated from a plugin ; u vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(uNode, LDS_VECTORPLUGIN);
  }

  /** \fn inline SimpleVector getUVector()
   *   \brief Return u vector of the LinearSystemDSXML
   *   \return SimpleVector : u of LinearSystemDSXML
   *  \exception XMLException
   */
  inline /*SiconosVector*/SimpleVector getUVector()
  {
    if (isUPlugin())
      XMLException::selfThrow("LinearSystemDSXML - getUVector : u vector is not given ; u is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(uNode);
  }

  /** \fn inline void setUVector(SiconosVector *v)
   *   \brief allows to save the u vector of the LinearSystemDSXML
   *   \param SiconosVector *u : SiconosVector U to save
   */
  inline void setUVector(SiconosVector *v)
  {
    if (uNode != NULL)
    {
      //SiconosDOMTreeTools::setSiconosVectorValue(uNode, *v);
      SiconosDOMTreeTools::setSiconosVectorNodeValue(uNode, *v);
    }
    else uNode = SiconosDOMTreeTools::createVectorNode(rootDSXMLNode, LDS_U, *v);
  }

  /** \fn inline string getFPlugin()
   *   \brief Return the f Plugin name of the LinearSystemDSXML
   *   \return The f Plugin name of the LinearSystemDSXML
   *  \exception XMLException
   */
  inline std::string getFPlugin()
  {
    if (!isFPlugin())
      XMLException::selfThrow("LinearSystemDSXML - getUPlugin : f is not calculated from a plugin ; f vector is given");
    return  SiconosDOMTreeTools::getStringAttributeValue(fNode, LDS_VECTORPLUGIN);
  }

  /** \fn inline SimpleVector getFVector()
   *   \brief Return f vector of the LinearSystemDSXML
   *   \return SimpleVector : value of f of LinearSystemDSXML
   *  \exception XMLException
   */
  inline /*SiconosVector*/ SimpleVector getFVector()
  {
    if (isFPlugin())
      XMLException::selfThrow("LinearSystemDSXML - getFVector : f vector is not given ; f is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(fNode);
  }

  /** \fn inline void setFVector(SiconosVector *v)
   *   \brief allows to save the f vector of the LinearSystemDSXML
   *   \return The f SimpleVector to save
   */
  inline void setFVector(SiconosVector *v)
  {
    if (fNode != NULL)
    {
      //SiconosDOMTreeTools::setSiconosVectorValue(fNode, *v);
      SiconosDOMTreeTools::setSiconosVectorNodeValue(fNode, *v);
    }
    else fNode = SiconosDOMTreeTools::createVectorNode(rootDSXMLNode, LDS_F, *v);
  }


  /** \fn bool isUPlugin()
   *   \brief Return true if u is calculated from a plugin
   *   \return True if u is calculated from plugin
   */
  inline bool isUPlugin()
  {
    return xmlHasProp((xmlNodePtr)uNode, (xmlChar *) LDS_VECTORPLUGIN.c_str());
  }

  /** \fn bool isFPlugin()
   *   \brief Return true if f is calculated from a plugin
   *   \return True if f is calculated from plugin
   */
  inline bool isFPlugin()
  {
    return xmlHasProp((xmlNodePtr)fNode, (xmlChar *) LDS_VECTORPLUGIN.c_str());
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
  xmlNode * BNode;
  xmlNode * uNode;
  xmlNode * fNode;


  //Methods

  /** \fn loadLinearSystemDSProperties()
   *   \brief load the different properties of a Linear System DS
   *   \exception XMLException : if a property of the LinearSystemDS lacks in the DOM tree
   */
  void loadLinearSystemDSProperties();

};

#endif
