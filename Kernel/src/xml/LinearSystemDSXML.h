
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


#include <vector>
#include <string>
#include <libxml/tree.h>

#include "DSXML.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"

#include "SiconosDOMTreeTools.h"


using namespace std;


const string LSDS_A = "A";
const string LSDS_B = "B";

const string LSDS_U = "u";
const string LSDS_F = "f";

const string LSDS_MATRIXPLUGIN = "matrixPlugin";
const string LSDS_VECTORPLUGIN = "vectorPlugin";


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
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->ANode);
  }

  /** \fn SiconosMatrix getB()
  *   \brief Return the B of the LinearSystemDSXML
  *   \return The B SiconosMatrix of the LinearSystemDSXML
  */
  inline SiconosMatrix getB()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->BNode);
  }

  /** \fn void setA(SiconosMatrix *m)
  *   \brief allows to save the A of the LinearSystemDSXML
  *   \return The A SiconosMatrix to save
  */
  inline void setA(SiconosMatrix *m)
  {
    if (this->ANode != NULL)
    {
      //SiconosDOMTreeTools::setSiconosMatrixValue(this->ANode, *m);
      SiconosDOMTreeTools::setSiconosMatrixValue(this->ANode, m);
    }
    else
    {
      this->ANode = SiconosDOMTreeTools::createMatrixNode(this->rootDSXMLNode, LSDS_A, m);
    }
  }

  /** \fn void setB(SiconosMatrix *m)
  *   \brief allows to save the B of the LinearSystemDSXML
  *   \return The B SiconosMatrix to save
  */
  inline void setB(SiconosMatrix *m)
  {
    if (this->BNode != NULL)
    {
      //SiconosDOMTreeTools::setSiconosMatrixValue(this->BNode, *m);
      SiconosDOMTreeTools::setSiconosMatrixValue(this->BNode, m);
    }
    else this->BNode = SiconosDOMTreeTools::createMatrixNode(this->rootDSXMLNode, LSDS_B, m);
  }
  /////////////////////////////

  /** \fn inline string getUPlugin()
  *   \brief Return the u Plugin name of the LinearSystemDSXML
  *   \return The u Plugin name of the LinearSystemDSXML
  *  \exception XMLException
  */
  inline string getUPlugin()
  {
    if (this->isUPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->uNode, LSDS_VECTORPLUGIN);
    XMLException::selfThrow("LinearSystemDSXML - getUPlugin : u is not calculated from a plugin ; u vector is given");
  }

  /** \fn inline SimpleVector getUVector()
  *   \brief Return u vector of the LinearSystemDSXML
  *   \return SimpleVector : u of LinearSystemDSXML
  *  \exception XMLException
  */
  inline /*SiconosVector*/SimpleVector getUVector()
  {
    if (this->isUPlugin())
      XMLException::selfThrow("LinearSystemDSXML - getUVector : u vector is not given ; u is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(this->uNode);
  }

  /** \fn inline void setUVector(SiconosVector *v)
  *   \brief allows to save the u vector of the LinearSystemDSXML
  *   \param SiconosVector *u : SiconosVector U to save
  */
  inline void setUVector(SiconosVector *v)
  {
    if (this->uNode != NULL)
    {
      //SiconosDOMTreeTools::setSiconosVectorValue(this->uNode, *v);
      SiconosDOMTreeTools::setSiconosVectorValue(this->uNode, v);
    }
    else this->uNode = SiconosDOMTreeTools::createVectorNode(this->rootDSXMLNode, LSDS_U, v);
  }

  /** \fn inline string getFPlugin()
  *   \brief Return the f Plugin name of the LinearSystemDSXML
  *   \return The f Plugin name of the LinearSystemDSXML
  *  \exception XMLException
  */
  inline string getFPlugin()
  {
    if (this->isFPlugin())
      return  SiconosDOMTreeTools::getStringAttributeValue(this->fNode, LSDS_VECTORPLUGIN);
    XMLException::selfThrow("LinearSystemDSXML - getUPlugin : f is not calculated from a plugin ; f vector is given");
  }

  /** \fn inline SimpleVector getFVector()
  *   \brief Return f vector of the LinearSystemDSXML
  *   \return SimpleVector : value of f of LinearSystemDSXML
  *  \exception XMLException
  */
  inline /*SiconosVector*/ SimpleVector getFVector()
  {
    if (this->isFPlugin())
      XMLException::selfThrow("LinearSystemDSXML - getFVector : f vector is not given ; f is calculated from a plugin");

    return  SiconosDOMTreeTools::getSiconosVectorValue(this->fNode);
  }

  /** \fn inline void setFVector(SiconosVector *v)
  *   \brief allows to save the f vector of the LinearSystemDSXML
  *   \return The f SimpleVector to save
  */
  inline void setFVector(SiconosVector *v)
  {
    if (this->fNode != NULL)
    {
      //SiconosDOMTreeTools::setSiconosVectorValue(this->fNode, *v);
      SiconosDOMTreeTools::setSiconosVectorValue(this->fNode, v);
    }
    else this->fNode = SiconosDOMTreeTools::createVectorNode(this->rootDSXMLNode, LSDS_F, v);
  }


  /** \fn bool isUPlugin()
  *   \brief Return true if u is calculated from a plugin
  *   \return True if u is calculated from plugin
  */
  inline bool isUPlugin()
  {
    return xmlHasProp((xmlNodePtr)uNode, (xmlChar *) LSDS_VECTORPLUGIN.c_str());
  }

  /** \fn bool isFPlugin()
  *   \brief Return true if f is calculated from a plugin
  *   \return True if f is calculated from plugin
  */
  inline bool isFPlugin()
  {
    return xmlHasProp((xmlNodePtr)fNode, (xmlChar *) LSDS_VECTORPLUGIN.c_str());
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
//$Log: LinearSystemDSXML.h,v $
//Revision 1.16  2004/09/27 08:24:26  charlety
//
//_ Modifications in doxygen comments.
//
//Revision 1.15  2004/09/10 11:26:27  charlety
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
//Revision 1.14  2004/08/23 14:30:02  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.13  2004/08/20 15:26:45  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NSDS and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.12  2004/07/29 14:25:43  jbarbier
