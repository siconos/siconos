//$Id: LagrangianTIDSXML.h,v 1.9 2004/09/10 11:26:27 charlety Exp $

/** \class LagrangianTIDSXML
*   \brief This class manages Lagrangian TIDS data
*   \author J. Blanc-Tranchant
*   \version 1.0
*   \date 05/11/2004
*
*
* $Date: 2004/09/10 11:26:27 $
* $Revision: 1.9 $
* $Author: charlety $
* $Source: /CVS/Siconos/SICONOS/src/xml/LagrangianTIDSXML.h,v $
*
* LagrangianTIDSXML allows to manage data of a LagrangianTIDS DOM tree.
*/


#ifndef __LAGRANGIANTIDSXML__
#define __LAGRANGIANTIDSXML__


#include <vector>
#include <string>
#include <libxml/tree.h>

#include "LagrangianNLDSXML.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SiconosMatrix.h"

#include "SiconosDOMTreeTools.h"


using namespace std;


const string LTIDS_K = "K";
const string LTIDS_C = "C";


class LagrangianTIDSXML : public LagrangianNLDSXML
{
public:
  LagrangianTIDSXML();

  /** \fn LagrangianTIDSXML(xmlNode * LagrangianTIDSNode, int number)
  *   \brief Build a LagrangianTIDSXML object from a DOM tree describing a LagrangianTIDS
  *   \param LagrangianTIDSNode : the LagrangianTIDS DOM tree
  *   \param bool isBVP : if NSDS is BVP LagrangianTIDS have boundary condition
  */
  LagrangianTIDSXML(xmlNode * LagrangianTIDSNode, bool isBVP);

  /** \fn SiconosMatrix getK()
  *   \brief Return the K of the LagrangianTIDSXML
  *   \return The K SiconosMatrix of the LagrangianTIDSXML
  */
  inline SiconosMatrix getK()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->KNode);
  }

  /** \fn void setK(SiconosMatrix *m)
  *   \brief allows to save the K of the LagrangianTIDSXML
  *   \return The K SiconosMatrix to save
  */
  inline void setK(SiconosMatrix *m)
  {
    if (this->KNode == NULL)
    {
      this->KNode = SiconosDOMTreeTools::createMatrixNode(this->rootDSXMLNode, LTIDS_K, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->KNode, m);
  }

  /** \fn SiconosMatrix getC()
  *   \brief Return the C of the LagrangianTIDSXML
  *   \return The C SiconosMatrix of the LagrangianTIDSXML
  */
  inline SiconosMatrix getC()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->CNode);
  }

  /** \fn void setC(SiconosMatrix *m)
  *   \brief allows to save the C of the LagrangianTIDSXML
  *   \return The C SiconosMatrix to save
  */
  inline void setC(SiconosMatrix *m)
  {
    if (this->CNode == NULL)
    {
      this->CNode = SiconosDOMTreeTools::createMatrixNode(this->rootDSXMLNode, LTIDS_C, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->CNode, m);
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

  xmlNode * KNode;
  xmlNode * CNode;


  //Methods

  /** \fn loadLagrangianTIDSProperties( bool )
  *   \brief load the different properties of a LagrangianTIDS
  *   \param bool : value which determines if the Dynamical System is BVP or not
  *   \exception XMLException : if a property of the LagrangianTIDS lacks in the DOM tree
  */
  void loadLagrangianTIDSProperties(bool);

};

#endif
//$Log: LagrangianTIDSXML.h,v $
//Revision 1.9  2004/09/10 11:26:27  charlety
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
//Revision 1.8  2004/08/23 14:30:02  jbarbier
//- All the dynamical systems can be created in a comand program and added to a
//NSDS. The save is OK, but the creation of the boundary conditions is not yet
//finished.
//
//Revision 1.7  2004/07/29 14:25:43  jbarbier
//- $Log: LagrangianTIDSXML.h,v $
//- Revision 1.9  2004/09/10 11:26:27  charlety
//-
//- _ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//-
//- _ All the tests which worked with the previous version of the vector are OK with the new version.
//-
//- _ Example SICONOS and bouncingBall are OK
//-
//- _ some comments have still to be adapted to NewSiconosVector .
//-
//- _ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//-
//- Revision 1.8  2004/08/23 14:30:02  jbarbier
//- - All the dynamical systems can be created in a comand program and added to a
//- NSDS. The save is OK, but the creation of the boundary conditions is not yet
//- finished.
//- and $Id: LagrangianTIDSXML.h,v 1.9 2004/09/10 11:26:27 charlety Exp $ added
//
