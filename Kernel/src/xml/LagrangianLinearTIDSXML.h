/** \class LagrangianLinearTIDSXML
*   \brief This class manages Lagrangian TIDS data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/11/2004
*
*
*
* LagrangianLinearTIDSXML allows to manage data of a LagrangianLinearTIDS DOM tree.
*/


#ifndef __LAGRANGIANTIDSXML__
#define __LAGRANGIANTIDSXML__

#include "LagrangianDSXML.h"

const std::string LTIDS_K = "K";
const std::string LTIDS_C = "C";


class LagrangianLinearTIDSXML : public LagrangianDSXML
{
public:
  LagrangianLinearTIDSXML();

  /** \fn LagrangianLinearTIDSXML(xmlNode * LagrangianLinearTIDSNode, int number)
  *   \brief Build a LagrangianLinearTIDSXML object from a DOM tree describing a LagrangianLinearTIDS
  *   \param LagrangianLinearTIDSNode : the LagrangianLinearTIDS DOM tree
  *   \param bool isBVP : if NSDS is BVP LagrangianLinearTIDS have boundary condition
  */
  LagrangianLinearTIDSXML(xmlNode * LagrangianLinearTIDSNode, bool isBVP);

  /** \fn SiconosMatrix getK()
  *   \brief Return the K of the LagrangianLinearTIDSXML
  *   \return The K SiconosMatrix of the LagrangianLinearTIDSXML
  */
  inline SiconosMatrix getK()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->KNode);
  }

  /** \fn void setK(SiconosMatrix *m)
  *   \brief allows to save the K of the LagrangianLinearTIDSXML
  *   \return The K SiconosMatrix to save
  */
  inline void setK(SiconosMatrix *m)
  {
    if (this->KNode == NULL)
    {
      this->KNode = SiconosDOMTreeTools::createMatrixNode(this->rootDSXMLNode, LTIDS_K, *m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(this->KNode, *m);
  }

  /** \fn SiconosMatrix getC()
  *   \brief Return the C of the LagrangianLinearTIDSXML
  *   \return The C SiconosMatrix of the LagrangianLinearTIDSXML
  */
  inline SiconosMatrix getC()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->CNode);
  }

  /** \fn void setC(SiconosMatrix *m)
  *   \brief allows to save the C of the LagrangianLinearTIDSXML
  *   \return The C SiconosMatrix to save
  */
  inline void setC(SiconosMatrix *m)
  {
    if (this->CNode == NULL)
    {
      this->CNode = SiconosDOMTreeTools::createMatrixNode(this->rootDSXMLNode, LTIDS_C, *m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(this->CNode, *m);
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

  /** \fn loadLagrangianLinearTIDSProperties( bool )
  *   \brief load the different properties of a LagrangianLinearTIDS
  *   \param bool : value which determines if the Dynamical System is BVP or not
  *   \exception XMLException : if a property of the LagrangianLinearTIDS lacks in the DOM tree
  */
  void loadLagrangianLinearTIDSProperties(bool);

};

#endif
