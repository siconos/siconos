
/** \class LinearBCXML
*   \brief This class manages Linear BC data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/25/2004
*
*
*
* LinearBCXML allows to manage data of a LinearBC DOM tree.
*/


#ifndef __LINEARBCXML__
#define __LINEARBCXML__


#include <libxml/tree.h>

//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include "SiconosDOMTreeTools.h"

#include "BoundaryConditionXML.h"


//using namespace std;

//Tags
const string LINEARBC_OMEGA = "Omega";
const string LINEARBC_OMEGA0 = "Omega0";
const string LINEARBC_OMEGAT = "OmegaT";


class LinearBCXML : public BoundaryConditionXML
{
public:

  LinearBCXML();

  /** \fn LinearBCXML(xmlNode * LinearBCNode)
  *   \brief Build a LinearBCXML object from a DOM tree describing a LinearBC
  *   \param xmlNode * LinearBCNode : the LinearBC DOM tree
  */
  LinearBCXML(xmlNode * LinearBCNode);

  ~LinearBCXML();

  /** \fn SimpleVector getOmega()
  *   \brief Return Omega of the LinearBCXML
  *   \return SimpleVector : Omega of LinearBCXML
  */
  inline /*SiconosVector*/SimpleVector getOmega()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->omegaNode);
  }

  /** \fn SiconosMatrix getOmega0()
  *   \brief Return the Omega0 of the LinearBCXML
  *   \return The Omega0 SiconosMatrix of the LinearBCXML
  */
  inline SiconosMatrix getOmega0()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->omega0Node);
  }

  /** \fn SiconosMatrix getOmegaT()
  *   \brief Return the OmegaT of the LinearBCXML
  *   \return The OmegaT SiconosMatrix of the LinearBCXML
  */
  inline SiconosMatrix getOmegaT()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->omegaTNode);
  }

  /** \fn void setOmega(SiconosVector *v)
  *   \brief allows to save the Omega of the LinearBCXML
  *   \param The Omega SiconosVector to save
  */
  inline void setOmega(SiconosVector *v)
  {
    if (this->omegaNode == NULL)
    {
      this->omegaNode = SiconosDOMTreeTools::createVectorNode(this->rootBCNode, LINEARBC_OMEGA, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->omegaNode, v);
  }

  /** \fn void setOmega0(SiconosMatrix *m)
  *   \brief allows to save the Omega0 of the LinearBCXML
  *   \param The Omega0 SiconosMatrix to save
  */
  inline void setOmega0(SiconosMatrix *m)
  {
    if (this->omega0Node == NULL)
    {
      this->omega0Node = SiconosDOMTreeTools::createMatrixNode(this->rootBCNode, LINEARBC_OMEGA0, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->omega0Node, m);
  }

  /** \fn void setOmegaT(SiconosMatrix *m)
  *   \brief allows to save the OmegaT of the LinearBCXML
  *   \param The OmegaT SiconosMatrix to save
  */
  inline void setOmegaT(SiconosMatrix *m)
  {
    if (this->omegaTNode == NULL)
    {
      this->omegaTNode = SiconosDOMTreeTools::createMatrixNode(this->rootBCNode, LINEARBC_OMEGAT, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->omegaTNode, m);
  }


  /** \fn void updateBoundaryConditionXML( xmlNode* node)//, BoundaryCondition* bc )
   *  \brief makes the operations to add a BoundaryCondition to the DynamicalSystemXML
   *  \param xmlNode* : the root node of this BoundaryCondition
  //     *  \param BoundaryCondition* : the BoundaryCondition of the DS
   */
  void updateBoundaryConditionXML(xmlNode* node/*, BoundaryCondition* bc*/);


private:

  //Nodes
  xmlNode * omegaNode;
  xmlNode * omega0Node;
  xmlNode * omegaTNode;

  //Methods
  /** \fn loadLinearBCProperties(xmlNode * LinearBCnode)
  *   \brief load the different properties of a LinearBC
  *   \exception XMLException : if a property of the LinearBC lacks in the DOM tree
  */
  void loadLinearBCProperties();
};


#endif
//$Log: LinearBCXML.h,v $
//Revision 1.10  2004/09/27 08:24:26  charlety
//
//_ Modifications in doxygen comments.
//
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
//Revision 1.8  2004/09/10 08:04:50  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.7  2004/07/29 14:25:43  jbarbier
