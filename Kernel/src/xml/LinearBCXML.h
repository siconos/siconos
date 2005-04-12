
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
