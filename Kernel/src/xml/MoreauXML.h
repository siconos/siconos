
/** \class MoreauXML
*   \brief This class manages Moreau data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/17/2004
*
*
* MoreauXML allows to manage data of a Moreau DOM tree.
*/

#ifndef __MOREAUXML__
#define __MOREAUXML__


#include <libxml/tree.h>

#include "OneStepIntegratorXML.h"


using namespace std;


const string MOREAU_R = "r";
const string MOREAU_W = "W";
const string MOREAU_THETA = "Theta";

class MoreauXML : public OneStepIntegratorXML
{
public:
  MoreauXML();

  /** \fn MoreauXML(xmlNode * MoreauNode)
  *   \brief Build a MoreauXML object from a DOM tree describing Moreau OneStepIntegrator
  *   \param xmlNode * MoreauNode : the Moreau DOM tree
  *   \param map<int, bool> definedDSNumbers : to know if DS numbers are not used by another OneStepIntegrator
  *   \exception XMLException : if the W property of the Moreau lacks in the DOM tree
  */
  MoreauXML(xmlNode * MoreauNode, map<int, bool> definedDSNumbers);

  /** \fn SiconosMatrix getW()
  *   \brief Return the W of the MoreauXML
  *   \return The W SiconosMatrix of the MoreauXML
  */

  ~MoreauXML();


  /** \fn bool hasW()
   *  \brief return true if wNode is defined
   *  \return true if wNode is defined
   */
  inline bool hasW()
  {
    return (this->WNode != NULL);
  }

  /** \fn SiconosMatrix getW()
  *   \brief Return the w of the OneStepIntegratorXML
  *   \return SiconosMatrix : the w of the OneStepIntegratorXML
  */
  inline SiconosMatrix getW()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->WNode);
  }

  /** \fn void setW(SiconosMatrix *m)
  *   \brief allows to save the w of the OneStepIntegratorXML
  *   \param SiconosMatrix* : the w to save
  */
  inline void setW(SiconosMatrix *m)
  {
    if (this->hasW() == false)
    {
      this->WNode = SiconosDOMTreeTools::createMatrixNode(this->rootIntegratorXMLNode, MOREAU_W, m);
    }
    //else SiconosDOMTreeTools::setSiconosMatrixValue(this->WNode, *m);
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->WNode, m);
  }

  /** \fn bool hasTheta()
   *  \brief return true if ThetaNode is defined
   *  \return true if ThetaNode is defined
   */
  inline bool hasTheta()
  {
    return (this->ThetaNode != NULL);
  }

  /** \fn SiconosMatrix getTheta()
  *   \brief Return the theta of the OneStepIntegratorXML
  *   \return SiconosMatrix : the theta of the OneStepIntegratorXML
  */
  inline double getTheta()
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(this->ThetaNode);
  }

  /** \fn void setTheta(double t)
  *   \brief allows to save  Theta of the OneStepIntegratorXML
  *   \param double t : the Theta to save
  */
  inline void setTheta(double t)
  {
    if (this->hasTheta() == false)
    {
      this->ThetaNode = SiconosDOMTreeTools::createDoubleNode(this->rootIntegratorXMLNode, MOREAU_THETA, t);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->ThetaNode, t);
  }


private:

  //Nodes
  xmlNode * WNode;
  xmlNode * ThetaNode;


};


#endif
//$Log: MoreauXML.h,v $
//Revision 1.19  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.18  2004/09/14 13:49:58  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.17  2004/08/09 15:00:55  jbarbier
//- changes in the cardinality of some attributes of the DynamicalSystem,
//OneStepIntegrator
//
//- modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//
//- corrections in the test xml files
//
//Revision 1.16  2004/07/12 13:04:34  jbarbier
