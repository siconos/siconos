//$Id: QPXML.h,v 1.15 2004/09/27 13:27:14 jbarbier Exp $

/** \class QPXML
*   \brief This class manages Lagrangian QP data
*   \author J. Blanc-Tranchant
*   \version 1.0
*   \date 05/18/2004
*
*
* $Date: 2004/09/27 13:27:14 $
* $Revision: 1.15 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/QPXML.h,v $
*
* QPXML allows to manage data of a QP DOM tree.
*/


#ifndef __QPXMLDEF__
#define __QPXMLDEF__



#include <string>
#include <libxml/tree.h>

#include "OneStepNSProblemXML.h"
#include "OneStepNSProblem.h"

#include "SiconosDOMTreeTools.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"

using namespace std;

class OneStepNSProblem;

const string QP_Q = "Q";
const string QP_P = "p";


class QPXML : public OneStepNSProblemXML
{
public:
  QPXML();

  /** \fn QPXML(xmlNode * QPNode)
  *   \brief Build a QPXML object from a DOM tree describing a QP
  *   \param QPNode : the QP DOM tree
  *   \param vector<int> definedInteractionNumbers : the Interaction numbers effectivly defined in the model
  *   \exception XMLException : if a property of the QP lacks in the DOM tree
  */
  QPXML(xmlNode * QPNode, vector<int> definedInteractionNumbers);

  /** \fn SiconosMatrix getQ()
  *   \brief Return Q
  *   \return The Q SiconosMatrix of the QP
  */
  inline SiconosMatrix getQ()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->QNode);
  }

  /** \fn SimpleVector getP()
  *   \brief Return p
  *   \return SimpleVector :  vector p of the QP
  */
  inline /*SiconosVector*/SimpleVector getP()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(this->pNode);
  }

  /** \fn void setQ(SiconosMatrix *m)
  *   \brief allows  to save Q
  *   \param The Q SiconosMatrix to save
  */
  inline void setQ(SiconosMatrix *m)
  {
    if (this->hasQ() == false)
    {
      this->QNode = SiconosDOMTreeTools::createMatrixNode(this->rootNSProblemXMLNode, QP_Q, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->QNode, m);
  }

  /** \fn void setP(SiconosVector *v)
  *   \brief allows to save p
  *   \param SimpleVector* : vector p to save
  */
  inline void setP(SiconosVector *v)
  {
    if (this->hasP() == false)
    {
      this->pNode = SiconosDOMTreeTools::createVectorNode(this->rootNSProblemXMLNode, QP_P, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->pNode, v);
  }

  /** \fn bool hasP()
   *  \brief returns true if pNode is defined
   *  \return true if pNode is defined
   */
  inline bool hasP()
  {
    return (this->pNode != NULL);
  }

  /** \fn bool hasQ()
   *  \brief returns true if QNode is defined
   *  \return true if QNode is defined
   */
  inline bool hasQ()
  {
    return (this->QNode != NULL);
  }

  /** \fn void updateOneStepNSProblemXML( xmlNode* node, OneStepNSProblemXML* str )
  *   \brief makes the operations to create a OneStepNSProblemXML to the StrategyXML
  *   \param xmlNode* : the root node of the OneStepNSProblemXML
  *   \param OneStepNSProblem* : the OneStepNSProblem of this OneStepNSProblemXML
  */
  void updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb);


private:
  //Nodes
  xmlNode * QNode;
  xmlNode * pNode;

};

#endif
//$Log: QPXML.h,v $
//Revision 1.15  2004/09/27 13:27:14  jbarbier
//
//- Siconos schema renamed : SiconosModelSchema-V1.0.xsd
//
//- new required tags of the model : title, author, description, date, xmlSchema.
//They replace previous attributes author, description and date of the Model.
//
//Revision 1.14  2004/09/27 08:24:26  charlety
//
//_ Modifications in doxygen comments.
//
//Revision 1.13  2004/09/14 13:49:59  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.12  2004/09/10 11:26:29  charlety
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
//Revision 1.11  2004/07/29 14:25:45  jbarbier
//- $Log: QPXML.h,v $
//- Revision 1.15  2004/09/27 13:27:14  jbarbier
//-
//- - Siconos schema renamed : SiconosModelSchema-V1.0.xsd
//-
//- - new required tags of the model : title, author, description, date, xmlSchema.
//- They replace previous attributes author, description and date of the Model.
//-
//- Revision 1.14  2004/09/27 08:24:26  charlety
//-
//- _ Modifications in doxygen comments.
//-
//- Revision 1.13  2004/09/14 13:49:59  jbarbier
//- - files added in sample/ to run run the main_siconos test program
//-
//- - all the platform can now be saved in an XML file when it is created manually
//-
//- Revision 1.12  2004/09/10 11:26:29  charlety
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
//- and $Id: QPXML.h,v 1.15 2004/09/27 13:27:14 jbarbier Exp $ added
//
