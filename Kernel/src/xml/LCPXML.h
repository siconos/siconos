//$Id: LCPXML.h,v 1.15 2004/09/27 13:27:14 jbarbier Exp $
/** \class LCPXML
*   \brief This class manages Lagrangian LCP data
*   \author J. Blanc-Tranchant
*   \version 1.0
*   \date 05/18/2004
*
* $Date: 2004/09/27 13:27:14 $
* $Revision: 1.15 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/LCPXML.h,v $
*
*
* LCPXML allows to manage data of a LCP DOM tree.
*/


#ifndef __LCPXML__
#define __LCPXML__



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

const string LCP_M = "M";
const string LCP_Q = "q";


class LCPXML : public OneStepNSProblemXML
{
public:
  LCPXML();

  /** \fn LCPXML(xmlNode * LCPNode)
  *   \brief Build a LCPXML object from a DOM tree describing a LCP
  *   \param LCPNode : the LCP DOM tree
  *   \param vector<int> definedInteractionNumbers : the Interaction numbers effectivly defined in the model
  *   \exception XMLException : if a property of the LCP lacks in the DOM tree
  */
  LCPXML(xmlNode * LCPNode, vector<int> definedInteractionNumbers);

  /** \fn SiconosMatrix getM()
  *   \brief Return M
  *   \return The M SiconosMatrix of the LCP
  */
  inline SiconosMatrix getM()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->MNode);
  }

  /** \fn SimpleVector getQ()
  *   \brief Return vector q
  *   \return SimpleVector : q vector of the LCP
  */
  inline /*SiconosVector*/SimpleVector getQ()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(this->qNode);
  }

  /** \fn void setM(SiconosMatrix *m)
  *   \brief allows to save M
  *   \return The M SiconosMatrix to save
  */
  inline void setM(SiconosMatrix *m)
  {
    if (this->hasM() == false)
    {
      this->MNode = SiconosDOMTreeTools::createMatrixNode(this->rootNSProblemXMLNode, LCP_M, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixValue(this->MNode, m);
  }

  /** \fn void setQ(SiconosVector *v)
  *   \brief allows to save q
  *   \return The q SiconosVector to save
  */
  inline void setQ(SiconosVector *v)
  {
    if (this->hasQ() == false)
    {
      this->qNode = SiconosDOMTreeTools::createVectorNode(this->rootNSProblemXMLNode, LCP_Q, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorValue(this->qNode, v);
  }

  /** \fn bool hasM()
   *  \brief returns true if MNode is defined
   *  \return true if MNode is defined
   */
  inline bool hasM()
  {
    return (this->MNode != NULL);
  }

  /** \fn bool hasQ()
   *  \brief returns true if qNode is defined
   *  \return true if qNode is defined
   */
  inline bool hasQ()
  {
    return (this->qNode != NULL);
  }

  /** \fn void updateOneStepNSProblemXML( xmlNode* node, OneStepNSProblemXML* str )
  *   \brief makes the operations to create a OneStepNSProblemXML to the StrategyXML
  *   \param xmlNode* : the root node of the OneStepNSProblemXML
  *   \param OneStepNSProblem* : the OneStepNSProblem of this OneStepNSProblemXML
  */
  void updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb);


private:

  //Nodes
  xmlNode * MNode;
  xmlNode * qNode;
};


#endif
//$Log: LCPXML.h,v $
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
//Revision 1.13  2004/09/14 13:49:57  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.12  2004/09/10 11:26:26  charlety
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
//Revision 1.11  2004/07/29 14:25:42  jbarbier
//- $Log: LCPXML.h,v $
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
//- Revision 1.13  2004/09/14 13:49:57  jbarbier
//- - files added in sample/ to run run the main_siconos test program
//-
//- - all the platform can now be saved in an XML file when it is created manually
//-
//- Revision 1.12  2004/09/10 11:26:26  charlety
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
//- and $Id: LCPXML.h,v 1.15 2004/09/27 13:27:14 jbarbier Exp $ added
//
//Revision 1.10  2004/07/29 14:04:00  jbarbier
//- new test on SiconosMemoryXML
//
//- last functions hasAttribute() in the XML part added
//
//Revision 1.9  2004/07/07 08:14:54  jbarbier
//-modifications on the test after renamming
//
//-modification of the XML schema, attributs row, col and size of Matrices and
//Vector changed from 'positiveInteger' to 'nonNegativeInteger'
//
//-function setSiconosVector/Matrix which take a SiconosVector/Matrix* in parameter to avoid
//unnecessary vector and matrix copies
//
//-CVS $id and $log added in the xml files
//
//Revision 1.8  2004/06/29 15:12:02  acary
//Change in the naming comvention for the LCP
//The LCP Matrix is now denoted by M.
//The LCP Vector is now denoted by q.
//
