
/** \class LagrangianLinearRXML
*   \brief This class manages LagrangianLinear Relation data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/25/2004
*
*
*
* LagrangianLinearRXML allows to manage data of a LLRelation DOM tree.
*/


#ifndef __LLRelationXML__
#define __LLRelationXML__


#include <libxml/tree.h>

#include "RelationXML.h"
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosDOMTreeTools.h"


using namespace std;


const string LLR_H = "H";
const string LLR_B = "b";
//#include "XMLTagsName.h"


class LagrangianLinearRXML : public RelationXML
{
public:

  LagrangianLinearRXML();

  /** \fn LagrangianLinearRXML(xmlNode * LLRelationNode)
  *   \brief Build a LagrangianLinearRXML object from a DOM tree describing a Relation with LL type
  *   \param LagrangianLinearRXML : the LagrangianLinearR DOM tree
  *   \exception XMLException : if a property of the LagrangianLinear Relation lacks in the DOM tree
  */
  LagrangianLinearRXML(xmlNode * LLRelationNode);

  ~LagrangianLinearRXML();

  /** \fn SiconosMatrix getH()
  *   \brief Return the H of the LLRelationXML
  *   \return The H SiconosMatrix of the LLRelationXML
  */
  inline SiconosMatrix getH()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->HNode);
  }


  /** \fn SimpleVector getB()
  *   \brief Return b vector of the LLRelationXML
  *   \return SimpleVector : b vector of the LLRelationXML
  */
  inline /*SiconosVector*/SimpleVector getB()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->bNode);
  }

  /** \fn void setH(SiconosMatrix *matrix)
  *   \brief Change the H matrix value (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for H matrix
  */
  void setH(SiconosMatrix *matrix);

  /** \fn void setB(SiconosVector *vector)
  *   \brief Change the b vector value (in xml file or external data file switch his origin position)
  *   \param SiconosVector vector : the new value for b vector
  */
  void setB(SiconosVector *vector);


private:


  //Nodes
  xmlNode * HNode;
  xmlNode * bNode;

};


#endif
//$Log: LagrangianLinearRXML.h,v $
//Revision 1.11  2005/03/08 12:41:38  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.10  2005/02/24 15:50:21  jbarbier
//- LCP prepared to changes needed for several interactions
//
//- new function for the SiconosMatrices to copy a block matrix into another matrix
//
//- tag names of BoundaryConditionXML, DSInputOutputXML, DSXML, InteractionXML, LagrangianLinearRXML, LagrangianNLDSXML put in XMLTagNames.h
//
//Revision 1.9  2005/01/18 10:35:17  jbarbier
//- attribute "r" no longer used for Moreau integrator
//
//- modificatoin in the tests for Moreau integrator
//
//- file XMLTagsName.h for further use to regroup all xml tags name...
//
//Revision 1.8  2004/09/27 08:24:26  charlety
//
//_ Modifications in doxygen comments.
//
//Revision 1.7  2004/09/14 13:49:58  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.6  2004/09/10 11:26:27  charlety
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
//Revision 1.5  2004/07/29 14:25:42  jbarbier
