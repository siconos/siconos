//$Id: LsodarXML.h,v 1.4 2004/09/23 14:09:24 jbarbier Exp $

/** \class LsodarXML
*   \brief This class manages Lsodar data part
*   \author J. Blanc-Tranchant
*   \version 1.0
*   \date 05/17/2004
*
*
* $Date: 2004/09/23 14:09:24 $
* $Revision: 1.4 $
* $Author: jbarbier $
* $Source: /CVS/Siconos/SICONOS/src/xml/LsodarXML.h,v $
*
* LsodarXML allows to manage data of a Lsodar DOM tree.
*/


#ifndef __LsodarXMLDEF__
#define __LsodarXMLDEF__


#include <libxml/tree.h>
#include "OneStepIntegratorXML.h"


using namespace std;

const string LSODAR_R = "r";


class LsodarXML : public OneStepIntegratorXML
{
public:

  LsodarXML();

  /** \fn LsodarXML(xmlNode * LsodarNode)
  *   \brief Build a LsodarXML object from a DOM tree describing Lsodar OneStepIntegrator
  *   \param LsodarNode : the Lsodar DOM tree
  *   \param map<int, bool> definedDSNumbers : to know if DS numbers are not used by another OneStepIntegrator
  */
  LsodarXML(xmlNode * LsodarNode,  map<int, bool> definedDSNumbers);


private:

  //Nodes
};


#endif
//$Log: LsodarXML.h,v $
//Revision 1.4  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.3  2004/09/14 13:49:58  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.2  2004/08/09 15:00:55  jbarbier
//- changes in the cardinality of some attributes of the DynamicalSystem,
//OneStepIntegrator
//
//- modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//
//- corrections in the test xml files
//
//Revision 1.1  2004/08/05 14:35:59  charlety
//
//_ LSODAR --> Lsodar (chapter 2)
//
//Revision 1.6  2004/07/29 14:25:42  jbarbier
//- $Log: LsodarXML.h,v $
//- Revision 1.4  2004/09/23 14:09:24  jbarbier
//- - modification of the integrators, the attribute r is always optional.
//-
//- - modification of the LagrangianNonLinearR. computeInput and computeOutput are
//- required.
//-
//- Revision 1.3  2004/09/14 13:49:58  jbarbier
//- - files added in sample/ to run run the main_siconos test program
//-
//- - all the platform can now be saved in an XML file when it is created manually
//-
//- Revision 1.2  2004/08/09 15:00:55  jbarbier
//- - changes in the cardinality of some attributes of the DynamicalSystem,
//- OneStepIntegrator
//-
//- - modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//-
//- - corrections in the test xml files
//-
//- Revision 1.1  2004/08/05 14:35:59  charlety
//-
//- _ LSODAR --> Lsodar (chapter 2)
//- and $Id: LsodarXML.h,v 1.4 2004/09/23 14:09:24 jbarbier Exp $ added
//
