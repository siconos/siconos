//$Id: AdamsXML.cpp,v 1.11 2004/09/23 14:09:24 jbarbier Exp $
#include "AdamsXML.h"

AdamsXML::AdamsXML() : OneStepIntegratorXML()
{}

AdamsXML::AdamsXML(xmlNode * AdamsNode,  map<int, bool> definedDSNumbersMap)
  : OneStepIntegratorXML(AdamsNode, definedDSNumbersMap)
{}

//$Log: AdamsXML.cpp,v $
//Revision 1.11  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.10  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
//
//Revision 1.9  2004/09/14 13:49:55  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.8  2004/08/09 15:00:54  jbarbier
//- changes in the cardinality of some attributes of the DynamicalSystem,
//OneStepIntegrator
//
//- modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//
//- corrections in the test xml files
//
//Revision 1.7  2004/07/29 14:25:41  jbarbier
//- $Log: AdamsXML.cpp,v $
//- Revision 1.11  2004/09/23 14:09:24  jbarbier
//- - modification of the integrators, the attribute r is always optional.
//-
//- - modification of the LagrangianNonLinearR. computeInput and computeOutput are
//- required.
//-
//- Revision 1.10  2004/09/15 13:23:13  jbarbier
//- - corrections in the OneStepNSProblem, for the XML save. The list of interaction
//- linked to the onestepnsproblem is now saved correctly. It is updated before
//- during the creation process.
//-
//- Revision 1.9  2004/09/14 13:49:55  jbarbier
//- - files added in sample/ to run run the main_siconos test program
//-
//- - all the platform can now be saved in an XML file when it is created manually
//-
//- Revision 1.8  2004/08/09 15:00:54  jbarbier
//- - changes in the cardinality of some attributes of the DynamicalSystem,
//- OneStepIntegrator
//-
//- - modifications in classes Moreau, Lsodar, Adams for these new cardinalities
//-
//- - corrections in the test xml files
//- and $Id: AdamsXML.cpp,v 1.11 2004/09/23 14:09:24 jbarbier Exp $ added
//
