
#include "LsodarXML.h"

LsodarXML::LsodarXML() : OneStepIntegratorXML()
{}

LsodarXML::LsodarXML(xmlNode * LsodarNode,  map<int, bool> definedDSNumbers)
  : OneStepIntegratorXML(LsodarNode, definedDSNumbers)
{}

//$Log: LsodarXML.cpp,v $
//Revision 1.5  2004/09/23 14:09:24  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.4  2004/09/15 13:23:13  jbarbier
//- corrections in the OneStepNSProblem, for the XML save. The list of interaction
//linked to the onestepnsproblem is now saved correctly. It is updated before
//during the creation process.
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
//Revision 1.7  2004/07/29 14:25:42  jbarbier
