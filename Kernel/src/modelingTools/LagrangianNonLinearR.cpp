
#include "LagrangianNonLinearR.h"
#include "check.h"

LagrangianNonLinearR::LagrangianNonLinearR()
{
  this->relationType = LAGRANGIANNONLINEARRELATION; //"LagrangianNonLinearR";
}
LagrangianNonLinearR::~LagrangianNonLinearR()
{}


void LagrangianNonLinearR::computeJacobian(void)
{
  if (computeJacobianPtr == NULL)
    RuntimeException::selfThrow("computeJacobian() is not linked to a plugin function");
  this->computeJacobianPtr(NULL, NULL, NULL, NULL);
}

void LagrangianNonLinearR::saveRelationToXML()
{
  IN("LagrangianNonLinearR::saveRelationToXML\n");
  Relation::saveRelationToXML();
  if (this->relationxml != NULL)
  {
    this->relationxml->setComputeInputPlugin(this->computeInputName);
    this->relationxml->setComputeOutputPlugin(this->computeOutputName);
  }
  else RuntimeException::selfThrow("LagrangianNonLinearR::saveRelationToXML - object RelationXML does not exist");
  OUT("LagrangianNonLinearR::saveRelationToXML\n");
}

void LagrangianNonLinearR::createRelation(LagrangianNonLinearRXML * relationXML,
    string computeInput, string computeOutput)//, Interaction * interaction)
{
  if (relationXML != NULL)
  {
    //    this->init();
    this->relationxml = relationXML;
    this->relationType = LAGRANGIANNONLINEARRELATION; //"LagrangianNonLinearR";
    this->fillRelationWithRelationXML();
  }
  else
  {
    this->relationxml = relationXML;

    // computeInput
    this->setComputeInputFunction(this->cShared.getPluginName(computeInput), this->cShared.getPluginFunctionName(computeInput));

    // computeOutput
    this->setComputeOutputFunction(this->cShared.getPluginName(computeOutput), this->cShared.getPluginFunctionName(computeOutput));
  }
}


LagrangianNonLinearR* LagrangianNonLinearR::convert(Relation *r)
{
  cout << "LagrangianNonLinearR::convert (Relation *r)" << endl;
  LagrangianNonLinearR* lnlr = dynamic_cast<LagrangianNonLinearR*>(r);
  return lnlr;
}


//$Log: LagrangianNonLinearR.cpp,v $
//Revision 1.13  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.12  2005/03/07 13:17:19  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.11  2005/02/01 11:08:41  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.10  2005/01/31 16:26:20  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.9  2004/12/08 12:49:37  jbarbier
//- changes in the XML Schema, respect of the recommandations of the W3C
//version 1.1
//
//- changes in all balises DS, Relation, NSLaw, OneStepIntegrator, OneStepNSProblem
//in the XML files into specific names like LagrangianDS, LinearSystemDS, ...
//for the DS
//
//Revision 1.8  2004/09/23 14:45:06  charlety
//
//_ Added a header file to main_siconos.cpp
//_ modified plugin functions signatures in model formalisation
//
//Revision 1.7  2004/09/23 14:09:23  jbarbier
//- modification of the integrators, the attribute r is always optional.
//
//- modification of the LagrangianNonLinearR. computeInput and computeOutput are
//required.
//
//Revision 1.6  2004/09/03 14:41:42  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NSDS
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.5  2004/08/17 15:12:38  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.4  2004/08/12 11:55:14  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.3  2004/07/29 14:25:36  jbarbier
