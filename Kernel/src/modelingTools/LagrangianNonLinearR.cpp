
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


