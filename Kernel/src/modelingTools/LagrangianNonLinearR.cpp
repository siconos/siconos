
#include "LagrangianNonLinearR.h"
using namespace std;



// Default constructor
LagrangianNonLinearR::LagrangianNonLinearR(): Relation()
{
  relationType = LAGRANGIANNONLINEARRELATION;
}

// xml constructor
LagrangianNonLinearR::LagrangianNonLinearR(RelationXML* relxml, Interaction* inter):
  Relation(relxml, inter), computeJacobianPtr(NULL), computeHPtr(NULL)
{
  relationType = LAGRANGIANNONLINEARRELATION;
}

// constructor from a set of data
LagrangianNonLinearR::LagrangianNonLinearR(const string& computeInput, const string& computeOutput, Interaction* inter):
  Relation(inter), computeJacobianPtr(NULL), computeHPtr(NULL)
{
  relationType = LAGRANGIANNONLINEARRELATION;
  // computeInput
  setComputeInputFunction(cShared.getPluginName(computeInput), cShared.getPluginFunctionName(computeInput));
  // computeOutput
  setComputeOutputFunction(cShared.getPluginName(computeOutput), cShared.getPluginFunctionName(computeOutput));
}

// copy constructor (inter is optional)
LagrangianNonLinearR::LagrangianNonLinearR(const Relation & newLNLR, Interaction* inter):
  Relation(newLNLR, inter)
{
  if (relationType !=  LAGRANGIANNONLINEARRELATION)
    RuntimeException::selfThrow("LagrangianNonLinearR:: copy constructor, inconsistent relation types for copy");

  //const LagrangianNonLinearR * lnlr = static_cast<const LagrangianNonLinearR*>(&newLNLR);
  // \todo

}


LagrangianNonLinearR::~LagrangianNonLinearR()
{}

void LagrangianNonLinearR::computeJacobian()
{
  if (computeJacobianPtr == NULL)
    RuntimeException::selfThrow("computeJacobian() is not linked to a plugin function");
  computeJacobianPtr(NULL, NULL, NULL, NULL);
}

void LagrangianNonLinearR::saveRelationToXML()
{
  IN("LagrangianNonLinearR::saveRelationToXML\n");
  if (relationxml != NULL)
  {
    relationxml->setComputeInputPlugin(computeInputName);
    relationxml->setComputeOutputPlugin(computeOutputName);
  }
  else RuntimeException::selfThrow("LagrangianNonLinearR::saveRelationToXML - object RelationXML does not exist");
  OUT("LagrangianNonLinearR::saveRelationToXML\n");
}

LagrangianNonLinearR* LagrangianNonLinearR::convert(Relation *r)
{
  cout << "LagrangianNonLinearR::convert (Relation *r)" << endl;
  LagrangianNonLinearR* lnlr = dynamic_cast<LagrangianNonLinearR*>(r);
  return lnlr;
}

