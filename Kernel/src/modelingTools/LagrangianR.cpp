
#include "LagrangianR.h"
using namespace std;



// Default constructor with optional interaction parameter
LagrangianR::LagrangianR(Interaction* inter): Relation(inter)
{
  relationType = LAGRANGIANRELATION;
}

// xml constructor
LagrangianR::LagrangianR(RelationXML* relxml, Interaction* inter):
  Relation(relxml, inter), computeJacobianPtr(NULL), computeHPtr(NULL)
{
  relationType = LAGRANGIANRELATION;
}

// constructor from a set of data
LagrangianR::LagrangianR(const string& computeInput, const string& computeOutput, Interaction* inter):
  Relation(inter), computeJacobianPtr(NULL), computeHPtr(NULL)
{
  relationType = LAGRANGIANRELATION;
  // computeInput
  setComputeInputFunction(cShared.getPluginName(computeInput), cShared.getPluginFunctionName(computeInput));
  // computeOutput
  setComputeOutputFunction(cShared.getPluginName(computeOutput), cShared.getPluginFunctionName(computeOutput));
}

// copy constructor (inter is optional)
LagrangianR::LagrangianR(const Relation & newLNLR, Interaction* inter):
  Relation(newLNLR, inter)
{
  if (relationType !=  LAGRANGIANRELATION || relationType !=  LAGRANGIANLINEARRELATION)
    RuntimeException::selfThrow("LagrangianR:: copy constructor, inconsistent relation types for copy");

  //const LagrangianR * lnlr = static_cast<const LagrangianR*>(&newLNLR);
  // \todo

}


LagrangianR::~LagrangianR()
{}

void LagrangianR::computeJacobian()
{
  if (computeJacobianPtr == NULL)
    RuntimeException::selfThrow("computeJacobian() is not linked to a plugin function");
  computeJacobianPtr(NULL, NULL, NULL, NULL);
}

void LagrangianR::saveRelationToXML()
{
  IN("LagrangianR::saveRelationToXML\n");
  if (relationxml != NULL)
  {
    relationxml->setComputeInputPlugin(computeInputName);
    relationxml->setComputeOutputPlugin(computeOutputName);
  }
  else RuntimeException::selfThrow("LagrangianR::saveRelationToXML - object RelationXML does not exist");
  OUT("LagrangianR::saveRelationToXML\n");
}

LagrangianR* LagrangianR::convert(Relation *r)
{
  cout << "LagrangianR::convert (Relation *r)" << endl;
  LagrangianR* lnlr = dynamic_cast<LagrangianR*>(r);
  return lnlr;
}

