#include "LagrangianNonLinearRXML.h"
using namespace std;

LagrangianNonLinearRXML::LagrangianNonLinearRXML()
  : RelationXML()
{}

LagrangianNonLinearRXML::LagrangianNonLinearRXML(xmlNode * LNLRelationNode)
  : RelationXML(LNLRelationNode)
{}

LagrangianNonLinearRXML::~LagrangianNonLinearRXML()
{}

string  LagrangianNonLinearRXML::getComputeInputPlugin() const
{
  return  SiconosDOMTreeTools::getStringAttributeValue(computeInputNode, COMPUTE_INPUT_TAG);
}

string  LagrangianNonLinearRXML::getComputeOutputPlugin() const
{
  return  SiconosDOMTreeTools::getStringAttributeValue(computeOutputNode, COMPUTE_OUTPUT_TAG);
}
