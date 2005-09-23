#include "LagrangianRXML.h"
using namespace std;

LagrangianRXML::LagrangianRXML()
  : RelationXML()
{}

LagrangianRXML::LagrangianRXML(xmlNode * LNLRelationNode)
  : RelationXML(LNLRelationNode)
{}

LagrangianRXML::~LagrangianRXML()
{}

string  LagrangianRXML::getComputeInputPlugin() const
{
  return  SiconosDOMTreeTools::getStringAttributeValue(computeInputNode, COMPUTE_INPUT_TAG);
}

string  LagrangianRXML::getComputeOutputPlugin() const
{
  return  SiconosDOMTreeTools::getStringAttributeValue(computeOutputNode, COMPUTE_OUTPUT_TAG);
}
