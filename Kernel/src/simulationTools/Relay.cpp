
#include "Relay.h"

Relay::Relay(): OneStepNSProblem()
{
  this->nspbType = RELAY_OSNSP;
}

Relay::Relay(OneStepNSProblemXML* osnspbxml): OneStepNSProblem(osnspbxml)
{
  this->nspbType = RELAY_OSNSP;
}

Relay::~Relay()
{}

void Relay::fillNSProblemWithNSProblemXML()
{
}

void Relay::saveNSProblemToXML()
{
  OneStepNSProblem::saveNSProblemToXML();
}

void Relay::createOneStepNSProblem(OneStepNSProblemXML * osnspbXML, Strategy * strategy)
{
  if (osnspbXML != NULL)
  {
    this->onestepnspbxml = osnspbXML;
    this->nspbType = RELAY_OSNSP;

    this->fillNSProblemWithNSProblemXML();
  }
  else
  {
    this->strategy = strategy;
    this->nspbType = RELAY_OSNSP;
    this->fillInteractionVector();
  }
  this->init();
}

Relay* Relay::convert(OneStepNSProblem* osnsp)
{
  cout << "Relay::convert (DynamicalSystem* osnsp)" << endl;
  Relay* r = dynamic_cast<Relay*>(osnsp);
  return r;
}

