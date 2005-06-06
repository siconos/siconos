
#include "Relay.h"
using namespace std;

Relay::Relay(): OneStepNSProblem()
{
  nspbType = RELAY_OSNSP;
}

Relay::Relay(OneStepNSProblemXML* osnspbxml, Strategy * newStrat):
  OneStepNSProblem(osnspbxml, newStrat)
{
  nspbType = RELAY_OSNSP;

}

Relay::~Relay()
{}

void Relay::saveNSProblemToXML()
{}

Relay* Relay::convert(OneStepNSProblem* osnsp)
{
  cout << "Relay::convert (DynamicalSystem* osnsp)" << endl;
  Relay* r = dynamic_cast<Relay*>(osnsp);
  return r;
}

