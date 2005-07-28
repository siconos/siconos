
#include "TimeStepping.h"
using namespace std;

// --- Default constructor ---
TimeStepping::TimeStepping(): Strategy()
{
  strategyType = TIMESTEPPING_STRATEGY;
}

// --- XML constructor ---
TimeStepping::TimeStepping(StrategyXML* strxml, Model *newModel): Strategy(strxml, newModel)
{
  strategyType = TIMESTEPPING_STRATEGY;
}

// --- Destructor ---
TimeStepping::~TimeStepping()
{}

void TimeStepping::createStrategy(StrategyXML * newStrategyXML, Model * newModel)
{
  IN("TimeStepping::createStrategy(StrategyXML * strategyXML, Model * model)\n");
  strategyType = TIMESTEPPING_STRATEGY;
  strategyxml = NULL;
  model = newModel;
  OUT("TimeStepping::createStrategy(StrategyXML * strategyXML, Model * model)\n");
}

TimeStepping* TimeStepping::convert(Strategy *str)
{
  cout << "TimeStepping::convert (Strategy *str)" << endl;
  TimeStepping* ts = dynamic_cast<TimeStepping*>(str);
  return ts;
}

