
#include "EventDriven.h"
using namespace std;

// --- Default constructor ---
EventDriven::EventDriven(): Strategy()
{
  strategyType = EVENTDRIVEN_STRATEGY;
}

// --- XML constructor ---
EventDriven::EventDriven(StrategyXML* strxml, Model *newModel): Strategy(strxml, newModel)
{
  strategyType = EVENTDRIVEN_STRATEGY;
}

// --- Destructor ---
EventDriven::~EventDriven()
{}

void EventDriven::createStrategy(StrategyXML * newStrategyXML, Model * newModel)
{
  strategyxml = NULL;
  strategyType = EVENTDRIVEN_STRATEGY;
  model = newModel;
}


EventDriven* EventDriven::convert(Strategy *str)
{
  cout << "EventDriven::convert (Strategy *str)" << endl;
  EventDriven* ed = dynamic_cast<EventDriven*>(str);
  return ed;
}

