
#include "EventDriven.h"

EventDriven::EventDriven()
{
  this->strategyType = EVENTDRIVEN_STRATEGY;
}

EventDriven::EventDriven(StrategyXML* strxml, Model *model): Strategy(strxml,
      model)
{
  this->strategyType = EVENTDRIVEN_STRATEGY;
}

EventDriven::~EventDriven()
{}

void EventDriven::createStrategy(StrategyXML * strategyXML, Model * model)//, TimeDiscretisation * timediscretisation)
{
  if (strategyXML != NULL)
  {
    //this->timeDiscretisation = NULL;
    this->timeDiscretisation = timeDiscretisation;
    if (timeDiscretisation != NULL) this->timeDiscretisation->setStrategy(this);

    this->integratorVector.clear();
    this->nsProblem = NULL;

    this->strategyxml = strategyXML;
    //  this->nsds = nsds;
    this->model = model;

    this->fillStrategyWithStrategyXML();
    this->linkStrategyXML();
  }
  else
  {
    this->strategyxml = NULL;
    this->strategyType = EVENTDRIVEN_STRATEGY;
    this->model = model;
  }
}


EventDriven* EventDriven::convert(Strategy *str)
{
  cout << "EventDriven::convert (Strategy *str)" << endl;
  EventDriven* ed = dynamic_cast<EventDriven*>(str);
  return ed;
}

