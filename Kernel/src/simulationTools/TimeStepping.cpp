
#include "TimeStepping.h"

#include "check.h"

TimeStepping::TimeStepping()
{
  IN("TimeStepping::TimeStepping()\n");
  this->strategyType = TIMESTEPPING_STRATEGY;
  OUT("TimeStepping::TimeStepping()\n");
}

TimeStepping::TimeStepping(StrategyXML* strxml, Model *model): Strategy(strxml,
      model)
{
  IN("TimeStepping::TimeStepping(StrategyXML* strxml, Model *model)\n");
  this->strategyType = TIMESTEPPING_STRATEGY;
  OUT("TimeStepping::TimeStepping(StrategyXML* strxml, Model *model)\n");
}

TimeStepping::~TimeStepping()
{}


void TimeStepping::createStrategy(StrategyXML * strategyXML, Model * model)//,  TimeDiscretisation * timediscretisation)
{
  IN("TimeStepping::createStrategy(StrategyXML * strategyXML, Model * model)\n");
  if (strategyXML != NULL)
  {
    // this->timeDiscretisation = NULL;
    // \warning where does ¨timeDiscretisation comes from ????  => comment the 2 next lines
    // this->timeDiscretisation = timeDiscretisation;
    // if( timeDiscretisation != NULL ) this->timeDiscretisation->setStrategy( this );

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
    this->strategyType = TIMESTEPPING_STRATEGY;
    this->strategyxml = NULL;
    this->model = model;
  }
  OUT("TimeStepping::createStrategy(StrategyXML * strategyXML, Model * model)\n");
}


TimeStepping* TimeStepping::convert(Strategy *str)
{
  cout << "TimeStepping::convert (Strategy *str)" << endl;
  TimeStepping* ts = dynamic_cast<TimeStepping*>(str);
  return ts;
}

