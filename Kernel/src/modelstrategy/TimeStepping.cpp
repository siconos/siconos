
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

//$Log: TimeStepping.cpp,v $
//Revision 1.20  2005/03/08 14:23:44  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.19  2005/02/01 11:08:42  charlety
//
//_ some displays of values during computations suppressed.
//
//Revision 1.18  2005/01/31 16:26:27  charlety
//
//_ Added a method named "convert" to classes which inherits from another. This is necessary for Python interface, in order to be able to use down-casting mechanism.
//
//Revision 1.17  2004/12/06 10:10:34  jbarbier
//- integration of Numerics and use of Numerics on the bouncing ball sample
//
//- Numerics is now needed to run the bouncing ball sample!
//
//Revision 1.16  2004/09/16 11:35:25  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.15  2004/09/07 07:31:12  jbarbier
//- create strategies methods of the Model done
//
//Revision 1.14  2004/08/18 14:37:19  jbarbier
//- creation of Model, NSDS, Strategy(TimeStepping and EventDriven) and
//DynamicalSystem available when the creation is in a command program
//
//Revision 1.13  2004/08/12 11:55:19  jbarbier
//- new methods createModel, createNSDS, createStrategy, ...
//they now allow to make the link with upper objects of the platform
//it will be used for the creation of the platform without XML input file
//
//- the createModel method is finished but the attributes of the other objects
//of the platform are missing for the conctruction
//
//Revision 1.12  2004/08/03 12:07:12  jbarbier
//- all test on th eModel are successfull
//
//- new tests on the Model with the opening of XML file
//
//- link TimeDiscretisation -> Strategy
//
//- attribute T of the Model is now optional
//
//Revision 1.11  2004/07/29 14:25:41  jbarbier
