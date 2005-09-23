#ifndef EVENTDRIVEN_H
#define EVENTDRIVEN_H

#include "Strategy.h"

/** \class EventDriven
 *  \brief It's a way to manage a simulation, where the event a more important than the time
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 *
 */

class EventDriven : public Strategy
{
public:

  /** \fn EventDriven()
   *  \brief defaut constructor
   */
  EventDriven();

  /** \fn EventDriven(StrategyXML*, Model*)
   *  \brief constructor with XML object of the EventDriven
   *  \param StrategyXML* : the XML object corresponding
   *  \param Model* : the Model which contains the Strategy
   */
  EventDriven(StrategyXML*, Model*);

  ~EventDriven();

  /** \fn void createStrategy(StrategyXML * strategyXML, Model * model, TimeDiscretisation * timediscretisation)
   *  \brief create the Strategy with an xml file, or the needed data
   *  \param StrategyXML* : the StrategyXML linked to this Strategy
   *  \param Model* : the Model which contains this Strategy
   *  \exception RuntimeException
   */
  void createStrategy(StrategyXML * strategyXML, Model * model);

  /** \fn EventDriven* convert (Strategy* str)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Strategy* : the Strategy which must be converted
   * \return a pointer on the Strategy if it is of the right type, NULL otherwise
   */
  static EventDriven* convert(Strategy* str);

};

#endif // EVENTDRIVEN_H
