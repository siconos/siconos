#ifndef TIMESTEPPING_H
#define TIMESTEPPING_H

#include "Strategy.h"
#include "StrategyXML.h"
#include <iostream>
#include <vector>
#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"


//using namespace std;

/** \class TimeStepping
 *  \brief It's a way to drive a simulation, where the resolution is only managed by the time
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 *
 *
 * détail de fonctionnalité de la classe, caractéristiques, ...
 */
class TimeStepping : public Strategy
{
public:

  /** \fn TimeStepping()
   *  \brief Default constructor
   */
  TimeStepping();

  /** \fn TimeStepping(StrategyXML*, Model*)
   *  \brief constructor with XML object of the TimeStepping
   *  \param StrategyXML* : the XML object corresponding
   *  \param Model* : the Model which contains the Strategy
   */
  TimeStepping(StrategyXML*, Model*);

  ~TimeStepping();

  /** \fn void createStrategy(StrategyXML * strategyXML, Model * model, TimeDiscretisation * timediscretisation)
   *  \brief allows to create the Strategy with an xml file, or the needed data
   *  \param StrategyXML* : the StrategyXML linked to this Strategy
   *  \param Model* : the Model which contains this Strategy
   *  \exception RuntimeException
   */
  void createStrategy(StrategyXML * strategyXML, Model * model);//, TimeDiscretisation * timediscretisation=NULL);


  /** \fn TimeStepping* convert (Strategy* str)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Strategy* : the Strategy which must be converted
   * \return a pointer on the Strategy if it is of the right type, NULL otherwise
   */
  static TimeStepping* convert(Strategy* str);
};

#endif // TIMESTEPPING_H
