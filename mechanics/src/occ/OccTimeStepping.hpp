#ifndef OccTimeStepping_hpp
#define OccTimeStepping_hpp

#include "MechanicsFwd.hpp"

#include <TimeStepping.hpp>

class OccTimeStepping : public TimeStepping
{

public:

  OccTimeStepping(SP::NonSmoothDynamicalSystem nsds, SP::TimeDiscretisation td) : TimeStepping(nsds,td) {};

  virtual void updateWorldFromDS();

};

#endif
