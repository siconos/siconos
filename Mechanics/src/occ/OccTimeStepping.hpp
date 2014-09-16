#ifndef OccTimeStepping_hpp
#define OccTimeStepping_hpp

#include "MechanicsFwd.hpp"

#include <TimeStepping.hpp>

class OccTimeStepping : public TimeStepping
{

public:

  OccTimeStepping(SP::TimeDiscretisation td) : TimeStepping(td) {};

  virtual void updateWorldFromDS();

};

#endif
