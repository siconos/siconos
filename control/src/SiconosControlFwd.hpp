#ifndef SiconosControlFwd_hpp
#define SiconosControlFwd_hpp

#include "SiconosSerialization.hpp"
#include "SiconosPointers.hpp"

/* Forward declarations */
DEFINE_SPTR(ControlSimulation)
DEFINE_SPTR(ControlZOHSimulation)
DEFINE_SPTR(ControlLsodarSimulation)
DEFINE_SPTR(ControlManager)

DEFINE_SPTR(Sensor)
DEFINE_SPTR(ControlSensor)
DEFINE_SPTR(LinearSensor)

DEFINE_SPTR(Actuator)
DEFINE_SPTR(CommonSMC)
DEFINE_SPTR(ExplicitLinearSMC)
DEFINE_SPTR(LinearSMC)
DEFINE_SPTR(LinearSMCOT2)
DEFINE_SPTR(LinearSMCimproved)
DEFINE_SPTR(PID)
DEFINE_SPTR(Twisting)
DEFINE_SPTR(ExplicitTwisting)
DEFINE_SPTR(RegularTwisting)

DEFINE_SPTR(Observer)
DEFINE_SPTR(LuenbergerObserver)
DEFINE_SPTR(SlidingReducedOrderObserver)

#endif
