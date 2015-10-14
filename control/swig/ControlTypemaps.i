

%define DEFINE_TYPEMAPS(TYPE)
%shared_ptr(TYPE)
%enddef

// generate this list with grep REGISTER *.i | sed -n 's/.*(\(.*\)).*/DEFINE_TYPEMAPS(\1);/p'
// and remove the last item

DEFINE_TYPEMAPS(Actuator);
DEFINE_TYPEMAPS(PID);
DEFINE_TYPEMAPS(CommonSMC);
DEFINE_TYPEMAPS(LinearSMC);
DEFINE_TYPEMAPS(ExplicitLinearSMC);
DEFINE_TYPEMAPS(LinearSMCOT2);
DEFINE_TYPEMAPS(LinearSMCimproved);
DEFINE_TYPEMAPS(Observer);
DEFINE_TYPEMAPS(LuenbergerObserver);
DEFINE_TYPEMAPS(SlidingReducedOrderObserver);
DEFINE_TYPEMAPS(Sensor);
DEFINE_TYPEMAPS(ControlSensor);
DEFINE_TYPEMAPS(LinearSensor);
DEFINE_TYPEMAPS(ControlSimulation);
DEFINE_TYPEMAPS(ControlLsodarSimulation);
DEFINE_TYPEMAPS(ControlZOHSimulation);
DEFINE_TYPEMAPS(ControlManager);
