// -*- c++ -*-
%module(package="siconos.control", directors="1", allprotected="1") controller

%include ControlBase.i

PY_REGISTER_WITHOUT_DIRECTOR(Actuator, Control);
%include Actuator.hpp
PY_FULL_REGISTER(PID, Control);
PY_FULL_REGISTER(CommonSMC, Control);
PY_FULL_REGISTER(LinearSMC, Control);
PY_FULL_REGISTER(ExplicitLinearSMC, Control);
PY_FULL_REGISTER(LinearSMCOT2, Control);
PY_FULL_REGISTER(LinearSMCimproved, Control);
PY_FULL_REGISTER(Twisting, Control);
PY_FULL_REGISTER(RegularTwisting, Control);
PY_FULL_REGISTER(ExplicitTwisting, Control);
