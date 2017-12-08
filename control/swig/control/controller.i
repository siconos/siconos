// -*- c++ -*-
%module(package="siconos.control", directors="1", allprotected="1") controller

%include ControlBase.i

PY_REGISTER_WITHOUT_DIRECTOR(Actuator);
%include Actuator.hpp
PY_FULL_REGISTER(PID);
PY_FULL_REGISTER(CommonSMC);
PY_FULL_REGISTER(LinearSMC);
PY_FULL_REGISTER(ExplicitLinearSMC);
PY_FULL_REGISTER(LinearSMCOT2);
PY_FULL_REGISTER(LinearSMCimproved);

