// -*- c++ -*-
%module(package="siconos.control", directors="1", allprotected="1") simulation

%include ControlBase.i

PY_FULL_REGISTER(ControlSimulation);
PY_FULL_REGISTER(ControlLsodarSimulation);
PY_FULL_REGISTER(ControlZOHSimulation);
PY_FULL_REGISTER(ControlManager);

