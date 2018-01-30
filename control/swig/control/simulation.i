// -*- c++ -*-
%module(package="siconos.control", directors="1", allprotected="1") simulation

%include ControlBase.i

PY_FULL_REGISTER(ControlSimulation, Control);
PY_FULL_REGISTER(ControlLsodarSimulation, Control);
PY_FULL_REGISTER(ControlZOHSimulation, Control);
PY_FULL_REGISTER(ControlManager, Control);

