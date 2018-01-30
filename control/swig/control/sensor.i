// -*- c++ -*-
%module(package="siconos.control", directors="1", allprotected="1") sensor

%include ControlBase.i

PY_REGISTER_WITHOUT_DIRECTOR(Sensor, Control);
%include Sensor.hpp
PY_FULL_REGISTER(ControlSensor, Control);
PY_FULL_REGISTER(LinearSensor, Control);

