// -*- c++ -*-
%module(directors="1", allprotected="1") Sensor

%include ControlBase.i

PY_REGISTER_WITHOUT_DIRECTOR(Sensor);
%include Sensor.hpp
PY_FULL_REGISTER(ControlSensor);
PY_FULL_REGISTER(LinearSensor);

