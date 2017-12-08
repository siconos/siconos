// -*- c++ -*-
%module(package="siconos.control", directors="1", allprotected="1") observer

%include ControlBase.i

PY_REGISTER_WITHOUT_DIRECTOR(Observer);
%include Observer.hpp
PY_FULL_REGISTER(LuenbergerObserver);
PY_FULL_REGISTER(SlidingReducedOrderObserver);

