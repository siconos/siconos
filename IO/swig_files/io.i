// -*- c++ -*-
// SWIG interface for Siconos IO
%module(package="siconos", directors="1", allprotected="1") io

%include start.i

%{
#include <SiconosKernel.hpp>
#include <SiconosRestart.hpp>
%}

%include handleException.i

%include sharedPointers.i

%include KernelTypes.i

%import kernel.i

%include "SiconosRestart.hpp"

#ifdef HAVE_SICONOS_MECHANICS
%include <MechanicsIO.hpp>
%{
#include <MechanicsIO.hpp>
%}
#endif
