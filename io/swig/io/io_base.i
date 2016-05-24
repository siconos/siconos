// -*- c++ -*-
// Base SWIG interface for Siconos IO
%module(package="siconos", directors="1", allprotected="1") io_base

%include start.i
// generated docstrings from doxygen xml output
%include io-docstrings.i
%include path.i
%{
#include <SiconosKernel.hpp>
#include <SiconosRestart.hpp>
%}

%include handleException.i

%include sharedPointers.i

%include KernelTypes.i

%import kernel.i

%include "SiconosRestart.hpp"

#ifdef WITH_MECHANICS
%include <MechanicsIO.hpp>
%{
#include <MechanicsIO.hpp>
%}
#endif
