// SWIG interface for Siconos IO
%module(directors="1", allprotected="1") IO

%include start.i

%{
#include <SiconosKernel.hpp>
#include <SiconosRestart.hpp>
%}

%include handleException.i

%import Kernel.i

%include "SiconosRestart.hpp"
