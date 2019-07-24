/* mbtb.i this file contains exported API of the MBTB library.*/
%module(package="siconos.mechanics.mechanisms") mbtb

%include start.i
#undef WITH_IO
#undef WITH_SERIALIZATION

%{
#include <SiconosKernel.hpp>
#include <MBTB_DATA.hpp>
#include <MBTB_internalTool.hpp>
#include <MBTB_PYTHON_API.hpp>
#include <ace.h>

  unsigned int sUseGravity = 0;
%}

%include handleException.i
%include sharedPointers.i
%include KernelTypes.i

%import kernel.i

%include <MBTB_DATA.hpp>
%include <MBTB_internalTool.hpp>
%include <MBTB_PYTHON_API.hpp>
%include <ace.h>

 //%include "MBTB_Body.hpp"
