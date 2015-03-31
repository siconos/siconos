/* mbtb.i this file contains exported API of the MBTB library.*/
%module mbtb

%include start.i
#undef WITH_IO
#undef WITH_SERIALIZATION

%include path.i

%{
#include <SiconosKernel.hpp>
#include <MBTB_PYTHON_API.hpp>
%}

%include handleException.i
%include sharedPointers.i
%include KernelTypes.i

%import Kernel/Kernel.i

%include <MBTB_PYTHON_API.hpp>
 //%include "MBTB_Body.hpp"
