// -*- c++ -*-
// Base SWIG interface for Siconos Control

%include start.i

#undef WITH_IO
#undef WITH_SERIALIZATION

#ifdef WITH_IO
%{
#include <SiconosFull.hpp>
%}
#endif

%include picklable.i

%include path.i

%include handleException.i

%include sharedPointers.i

%include KernelTypes.i

%{
#include <SiconosKernel.hpp>
%}
%import Kernel.i

%include pyRegister.i

%{
#include <SiconosControlFwd.hpp>
%}
%include <SiconosControlFwd.hpp>

// force the definition of SWIGTYPE_p_Interaction...
typedef Interaction Interaction;

%include ControlTypemaps.i
