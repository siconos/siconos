// -*- c++ -*-
// Base SWIG interface for Siconos Control

%include start.i
// generated docstrings from doxygen xml output
%include control-docstrings.i

#undef WITH_IO
#undef WITH_SERIALIZATION

#ifdef WITH_SERIALIZATION
%{
#include <SiconosFull.hpp>
%}
#endif

%include picklable.i

%include handleException.i

%include sharedPointers.i

%include stl.i

%include KernelTypes.i

%{
  #include <SiconosKernel.hpp>
%}
%import kernel.i

%{
#include <SiconosControlFwd.hpp>
%}
%include <SiconosControlFwd.hpp>

// force the definition of SWIGTYPE_p_Interaction...
typedef Interaction Interaction;

%include ControlTypemaps.i

%fragment("NumPy_Fragments");

// this does not work, why why why ???
//%include ../KernelRegistration.i
//%include pySharedPtr.i
//KERNEL_REGISTRATION();

%include pyRegister.i

