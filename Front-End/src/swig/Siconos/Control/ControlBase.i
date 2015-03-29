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

%include path.i

%include picklable.i

%include handleException.i

%include sharedPointers.i

%include stl.i

%include KernelTypes.i

%{
  #include <SiconosKernel.hpp>
%}
%import Kernel/Kernel.i

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

