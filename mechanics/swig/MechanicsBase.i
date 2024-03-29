// -*- c++ -*-
// Base SWIG interface for siconos mechanics

%include start.i
// generated docstrings from doxygen xml output
 //%include mechanics-docstrings.i

%include sharedPointers.i

%{
#include <MechanicsFwd.hpp>
%}
%include <MechanicsFwd.hpp>

%include serialization.i

%include picklable.i

%include handleException.i
%include stl.i
%include KernelTypes.i

%{
#include <SiconosKernel.hpp>
#include <boost/typeof/typeof.hpp>
#include <vector>
#include "SiconosFwd.hpp"
%}

%import kernel.i

%include pyRegister.i

 // force the definition of SWIGTYPE_p_Interaction...
typedef Interaction Interaction;
