// -*- c++ -*-
%module(directors="1", allprotected="1") Base

%include start.i

#undef WITH_IO

#ifdef WITH_IO
%{
#include <SiconosFull.hpp>
%}
#endif
%include picklable.i

%include path.i

%include handleException.i

%include sharedPointers.i

%include stl.i

%include KernelTypes.i
%{
#include <SiconosKernel.hpp>
#include <boost/typeof/typeof.hpp>
%}

%import Kernel.i

%include pyRegister.i

%fragment("NumPy_Fragments");


%{
#include <MechanicsFwd.hpp>
%}
%include <MechanicsFwd.hpp>

PY_FULL_REGISTER(SpaceFilter);                                             
PY_FULL_REGISTER(SiconosBodies);
PY_FULL_REGISTER(ExternalBody);
