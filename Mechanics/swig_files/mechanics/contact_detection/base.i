// -*- c++ -*-
%module(directors="1", allprotected="1") Base

%include start.i

%include sharedPointers.i

%{
#include <MechanicsFwd.hpp>
%}
%include <MechanicsFwd.hpp>

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

%include stl.i

%include KernelTypes.i
%{
#include <SiconosKernel.hpp>
#include <boost/typeof/typeof.hpp>
%}

%import Kernel/Kernel.i

%include pyRegister.i

%fragment("NumPy_Fragments");

// suppress warning
%ignore  STD11::enable_shared_from_this< Hashed >;
%template (sharedHashed) STD11::enable_shared_from_this< Hashed >;

PY_FULL_REGISTER(SpaceFilter);
PY_FULL_REGISTER(SiconosBodies);
PY_FULL_REGISTER(ExternalBody);
