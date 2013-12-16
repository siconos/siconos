// -*- c++ -*-
// SWIG interface for Siconos Mechanics/joints
%module(directors="1", allprotected="1") Joints

%include start.i

#undef WITH_IO

#ifdef WITH_IO
%{
#include <SiconosFull.hpp>
%}
#endif
%include picklable.i

%include path.i

%{
#include <SiconosKernel.hpp>
%}

// common declarations

%include handleException.i
%include sharedPointers.i
%include KernelTypes.i

%import Kernel.i

%include pyRegister.i

%{
#include <MechanicsFwd.hpp>
%}
%include <MechanicsFwd.hpp>


PY_FULL_REGISTER(KneeJointR);
PY_FULL_REGISTER(PivotJointR);
PY_FULL_REGISTER(PrismaticJointR);
