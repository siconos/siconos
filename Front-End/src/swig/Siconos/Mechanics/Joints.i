
// SWIG interface for Siconos Mechanics/joints
%module(directors="1", allprotected="1") Joints

%include start.i

%include path.i

%{
#include <SiconosKernel.hpp>
%}

// common declarations

%include handleException.i
%include sharedPointers.i
%include KernelTypes.i

%{
#include <SiconosKernel.hpp>
%}
%import Kernel.i

%include pyRegister.i

PY_FULL_REGISTER(KneeJointR);
PY_FULL_REGISTER(PivotJointR);
PY_FULL_REGISTER(PrismaticJointR);
