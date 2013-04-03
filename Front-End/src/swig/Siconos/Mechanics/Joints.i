
// SWIG interface for Siconos Mechanics/joints
%module(directors="1", allprotected="1") Joints

%include start.i

%include path.i

%{
#include <SiconosKernel.hpp>
%}

// common declarations

%include handleException.i

%import Kernel.i

%include pyRegister.i

PY_REGISTER(KneeJointR);
PY_REGISTER(PivotJointR);
PY_REGISTER(PrismaticJointR);
