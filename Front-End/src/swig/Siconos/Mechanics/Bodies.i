
// SWIG interface for Siconos Mechanics/bodies
%module(directors="1", allprotected="1") Bodies

%include start.i

%include path.i

%include handleException.i

%include sharedPointers.i

%include KernelTypes.i

%{
#include <SiconosKernel.hpp>
%}
%import Kernel.i

%include pyRegister.i

PY_REGISTER(Disk);
PY_REGISTER(Circle);
PY_REGISTER(DiskDiskR);
PY_REGISTER(DiskPlanR);
PY_REGISTER(DiskMovingPlanR);
PY_REGISTER(SphereLDS);
PY_REGISTER(SphereNEDS);
PY_REGISTER(SphereLDSPlanR);
PY_REGISTER(SphereNEDSPlanR);
PY_REGISTER(SphereLDSSphereLDSR);
PY_REGISTER(SphereNEDSSphereNEDSR);


