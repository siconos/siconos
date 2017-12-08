// -*- c++ -*-
// SWIG interface for Siconos Mechanics/bodies
%module(package="siconos.mechanics.collision", directors="1", allprotected="1") bodies

%include MechanicsBase.i

PY_FULL_REGISTER(CircularDS);
PY_FULL_REGISTER(CircularR);
PY_FULL_REGISTER(Disk);
PY_FULL_REGISTER(Circle);
PY_FULL_REGISTER(DiskDiskR);
PY_FULL_REGISTER(DiskPlanR);
PY_FULL_REGISTER(DiskMovingPlanR);
PY_FULL_REGISTER(SphereLDS);
PY_FULL_REGISTER(SphereNEDS);
PY_FULL_REGISTER(SphereLDSPlanR);
PY_FULL_REGISTER(SphereNEDSPlanR);
PY_FULL_REGISTER(SphereLDSSphereLDSR);
PY_FULL_REGISTER(SphereNEDSSphereNEDSR);
