// -*- c++ -*-
// SWIG interface for Siconos Mechanics/bodies
%module(package="siconos.mechanics.collision", directors="1", allprotected="1") bodies

%include MechanicsBase.i

PY_FULL_REGISTER(CircularDS, Mechanics);
PY_FULL_REGISTER(CircularR, Mechanics);
PY_FULL_REGISTER(Disk, Mechanics);
PY_FULL_REGISTER(Circle, Mechanics);
PY_FULL_REGISTER(DiskDiskR, Mechanics);
PY_FULL_REGISTER(DiskPlanR, Mechanics);
PY_FULL_REGISTER(DiskMovingPlanR, Mechanics);
PY_FULL_REGISTER(SphereLDS, Mechanics);
PY_FULL_REGISTER(SphereNEDS, Mechanics);
PY_FULL_REGISTER(SphereLDSPlanR, Mechanics);
PY_FULL_REGISTER(SphereNEDSPlanR, Mechanics);
PY_FULL_REGISTER(SphereLDSSphereLDSR, Mechanics);
PY_FULL_REGISTER(SphereNEDSSphereNEDSR, Mechanics);
