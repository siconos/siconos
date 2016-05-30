// -*- c++ -*-
// SWIG interface for Siconos Mechanics/proposed
%module(package="mechanics", directors="1", allprotected="1") proposed

%include MechanicsBase.i

PY_REGISTER_WITHOUT_HEADER(SiconosSphere);
PY_REGISTER_WITHOUT_HEADER(SiconosPlane);
PY_REGISTER_WITHOUT_HEADER(SiconosBox);
PY_REGISTER_WITHOUT_HEADER(SiconosConvexHull);
PY_FULL_REGISTER(SiconosShape);
PY_FULL_REGISTER(SiconosContactor);
PY_FULL_REGISTER(BodyDS);
PY_FULL_REGISTER(BodyTimeStepping);
PY_FULL_REGISTER(SiconosBroadphase);
PY_REGISTER_WITHOUT_HEADER(BulletOptions);
PY_REGISTER_WITHOUT_HEADER(BulletStatistics);
PY_FULL_REGISTER(BulletBroadphase);
