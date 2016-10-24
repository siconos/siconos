// -*- c++ -*-
// SWIG interface for Siconos Mechanics/collision
%module(package="collision", directors="1", allprotected="1") base

%include MechanicsBase.i

PY_REGISTER_WITHOUT_HEADER(SiconosSphere);
PY_REGISTER_WITHOUT_HEADER(SiconosPlane);
PY_REGISTER_WITHOUT_HEADER(SiconosBox);
PY_REGISTER_WITHOUT_HEADER(SiconosCylinder);
PY_REGISTER_WITHOUT_HEADER(SiconosConvexHull);
PY_FULL_REGISTER(SiconosShape);
PY_REGISTER_WITHOUT_HEADER(SiconosContactorSet);
PY_FULL_REGISTER(SiconosContactor);
PY_FULL_REGISTER(BodyDS);
PY_FULL_REGISTER(SiconosCollisionManager);
