// -*- c++ -*-
// SWIG interface for Siconos Mechanics/joints
%module(package="mechanics", directors="1", allprotected="1") joints

%include MechanicsBase.i

PY_FULL_REGISTER(KneeJointR);
PY_FULL_REGISTER(PivotJointR);
PY_FULL_REGISTER(PrismaticJointR);
PY_FULL_REGISTER(FixedJointR);
