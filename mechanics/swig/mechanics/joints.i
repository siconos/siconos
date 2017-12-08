// -*- c++ -*-
// SWIG interface for Siconos Mechanics/joints
%module(package="siconos.mechanics", directors="1", allprotected="1") joints

// Ignore some shadowed (redundant for Python) functions
%ignore JointFrictionR(SP::NewtonEulerJointR, unsigned int);

%include MechanicsBase.i

PY_FULL_REGISTER(NewtonEulerJointR); // Abstract
PY_FULL_REGISTER(KneeJointR);
PY_FULL_REGISTER(PivotJointR);
PY_FULL_REGISTER(PrismaticJointR);
PY_FULL_REGISTER(FixedJointR);
PY_FULL_REGISTER(CylindricalJointR);
PY_FULL_REGISTER(CouplerJointR);
PY_FULL_REGISTER(JointStopR);
PY_FULL_REGISTER(JointFrictionR);

%inline
%{
  // For converting interaction.relations() to known Relations
  SP::NewtonEulerJointR cast_NewtonEulerJointR(SP::Relation rel)
    { return std11::dynamic_pointer_cast<NewtonEulerJointR>(rel); }
  SP::KneeJointR cast_KneeJointR(SP::Relation rel)
    { return std11::dynamic_pointer_cast<KneeJointR>(rel); }
  SP::PivotJointR cast_PivotJointR(SP::Relation rel)
    { return std11::dynamic_pointer_cast<PivotJointR>(rel); }
  SP::PrismaticJointR cast_PrismaticJointR(SP::Relation rel)
    { return std11::dynamic_pointer_cast<PrismaticJointR>(rel); }
  SP::FixedJointR cast_FixedJointR(SP::Relation rel)
    { return std11::dynamic_pointer_cast<FixedJointR>(rel); }
  SP::CylindricalJointR cast_CylindricalJointR(SP::Relation rel)
    { return std11::dynamic_pointer_cast<CylindricalJointR>(rel); }
  SP::CouplerJointR cast_CouplerJointR(SP::Relation rel)
    { return std11::dynamic_pointer_cast<CouplerJointR>(rel); }
  SP::JointStopR cast_JointStopR(SP::Relation rel)
    { return std11::dynamic_pointer_cast<JointStopR>(rel); }
  SP::JointFrictionR cast_JointFrictionR(SP::Relation rel)
    { return std11::dynamic_pointer_cast<JointFrictionR>(rel); }
%}
