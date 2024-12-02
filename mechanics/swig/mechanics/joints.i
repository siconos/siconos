// -*- c++ -*-
// SWIG interface for Siconos Mechanics/joints
%module(package="siconos.mechanics", directors="1", allprotected="1") joints

// Ignore some shadowed (redundant for Python) functions
%ignore JointFrictionR(std::shared_ptr<siconos::joints::NewtonEulerJointR>, unsigned int);

%include MechanicsBase.i

PY_FULL_REGISTER(NewtonEulerJointR, Mechanics); // Abstract
PY_FULL_REGISTER(KneeJointR, Mechanics);
PY_FULL_REGISTER(PivotJointR, Mechanics);
PY_FULL_REGISTER(PrismaticJointR, Mechanics);
PY_FULL_REGISTER(FixedJointR, Mechanics);
PY_FULL_REGISTER(CylindricalJointR, Mechanics);
PY_FULL_REGISTER(CouplerJointR, Mechanics);
PY_FULL_REGISTER(JointStopR, Mechanics);
PY_FULL_REGISTER(JointFrictionR, Mechanics);

%inline
%{
  // For converting interaction.relations() to known Relations
  std::shared_ptr<siconos::joints::NewtonEulerJointR> cast_NewtonEulerJointR(std::shared_ptr<siconos::modeling::Relation> rel)
    { return std::dynamic_pointer_cast<NewtonEulerJointR>(rel); }
  std::shared_ptr<siconos::joints::KneeJointR> cast_KneeJointR(std::shared_ptr<siconos::modeling::Relation> rel)
    { return std::dynamic_pointer_cast<KneeJointR>(rel); }
  std::shared_ptr<siconos::joints::PivotJointR> cast_PivotJointR(std::shared_ptr<siconos::modeling::Relation> rel)
    { return std::dynamic_pointer_cast<PivotJointR>(rel); }
  std::shared_ptr<siconos::joints::PrismaticJointR> cast_PrismaticJointR(std::shared_ptr<siconos::modeling::Relation> rel)
    { return std::dynamic_pointer_cast<PrismaticJointR>(rel); }
  std::shared_ptr<siconos::joints::FixedJointR> cast_FixedJointR(std::shared_ptr<siconos::modeling::Relation> rel)
    { return std::dynamic_pointer_cast<FixedJointR>(rel); }
  std::shared_ptr<siconos::joints::CylindricalJointR> cast_CylindricalJointR(std::shared_ptr<siconos::modeling::Relation> rel)
    { return std::dynamic_pointer_cast<CylindricalJointR>(rel); }
  std::shared_ptr<siconos::joints::CouplerJointR> cast_CouplerJointR(std::shared_ptr<siconos::modeling::Relation> rel)
    { return std::dynamic_pointer_cast<CouplerJointR>(rel); }
  std::shared_ptr<siconos::joints::JointStopR> cast_JointStopR(std::shared_ptr<siconos::modeling::Relation> rel)
    { return std::dynamic_pointer_cast<JointStopR>(rel); }
  std::shared_ptr<siconos::joints::JointFrictionR> cast_JointFrictionR(std::shared_ptr<siconos::modeling::Relation> rel)
    { return std::dynamic_pointer_cast<JointFrictionR>(rel); }
%}
