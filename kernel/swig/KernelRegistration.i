// Warning: 
// You have to PY_REGISTER base classe before derivated classes
#undef PY_REGISTER
#define KERNEL_REGISTRATION()                                           \
  PY_REGISTER(SiconosMemory, Kernel)                                            \
  PY_REGISTER(NonSmoothLaw, Kernel);                                            \
  PY_REGISTER(NewtonImpactNSL, Kernel);                                         \
  PY_REGISTER(NewtonImpactFrictionNSL, Kernel);                                 \
  PY_REGISTER(NewtonImpactRollingFrictionNSL, Kernel);                  \
  PY_REGISTER(MixedComplementarityConditionNSL, Kernel);                        \
  PY_REGISTER(ComplementarityConditionNSL, Kernel);                             \
  PY_REGISTER(EqualityConditionNSL, Kernel);                                    \
  PY_REGISTER(MultipleImpactNSL, Kernel);                                       \
  PY_REGISTER(RelayNSL, Kernel);                                                \
  PY_REGISTER(NormalConeNSL, Kernel);                                           \
  PY_REGISTER(DynamicalSystem, Kernel);                                         \
  PY_REGISTER(NonSmoothDynamicalSystem, Kernel);                                \
  PY_REGISTER(Topology, Kernel);                                                \
  PY_REGISTER(SecondOrderDS, Kernel);                                   \
  PY_REGISTER(LagrangianDS, Kernel);                                            \
  PY_REGISTER(LagrangianLinearTIDS, Kernel);                                    \
  PY_REGISTER(LagrangianLinearDiagonalDS, Kernel);                              \
  PY_REGISTER(NewtonEulerDS, Kernel);                                           \
  PY_REGISTER(FirstOrderNonLinearDS, Kernel);                                   \
  PY_REGISTER(FirstOrderLinearDS, Kernel);                                      \
  PY_REGISTER(FirstOrderLinearTIDS, Kernel);                                    \
  PY_REGISTER(Relation, Kernel);                                                \
  PY_REGISTER(LagrangianR, Kernel);                                             \
  PY_REGISTER(LagrangianLinearTIR, Kernel);                                     \
  PY_REGISTER(LagrangianRheonomousR, Kernel);                                   \
  PY_REGISTER(LagrangianScleronomousR, Kernel);                                 \
  PY_REGISTER(LagrangianCompliantR, Kernel);                                    \
  PY_REGISTER(NewtonEulerR, Kernel);                                            \
  PY_REGISTER(NewtonEuler1DR, Kernel);                                  \
  PY_REGISTER(NewtonEuler3DR, Kernel);                                  \
  PY_REGISTER(NewtonEuler5DR, Kernel);                                  \
  PY_REGISTER(FirstOrderR, Kernel);                                             \
  PY_REGISTER(FirstOrderNonLinearR, Kernel);                                    \
  PY_REGISTER(FirstOrderType1R, Kernel);                                        \
  PY_REGISTER(FirstOrderType2R, Kernel);                                        \
  PY_REGISTER(FirstOrderLinearR, Kernel);                                       \
  PY_REGISTER(FirstOrderLinearTIR, Kernel);                                     \
  PY_REGISTER(Interaction, Kernel);                                             \
  PY_REGISTER(TimeDiscretisation, Kernel);                                      \
  PY_REGISTER(OneStepNSProblem, Kernel);                                        \
  PY_REGISTER(OneStepIntegrator, Kernel);                                       \
  PY_REGISTER(LinearOSNS, Kernel);                                              \
  PY_REGISTER(LsodarOSI, Kernel);                                               \
  PY_REGISTER(LCP, Kernel);                                                     \
  PY_REGISTER(AVI, Kernel);                                                     \
  PY_REGISTER(QP, Kernel);                                                      \
  PY_REGISTER(Relay, Kernel);                                                   \
  PY_REGISTER(MLCP, Kernel);                                                    \
  PY_REGISTER(MLCPProjectOnConstraints, Kernel);                                \
  PY_REGISTER(GenericMechanical, Kernel);                                       \
  PY_REGISTER(FrictionContact, Kernel);                                         \
  PY_REGISTER(GlobalFrictionContact, Kernel);                           \
  PY_REGISTER(RollingFrictionContact, Kernel);                          \
  PY_REGISTER(EulerMoreauOSI, Kernel);                                          \
  PY_REGISTER(MoreauJeanOSI, Kernel);                                           \
  PY_REGISTER(MoreauJeanBilbaoOSI, Kernel);                                     \
  PY_REGISTER(MoreauJeanCombinedProjectionOSI, Kernel);                         \
  PY_REGISTER(MoreauJeanDirectProjectionOSI, Kernel);                           \
  PY_REGISTER(MoreauJeanGOSI, Kernel);                                          \
  PY_REGISTER(ZeroOrderHoldOSI, Kernel);                                        \
  PY_REGISTER(Simulation, Kernel);                                              \
  PY_REGISTER(TimeStepping, Kernel);                                            \
  PY_REGISTER(TimeSteppingCombinedProjection, Kernel);                          \
  PY_REGISTER(TimeSteppingDirectProjection, Kernel);                            \
  PY_REGISTER(InteractionManager, Kernel);                                      \
  PY_REGISTER(EventDriven, Kernel);                                             \
  PY_REGISTER(EventsManager, Kernel);                                           \
  PY_REGISTER(Event, Kernel);                                                   \
  PY_REGISTER(BoundaryCondition, Kernel);                                       \
  PY_REGISTER(HarmonicBC, Kernel);                                              \
  PY_REGISTER(FixedBC, Kernel);                                                 \
  PY_REGISTER(OSNSMatrix, Kernel);                                              \
  PY_REGISTER(BlockCSRMatrix, Kernel);

%feature("nodirector") SecondOrderDS::dimension;
%feature("nodirector") SecondOrderDS::q;
%feature("nodirector") SecondOrderDS::q0;
%feature("nodirector") SecondOrderDS::velocity;
%feature("nodirector") SecondOrderDS::velocity0;
%feature("nodirector") SecondOrderDS::acceleration;
%feature("nodirector") SecondOrderDS::forces;
%feature("nodirector") SecondOrderDS::jacobianqForces;
%feature("nodirector") SecondOrderDS::jacobianvForces;
%feature("nodirector") SecondOrderDS::p;
%feature("nodirector") SecondOrderDS::mass;
%feature("nodirector") SecondOrderDS::inverseMass;
%feature("nodirector") SecondOrderDS::qMemory;
%feature("nodirector") SecondOrderDS::velocityMemory;
%feature("nodirector") SecondOrderDS::forcesMemory;
%feature("nodirector") SecondOrderDS::reactionToBoundaryConditions;


