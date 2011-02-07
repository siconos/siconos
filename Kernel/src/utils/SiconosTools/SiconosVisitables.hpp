#ifndef SiconosVisitables_hpp
#define SiconosVisitables_hpp

#undef REGISTER
#define SICONOS_VISITABLES()                                 \
  REGISTER(DynamicalSystem)                                  \
  REGISTER(DynamicalSystemXML)                               \
  REGISTER(DiskPlanR)                                        \
  REGISTER(DiskMovingPlanR)                                  \
  REGISTER(CircleCircleR)                                    \
  REGISTER(DiskDiskR)                                        \
  REGISTER(SphereNEDSPlanR)                                  \
  REGISTER(SphereNEDSSphereNEDSR)                            \
  REGISTER(SphereLDSSphereLDSR)                              \
  REGISTER(SphereLDSPlanR)                                   \
  REGISTER(NonSmoothLaw)                                     \
  REGISTER(MixedComplementarityConditionNSL)                 \
  REGISTER(EqualityConditionNSL)                             \
  REGISTER(ComplementarityConditionNSL)                      \
  REGISTER(RelayNSL)                                         \
  REGISTER(NewtonImpactNSL)                                  \
  REGISTER(MultipleImpactNSL)                                \
  REGISTER(NewtonImpactFrictionNSL)                          \
  REGISTER(Simulation)                                       \
  REGISTER(TimeStepping)                                     \
  REGISTER(TimeSteppingProjectOnConstraints)                 \
  REGISTER(EventDriven)                                      \
  REGISTER(OneStepNSProblem)                                 \
  REGISTER(LCP)                                              \
  REGISTER(OSNSMultipleImpact)                               \
  REGISTER(FrictionContact)                                  \
  REGISTER(SimpleVector)                                     \
  REGISTER(BlockVector)                                      \
  REGISTER(LagrangianDS)                                     \
  REGISTER(LagrangianLinearTIDS)                             \
  REGISTER(FirstOrderLinearDS)                               \
  REGISTER(FirstOrderNonLinearDS)                            \
  REGISTER(FirstOrderLinearTIDS)                             \
  REGISTER(FirstOrderType1R)                                 \
  REGISTER(FirstOrderType2R)                                 \
  REGISTER(FirstOrderLinearR)                                \
  REGISTER(FirstOrderLinearTIR)                              \
  REGISTER(LagrangianScleronomousR)                          \
  REGISTER(LagrangianRheonomousR)                            \
  REGISTER(LagrangianCompliantR)                             \
  REGISTER(LagrangianLinearTIR)                              \
  REGISTER(NewtonEulerDS)                                    \
  REGISTER(NewtonEulerR)                                     \
  REGISTER_BASE(ExternalBody, LagrangianDS)                  \
  REGISTER_BASE(Lmgc2DDSK, LagrangianDS)                     \
  REGISTER_BASE(Lmgc2DPOLYG, LagrangianDS)                   \
  REGISTER_BASE(Disk, LagrangianDS)                          \
  REGISTER_BASE(Circle, LagrangianDS)                        \
  REGISTER_BASE(SphereLDS, LagrangianDS)                     \
  REGISTER_BASE(SphereNEDS, NewtonEulerDS)                   \
 
#endif
