#ifndef SiconosVisitables_hpp
#define SiconosVisitables_hpp

#undef REGISTER
#define SICONOS_VISITABLES()                      \
  REGISTER(DynamicalSystem)                       \
  REGISTER(DynamicalSystemXML)                    \
  REGISTER(Relation)                              \
  REGISTER(DiskPlanR)                             \
  REGISTER(DiskMovingPlanR)                       \
  REGISTER(CircleCircleR)                         \
  REGISTER(DiskDiskR)                             \
  REGISTER(SphereNEDSPlanR)                       \
  REGISTER(SphereNEDSSphereNEDSR)                 \
  REGISTER(SphereLDSSphereLDSR)                   \
  REGISTER(SphereLDSPlanR)                        \
  REGISTER(NonSmoothLaw)                          \
  REGISTER(MixedComplementarityConditionNSL)      \
  REGISTER(EqualityConditionNSL)                  \
  REGISTER(ComplementarityConditionNSL)           \
  REGISTER(RelayNSL)                              \
  REGISTER(NewtonImpactNSL)                       \
  REGISTER(MultipleImpactNSL)                     \
  REGISTER(NewtonImpactFrictionNSL)               \
  REGISTER(Simulation)                            \
  REGISTER(TimeStepping)                          \
  REGISTER(TimeSteppingD1Minus)                   \
  REGISTER(TimeSteppingProjectOnConstraints)      \
  REGISTER(TimeSteppingCombinedProjection)        \
  REGISTER(EventDriven)                           \
  REGISTER(OneStepIntegrator)                     \
  REGISTER(Moreau)                                \
  REGISTER(MoreauCombinedProjectionOSI)           \
  REGISTER(Lsodar)                                \
  REGISTER(D1MinusLinear)                         \
  REGISTER(SchatzmanPaoli)                        \
  REGISTER(OneStepNSProblem)                      \
  REGISTER(LinearOSNS)                            \
  REGISTER(LCP)                                   \
  REGISTER(MLCP)                                  \
  REGISTER(MLCPProjectOnConstraints)              \
  REGISTER(OSNSMultipleImpact)                    \
  REGISTER(FrictionContact)                       \
  REGISTER(SimpleVector)                          \
  REGISTER(SimpleMatrix)                          \
  REGISTER(BlockVector)                           \
  REGISTER(BlockMatrix)                           \
  REGISTER(LagrangianDS)                          \
  REGISTER(LagrangianLinearTIDS)                  \
  REGISTER(FirstOrderLinearDS)                    \
  REGISTER(FirstOrderNonLinearDS)                 \
  REGISTER(FirstOrderLinearTIDS)                  \
  REGISTER(FirstOrderType1R)                      \
  REGISTER(FirstOrderType2R)                      \
  REGISTER(FirstOrderLinearR)                     \
  REGISTER(FirstOrderLinearTIR)                   \
  REGISTER(LagrangianScleronomousR)               \
  REGISTER(LagrangianRheonomousR)                 \
  REGISTER(LagrangianCompliantR)                  \
  REGISTER(LagrangianLinearTIR)                   \
  REGISTER(NewtonEulerDS)                         \
  REGISTER(NewtonEulerR)                          \
  REGISTER(BulletR)                               \
  REGISTER(BulletRImpact)                         \
  REGISTER(Lmgc2DR)                               \
  REGISTER(SpaceFilter)                           \
  REGISTER(BulletSpaceFilter)                     \
  REGISTER_STRUCT(DynamicalSystemsGraph)          \
  REGISTER_STRUCT(InteractionsGraph)              \
  REGISTER_STRUCT(DynamicalSystemsSubGraph)       \
  REGISTER_STRUCT(InteractionsSubGraph)           \
  REGISTER_BASE(ExternalBody, LagrangianDS)       \
  REGISTER_BASE_EXTERN(Lmgc2DDS, LagrangianDS)    \
  REGISTER_BASE_EXTERN(Lmgc2DDSK, LagrangianDS)   \
  REGISTER_BASE_EXTERN(Lmgc2DPOLYG, LagrangianDS) \
  REGISTER_BASE(Disk, LagrangianDS)               \
  REGISTER_BASE(CircularDS, LagrangianDS)         \
  REGISTER_BASE(Circle, LagrangianDS)             \
  REGISTER_BASE(SphereLDS, LagrangianDS)          \
  REGISTER_BASE(SphereNEDS, NewtonEulerDS)        \
  REGISTER_BASE(BulletDS, NewtonEulerDS)

#endif
