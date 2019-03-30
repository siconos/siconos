#ifndef SiconosVisitables_hpp
#define SiconosVisitables_hpp

#undef REGISTER
#undef REGISTER_STRUCT
#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN

#define KERNEL_CLASSES()                               \
  REGISTER(DynamicalSystem)                            \
  REGISTER(Relation)                                   \
  REGISTER(NonSmoothLaw)                               \
  REGISTER(MixedComplementarityConditionNSL)           \
  REGISTER(EqualityConditionNSL)                       \
  REGISTER(ComplementarityConditionNSL)                \
  REGISTER(RelayNSL)                                   \
  REGISTER(NormalConeNSL)                              \
  REGISTER(NewtonImpactNSL)                            \
  REGISTER(MultipleImpactNSL)                          \
  REGISTER(NewtonImpactFrictionNSL)                    \
  REGISTER(NewtonImpactRollingFrictionNSL)             \
  REGISTER(Simulation)                                 \
  REGISTER(TimeStepping)                               \
  REGISTER(TimeSteppingD1Minus)                        \
  REGISTER(TimeSteppingDirectProjection)               \
  REGISTER(TimeSteppingCombinedProjection)             \
  REGISTER(EventDriven)                                \
  REGISTER(OneStepIntegrator)                          \
  REGISTER(EulerMoreauOSI)                             \
  REGISTER(MoreauJeanOSI)                              \
  REGISTER(MoreauJeanBilbaoOSI)                        \
  REGISTER(MoreauJeanGOSI)                             \
  REGISTER(MoreauJeanDirectProjectionOSI)              \
  REGISTER(MoreauJeanCombinedProjectionOSI)            \
  REGISTER(LsodarOSI)                                  \
  REGISTER(Hem5OSI)                                    \
  REGISTER(NewMarkAlphaOSI)                            \
  REGISTER(D1MinusLinearOSI)                           \
  REGISTER(SchatzmanPaoliOSI)                          \
  REGISTER(ZeroOrderHoldOSI)                           \
  REGISTER(OneStepNSProblem)                           \
  REGISTER(LinearOSNS)                                 \
  REGISTER(LCP)                                        \
  REGISTER(MLCP)                                       \
  REGISTER(MLCPProjectOnConstraints)                   \
  REGISTER(OSNSMultipleImpact)                         \
  REGISTER(FrictionContact)                            \
  REGISTER(GlobalFrictionContact)                      \
  REGISTER(SiconosVector)                              \
  REGISTER(SimpleMatrix)                               \
  REGISTER(BlockVector)                                \
  REGISTER(BlockMatrix)                                \
  REGISTER(SecondOrderDS)                              \
  REGISTER(LagrangianDS)                               \
  REGISTER(LagrangianLinearTIDS)                       \
  REGISTER(LagrangianLinearDiagonalDS)                 \
  REGISTER(FirstOrderLinearDS)                         \
  REGISTER(FirstOrderNonLinearDS)                      \
  REGISTER(FirstOrderLinearTIDS)                       \
  REGISTER(FirstOrderType1R)                           \
  REGISTER(FirstOrderType2R)                           \
  REGISTER(FirstOrderLinearR)                          \
  REGISTER(FirstOrderLinearTIR)                        \
  REGISTER(LagrangianScleronomousR)                    \
  REGISTER(LagrangianRheonomousR)                      \
  REGISTER(LagrangianCompliantR)                       \
  REGISTER(LagrangianLinearTIR)                        \
  REGISTER(NewtonEulerDS)                              \
  REGISTER(NewtonEulerR)                               \
  REGISTER_STRUCT(DynamicalSystemsGraph)               \
  REGISTER_STRUCT(InteractionsGraph)                   \
  REGISTER_STRUCT(DynamicalSystemsSubGraph)            \
  REGISTER_STRUCT(InteractionsSubGraph)

#ifndef SICONOS_VISITABLES
#define SICONOS_VISITABLES()                    \
  KERNEL_CLASSES()
#endif

#endif
