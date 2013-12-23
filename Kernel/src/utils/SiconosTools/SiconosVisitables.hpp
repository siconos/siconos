#ifndef SiconosVisitables_hpp
#define SiconosVisitables_hpp

#undef REGISTER
#undef REGISTER_STRUCT
#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN

#define KERNEL_CLASSES()                          \
  REGISTER(DynamicalSystem)                       \
  REGISTER(DynamicalSystemXML)                    \
    REGISTER(Relation)                            \
    REGISTER(NonSmoothLaw)                        \
    REGISTER(MixedComplementarityConditionNSL)    \
    REGISTER(EqualityConditionNSL)                \
    REGISTER(ComplementarityConditionNSL)         \
    REGISTER(RelayNSL)                            \
    REGISTER(NewtonImpactNSL)                     \
    REGISTER(MultipleImpactNSL)                   \
    REGISTER(NewtonImpactFrictionNSL)             \
    REGISTER(Simulation)                          \
    REGISTER(TimeStepping)                        \
    REGISTER(TimeSteppingD1Minus)                 \
    REGISTER(TimeSteppingProjectOnConstraints)    \
    REGISTER(TimeSteppingCombinedProjection)      \
    REGISTER(EventDriven)                         \
    REGISTER(OneStepIntegrator)                   \
    REGISTER(MoreauJeanOSI)                              \
    REGISTER(MoreauJeanDirectProjectionOSI)       \
    REGISTER(MoreauJeanCombinedProjectionOSI)         \
    REGISTER(Lsodar)                              \
    REGISTER(Hem5)                                \
    REGISTER(NewMarkAlphaOSI)                     \
    REGISTER(D1MinusLinear)                       \
    REGISTER(SchatzmanPaoli)                       \
    REGISTER(ZeroOrderHold)                        \
    REGISTER(OneStepNSProblem)                     \
    REGISTER(LinearOSNS)                           \
    REGISTER(LCP)                                  \
    REGISTER(MLCP)                                 \
    REGISTER(MLCPProjectOnConstraints)             \
    REGISTER(OSNSMultipleImpact)                   \
    REGISTER(FrictionContact)                      \
    REGISTER(SiconosVector)                        \
    REGISTER(SimpleMatrix)                         \
    REGISTER(BlockVector)                          \
    REGISTER(BlockMatrix)                          \
    REGISTER(LagrangianDS)                         \
    REGISTER(LagrangianLinearTIDS)                 \
    REGISTER(FirstOrderLinearDS)                   \
    REGISTER(FirstOrderNonLinearDS)                \
    REGISTER(FirstOrderLinearTIDS)                 \
    REGISTER(FirstOrderType1R)                     \
    REGISTER(FirstOrderType2R)                     \
    REGISTER(FirstOrderLinearR)                    \
    REGISTER(FirstOrderLinearTIR)                  \
    REGISTER(LagrangianScleronomousR)              \
    REGISTER(LagrangianRheonomousR)                \
    REGISTER(LagrangianCompliantR)                 \
    REGISTER(LagrangianLinearTIR)                  \
    REGISTER(NewtonEulerDS)                        \
    REGISTER(NewtonEulerR)                         \
    REGISTER_STRUCT(DynamicalSystemsGraph)         \
    REGISTER_STRUCT(InteractionsGraph)             \
    REGISTER_STRUCT(DynamicalSystemsSubGraph)      \
    REGISTER_STRUCT(InteractionsSubGraph)          \


#define MECHANICS_CLASSES()                       \
  REGISTER(SpaceFilter)                           \
  REGISTER(DiskPlanR)                             \
  REGISTER(DiskMovingPlanR)                       \
    REGISTER(CircleCircleR)                       \
    REGISTER(DiskDiskR)                           \
    REGISTER(SphereNEDSPlanR)                     \
    REGISTER(SphereNEDSSphereNEDSR)               \
    REGISTER(SphereLDSSphereLDSR)                 \
    REGISTER(SphereLDSPlanR)                      \
    REGISTER_BASE(Disk, LagrangianDS)             \
    REGISTER_BASE(CircularDS, LagrangianDS)       \
    REGISTER_BASE(Circle, LagrangianDS)           \
    REGISTER_BASE(SphereLDS, LagrangianDS)        \
    REGISTER_BASE(SphereNEDS, NewtonEulerDS)      \
    
#define BULLET_CLASSES()                          \
  REGISTER(BulletSpaceFilter)                     \
    REGISTER(BulletR)                             \
    REGISTER_BASE(BulletDS, NewtonEulerDS)  


#define LMGC_CLASSES()                            \
  REGISTER_BASE_EXTERN(Lmgc2DDS, LagrangianDS)    \
  REGISTER_BASE_EXTERN(Lmgc2DDSK, LagrangianDS)   \
  REGISTER_BASE_EXTERN(Lmgc2DPOLYG, LagrangianDS) \
  REGISTER(Lmgc2DR)                               

#define SICONOS_VISITABLES()                    \
  KERNEL_CLASSES()                              \
  MECHANICS_CLASSES()                           \
  BULLET_CLASSES()                              \
  LMGC_CLASSES()


#endif
