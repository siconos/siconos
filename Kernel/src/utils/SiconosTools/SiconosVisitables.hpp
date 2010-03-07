#ifndef SiconosVisitables_hpp
#define SiconosVisitables_hpp

#undef REGISTER
#define SICONOS_VISITABLES() \
    REGISTER(DynamicalSystem)\
    REGISTER(Disk)\
    REGISTER(Circle)\
    REGISTER(DiskPlanR)\
    REGISTER(DiskMovingPlanR)\
    REGISTER(CircleCircleR)\
    REGISTER(DiskDiskR)\
    REGISTER(SphereLDS)\
    REGISTER(SphereLDSSphereLDSR)\
    REGISTER(SphereLDSPlanR)\
    REGISTER(ExternalBody)\
    REGISTER(NonSmoothLaw)\
    REGISTER(MixedComplementarityConditionNSL)\
    REGISTER(EqualityConditionNSL)\
    REGISTER(ComplementarityConditionNSL)\
    REGISTER(RelayNSL)\
    REGISTER(NewtonImpactNSL)\
    REGISTER(NewtonImpactFrictionNSL)\
    REGISTER(Simulation)\
    REGISTER(TimeStepping)\
    REGISTER(EventDriven)\
    REGISTER(OneStepNSProblem)\
    REGISTER(LCP)\
    REGISTER(FrictionContact)\
    REGISTER(Lmgc2DDSK)\
    REGISTER(SimpleVector)\
    REGISTER(BlockVector)\
    REGISTER(LagrangianDS)\
    REGISTER(LagrangianLinearTIDS)\
    REGISTER(FirstOrderLinearDS)\
    REGISTER(FirstOrderNonLinearDS)\
    REGISTER(FirstOrderType1R)\
    REGISTER(FirstOrderLinearR)\
    REGISTER(FirstOrderLinearTIR)\
    REGISTER(LagrangianScleronomousR)\
    REGISTER(LagrangianRheonomousR)\
    REGISTER(LagrangianCompliantR)\
    REGISTER(LagrangianLinearTIR)\
 

#endif
