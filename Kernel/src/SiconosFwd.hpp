#ifndef SiconosFwd_hpp
#define SiconosFwd_hpp

#include "SiconosSerialization.hpp"
#include "SiconosPointers.hpp"

/* Forward declarations */

// --- Numerics ---
TYPEDEF_SPTR(MixedLinearComplementarityProblem)
TYPEDEF_SPTR(NumericsOptions)
TYPEDEF_SPTR(SolverOptions)
TYPEDEF_SPTR(NumericsMatrix)
// ----------------

/* Kernel */
DEFINE_SPTR(BlockCSRMatrix)

DEFINE_SPTR(Interaction)

DEFINE_SPTR(Model)

DEFINE_SPTR(NonSmoothDynamicalSystem)

// --- Non-Smooth problems ---
DEFINE_SPTR(OneStepNSProblem)
DEFINE_SPTR(QP)
DEFINE_SPTR(LinearOSNS)
DEFINE_SPTR(FrictionContact)
DEFINE_SPTR(LCP)
DEFINE_SPTR(MLCP)
DEFINE_SPTR(MLCPProjectOnConstraints)
DEFINE_SPTR(Relay)
DEFINE_SPTR(Equality)
DEFINE_SPTR(GenericMechanical)
DEFINE_SPTR(OSNSMultipleImpact)
// ----------------------------


DEFINE_SPTR(OneStepIntegrator)

DEFINE_SPTR(Relation)

DEFINE_SPTR(Simulation)
DEFINE_SPTR(EventDriven)
DEFINE_SPTR(TimeStepping)
DEFINE_SPTR(EventsManager)

DEFINE_SPTR(RelayNSL)
DEFINE_SPTR(MixedComplementarityConditionNSL)

DEFINE_SPTR(TimeDiscretisation)

// Dynamical systems
DEFINE_SPTR(DynamicalSystem)
DEFINE_SPTR(LagrangianLinearTIDS)
DEFINE_SPTR(NewtonEulerDS)

DEFINE_SPTR(Event)
DEFINE_SPTR(NonSmoothLaw)
DEFINE_SPTR(DynamicalSystemsSet)

DEFINE_SPTR(MatrixIntegrator)
DEFINE_SPTR(PluggedObject)
DEFINE_SPTR(SubPluggedObject)

DEFINE_SPTR(SiconosMatrix)
DEFINE_SPTR(SimpleMatrix)
DEFINE_SPTR(BlockMatrix)
DEFINE_SPTR(SiconosVector)
DEFINE_SPTR(BlockVector)

DEFINE_SPTR(OSNSMatrix)

DEFINE_SPTR(SiconosMemory)

DEFINE_SPTR(NewtonEulerR)

TYPEDEF_SAPTR(integer)
TYPEDEF_SPTR(integer)
TYPEDEF_SAPTR(doublereal)
TYPEDEF_SPTR(doublereal)

#endif
