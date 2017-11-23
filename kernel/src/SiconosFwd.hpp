#ifndef SiconosFwd_hpp
#define SiconosFwd_hpp

#include "SiconosSerialization.hpp"
#include "SiconosPointers.hpp"

/* Forward declarations */

// --- Numerics ---

#include <MixedLinearComplementarityProblem.h>
#include <SolverOptions.h>
#include <NumericsMatrix.h>

TYPEDEF_SPTR(MixedLinearComplementarityProblem)
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
DEFINE_SPTR(AVI)
DEFINE_SPTR(MLCP)
DEFINE_SPTR(MLCPProjectOnConstraints)
DEFINE_SPTR(Relay)
DEFINE_SPTR(Equality)
DEFINE_SPTR(GenericMechanical)
DEFINE_SPTR(OSNSMultipleImpact)
// ----------------------------


DEFINE_SPTR(OneStepIntegrator)

DEFINE_SPTR(Relation)

DEFINE_SPTR(FirstOrderR)

DEFINE_SPTR(Simulation)
DEFINE_SPTR(EventDriven)
DEFINE_SPTR(TimeStepping)
DEFINE_SPTR(EventsManager)
DEFINE_SPTR(InteractionManager)

DEFINE_SPTR(RelayNSL)
DEFINE_SPTR(MixedComplementarityConditionNSL)
DEFINE_SPTR(NormalConeNSL)

DEFINE_SPTR(TimeDiscretisation)

// Dynamical systems
DEFINE_SPTR(DynamicalSystem)
DEFINE_SPTR(LagrangianLinearTIDS)
DEFINE_SPTR(LagrangianLinearDiagonalDS)
DEFINE_SPTR(NewtonEulerDS)

DEFINE_SPTR(Event)
DEFINE_SPTR(NonSmoothLaw)

DEFINE_SPTR(MatrixIntegrator)
DEFINE_SPTR(PluggedObject)
DEFINE_SPTR(SubPluggedObject)
DEFINE_SPTR_STRUCT(ExtraAdditionalTerms)

DEFINE_SPTR(SiconosMatrix)
DEFINE_SPTR(SimpleMatrix)
DEFINE_SPTR(BlockMatrix)
DEFINE_SPTR(SiconosVector)
DEFINE_SPTR(BlockVector)

DEFINE_SPTR(OSNSMatrix)

DEFINE_SPTR(SiconosMemory)
#include <vector>
typedef std::vector<SP::SiconosMemory> VectorOfMemories;

DEFINE_SPTR(NewtonEulerR)
DEFINE_SPTR(NewtonEulerFrom1DLocalFrameR)
DEFINE_SPTR(NewtonEulerFrom3DLocalFrameR)

// OSI
DEFINE_SPTR(EulerMoreauOSI)
DEFINE_SPTR(MoreauJeanOSI)
DEFINE_SPTR(MoreauJeanBilbaoOSI)
DEFINE_SPTR(MoreauJeanGOSI)
DEFINE_SPTR(MoreauJeanCombinedProjectionOSI)
DEFINE_SPTR(MoreauJeanDirectProjectionOSI)
DEFINE_SPTR(LsodarOSI)
DEFINE_SPTR(Hem5OSI)
DEFINE_SPTR(D1MinusLinearOSI)
DEFINE_SPTR(SchatzmanPaoliOSI)
DEFINE_SPTR(ZeroOrderHoldOSI)
DEFINE_SPTR(NewMarkAlphaOSI);


// Graph thing
DEFINE_SPTR_STRUCT(InteractionProperties)
DEFINE_SPTR_STRUCT(GraphProperties)
DEFINE_SPTR_STRUCT(DynamicalSystemsGraph)
DEFINE_SPTR_STRUCT(InteractionsGraph)

#ifndef _F2C_INCLUDE_H
typedef int integer;
typedef double doublereal;
#endif

TYPEDEF_SAPTR(integer)
TYPEDEF_SAPTR(doublereal)
TYPEDEF_SPTR(integer)
TYPEDEF_SPTR(doublereal)

#endif
