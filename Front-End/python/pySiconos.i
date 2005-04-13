// pySiconos.i - SWIG interface
%module pySiconos
%{
//#include "SiconosNumerics.h"

#include "SiconosException.h"
#include "RuntimeException.h"
#include "XMLException.h"

#include "SiconosSharedLibraryException.h"
#include "SiconosSharedLibrary.h"

#include "SiconosMatrixException.h"
#include "SiconosVectorException.h"
#include "CompositeVector.h"

#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"

#include "SiconosMemory.h"
#include "SiconosMemoryException.h"

#include "SiconosDOMTreeTools.h"
#include "SiconosMemoryXML.h"
#include "RelationXML.h"
#include "LagrangianLinearRXML.h"
#include "LagrangianNonLinearRXML.h"
#include "LinearTIRXML.h"
#include "NonSmoothLawXML.h"
#include "NewtonImpactLawNSLXML.h"
#include "NewtonImpactFrictionNSLXML.h"
#include "RelayNSLXML.h"
#include "ComplementarityConditionNSLXML.h"
#include "InteractionXML.h"
#include "BoundaryConditionXML.h"
#include "LinearBCXML.h"
#include "PeriodicBCXML.h"
#include "NLinearBCXML.h"
#include "DSXML.h"
#include "LagrangianDSXML.h"
#include "LagrangianLinearTIDSXML.h"
#include "LinearSystemDSXML.h"
#include "NSDSXML.h"
#include "TimeDiscretisationXML.h"
#include "OneStepNSProblemXML.h"
#include "QPXML.h"
#include "LCPXML.h"
#include "CFDXML.h"
#include "OneStepIntegratorXML.h"
#include "AdamsXML.h"
#include "MoreauXML.h"
#include "LsodarXML.h"
#include "StrategyXML.h"
#include "DSInputOutputXML.h"
#include "EqualityConstraintXML.h"
#include "LagrangianDSIOXML.h"
#include "LagrangianLinearDSIOXML.h"
#include "LagrangianLinearECXML.h"
#include "LagrangianECXML.h"
#include "LinearDSIOXML.h"
#include "LinearECXML.h"
#include "LinearTIECXML.h"
#include "SiconosModelXML.h"
#include "XMLTagsName.h"


#include "Relation.h"
#include "LagrangianLinearR.h"
#include "LagrangianNonLinearR.h"
#include "LinearTIR.h"
#include "NonSmoothLaw.h"
#include "RelayNSL.h"
#include "NewtonImpactLawNSL.h"
#include "NewtonImpactFrictionNSL.h"
#include "ComplementarityConditionNSL.h"
#include "BoundaryCondition.h"
#include "PeriodicBC.h"
#include "NLinearBC.h"
#include "LinearBC.h"
#include "Interaction.h"
#include "DynamicalSystem.h"
#include "LinearSystemDS.h"
#include "LagrangianDS.h"
#include "LagrangianLinearTIDS.h"
#include "DSInputOutput.h"
#include "LagrangianDSIO.h"
#include "LagrangianLinearDSIO.h"
#include "LinearDSIO.h"
#include "EqualityConstraint.h"
#include "LagrangianEC.h"
#include "LagrangianLinearEC.h"
#include "LinearEC.h"
#include "LinearTIEC.h"
#include "NonSmoothDynamicalSystem.h"


#include "TimeDiscretisation.h"
#include "OneStepNSProblem.h"
#include "Relay.h"
#include "LCP.h"
#include "CFD.h"
#include "QP.h"
#include "OneStepIntegrator.h"
#include "Moreau.h"
#include "Lsodar.h"
#include "Adams.h"
#include "Strategy.h"
#include "EventDriven.h"
#include "TimeStepping.h"


#include "KernelDefaultConfig.h"
#include "Model.h"
#include "SiconosConst.h"
%} 

%typemap(in) string  {
$1 = string(PyString_AsString($input));
}

# OPERATORS
%rename(assign) *::operator=;


// Parse the original header file

// strings
%include "std_string.i"
// vector of the C++ STL
%include "std_vector.i"



// -- Numerics ---
//%include "SiconosNumerics.h"

// --- Utils ---
%include "SiconosException.h"
%include "RuntimeException.h"
%include "XMLException.h"
%include "SiconosSharedLibraryException.h"
%include "SiconosSharedLibrary.h"
%include "SiconosMatrixException.h"
%include "SiconosVectorException.h"
%include "NewSiconosVector.h"
%include "SimpleVector.h"
%include "CompositeVector.h"
%include "SiconosMatrix.h"
%include "SiconosMemoryException.h"
%include "SiconosMemory.h"


// --- Xml ---
%include "SiconosDOMTreeTools.h"

// Xml - utils
%include "SiconosMemoryXML.h"

// Xml - formalisation
%include "RelationXML.h"
%include "LagrangianLinearRXML.h"
%include "LagrangianNonLinearRXML.h"
%include "LinearTIRXML.h"
// ---
%include "NonSmoothLawXML.h"
%include "NewtonImpactLawNSLXML.h"
%include "NewtonImpactFrictionNSLXML.h"
%include "RelayNSLXML.h"
%include "ComplementarityConditionNSLXML.h"
// ---
%include "InteractionXML.h"
// ---
%include "BoundaryConditionXML.h"
%include "LinearBCXML.h"
%include "PeriodicBCXML.h"
%include "NLinearBCXML.h"
// --- 
%include "DSXML.h"
%include "LagrangianDSXML.h"
%include "LagrangianLinearTIDSXML.h"
%include "LinearSystemDSXML.h"
// ---
%include "NSDSXML.h"

// Xml - strategy
%include "TimeDiscretisationXML.h"
// ---
%include "OneStepNSProblemXML.h"
%include "QPXML.h"
%include "LCPXML.h"
%include "CFDXML.h"
// ---
%include "OneStepIntegratorXML.h"
%include "AdamsXML.h"
%include "MoreauXML.h"
%include "LsodarXML.h"
// ---
%include "StrategyXML.h"

// Xml - EC and DSIO
%include "DSInputOutputXML.h"
%include "EqualityConstraintXML.h"
%include "LagrangianDSIOXML.h"
%include "LagrangianLinearDSIOXML.h"
%include "LagrangianLinearECXML.h"
%include "LagrangianECXML.h"
%include "LinearDSIOXML.h"
%include "LinearECXML.h"
%include "LinearTIECXML.h"

// Xml - model
%include "SiconosModelXML.h"
%include "XMLTagsName.h"


// --- ModelFormalisation ---

%include "Relation.h"
%include "LagrangianLinearR.h"
%include "LagrangianNonLinearR.h"
%include "LinearTIR.h"
// ---
%include "NonSmoothLaw.h"
%include "RelayNSL.h"
%include "NewtonImpactLawNSL.h"
%include "NewtonImpactFrictionNSL.h"
%include "ComplementarityConditionNSL.h"
// ---
%include "BoundaryCondition.h"
%include "PeriodicBC.h"
%include "NLinearBC.h"
%include "LinearBC.h"
// ---
%include "Interaction.h"
// ---
%include "DynamicalSystem.h"
%include "LinearSystemDS.h"
%include "LagrangianDS.h"
%include "LagrangianLinearTIDS.h"
// ---
%include "DSInputOutput.h"
%include "LagrangianDSIO.h"
%include "LagrangianLinearDSIO.h"
%include "LinearDSIO.h"
// ---
%include "EqualityConstraint.h"
%include "LagrangianEC.h"
%include "LagrangianLinearEC.h"
%include "LinearEC.h"
%include "LinearTIEC.h"
// ---
%include "NonSmoothDynamicalSystem.h"



// --- simulationTools ---
%include "TimeDiscretisation.h"
// ---
%include "OneStepNSProblem.h"
%include "Relay.h"
%include "LCP.h"
%include "CFD.h"
%include "QP.h"
// ---
%include "OneStepIntegrator.h"
%include "Moreau.h"
%include "Lsodar.h"
%include "Adams.h"
// ---
%include "Strategy.h"
%include "EventDriven.h"
%include "TimeStepping.h"


// --- Model ---
%include "KernelDefaultConfig.h"
%include "Model.h"
%include "SiconosConst.h"

namespace std {
  %template(intVector) vector<int>;
  %template(doubleVector) vector<double>;
  %template(dsioVector) vector<DSInputOutput*>;
  %template(dsVector) vector<DynamicalSystem*>;
  %template(ecVector) vector<EqualityConstraint*>;
  %template(osiVector) vector<OneStepIntegrator*>;
};
