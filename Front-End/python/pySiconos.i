// pySiconos.i - SWIG interface
%module pySiconos
%{
#include "../src/utils/SiconosException/SiconosException.h"
#include "../src/utils/SiconosException/RuntimeException.h"
#include "../src/utils/SiconosException/XMLException.h"

#include "../src/utils/SiconosSharedLibrary/SiconosSharedLibraryException.h"
#include "../src/utils/SiconosSharedLibrary/SiconosSharedLibrary.h"

#include "../src/utils/SiconosAlgebra/SiconosMatrixException.h"
#include "../src/utils/SiconosAlgebra/SiconosVectorException.h"
#include "../src/utils/SiconosAlgebra/CompositeVector.h"

#include "../src/utils/SiconosAlgebra/NewSiconosVector.h"
#include "../src/utils/SiconosAlgebra/SimpleVector.h"
#include "../src/utils/SiconosAlgebra/SiconosMatrix.h"

#include "../src/utils/SiconosMemory/SiconosMemory.h"
#include "../src/utils/SiconosMemory/SiconosMemoryException.h"


#include "../src/xml/SiconosDOMTreeTools.h"
#include "../src/xml/SiconosMemoryXML.h"
#include "../src/xml/RelationXML.h"
#include "../src/xml/LagrangianLinearRXML.h"
#include "../src/xml/LagrangianNonLinearRXML.h"
#include "../src/xml/LinearTIRXML.h"
#include "../src/xml/NonSmoothLawXML.h"
#include "../src/xml/NewtonImpactLawNSLXML.h"
#include "../src/xml/NewtonImpactFrictionNSLXML.h"
#include "../src/xml/RelayNSLXML.h"
#include "../src/xml/ComplementarityConditionNSLXML.h"
#include "../src/xml/InteractionXML.h"
#include "../src/xml/BoundaryConditionXML.h"
#include "../src/xml/LinearBCXML.h"
#include "../src/xml/PeriodicBCXML.h"
#include "../src/xml/NLinearBCXML.h"
#include "../src/xml/DSXML.h"
#include "../src/xml/LagrangianNLDSXML.h"
#include "../src/xml/LagrangianTIDSXML.h"
#include "../src/xml/LinearSystemDSXML.h"
#include "../src/xml/NSDSXML.h"
#include "../src/xml/TimeDiscretisationXML.h"
#include "../src/xml/OneStepNSProblemXML.h"
#include "../src/xml/QPXML.h"
#include "../src/xml/LCPXML.h"
#include "../src/xml/OneStepIntegratorXML.h"
#include "../src/xml/AdamsXML.h"
#include "../src/xml/MoreauXML.h"
#include "../src/xml/LsodarXML.h"
#include "../src/xml/StrategyXML.h"
#include "../src/xml/DSInputOutputXML.h"
#include "../src/xml/EqualityConstraintXML.h"
#include "../src/xml/LagrangianDSIOXML.h"
#include "../src/xml/LagrangianECXML.h"
#include "../src/xml/LinearDSIOXML.h"
#include "../src/xml/LinearECXML.h"
#include "../src/xml/LinearTIECXML.h"
#include "../src/xml/SiconosModelXML.h"
#include "../src/xml/XMLTagsName.h"

#include "../src/modelingTools/Relation.h"
#include "../src/modelingTools/LagrangianLinearR.h"
#include "../src/modelingTools/LagrangianNonLinearR.h"
#include "../src/modelingTools/LinearTIR.h"
#include "../src/modelingTools/NonSmoothLaw.h"
#include "../src/modelingTools/RelayNSL.h"
#include "../src/modelingTools/NewtonImpactLawNSL.h"
#include "../src/modelingTools/NewtonImpactFrictionNSL.h"
#include "../src/modelingTools/ComplementarityConditionNSL.h"
#include "../src/modelingTools/BoundaryCondition.h"
#include "../src/modelingTools/PeriodicBC.h"
#include "../src/modelingTools/NLinearBC.h"
#include "../src/modelingTools/LinearBC.h"
#include "../src/modelingTools/Interaction.h"
#include "../src/modelingTools/DynamicalSystem.h"
#include "../src/modelingTools/LinearSystemDS.h"
#include "../src/modelingTools/LagrangianNLDS.h"
#include "../src/modelingTools/LagrangianTIDS.h"
#include "../src/modelingTools/DSInputOutput.h"
#include "../src/modelingTools/LagrangianDSIO.h"
#include "../src/modelingTools/LinearDSIO.h"
#include "../src/modelingTools/EqualityConstraint.h"
#include "../src/modelingTools/LagrangianEC.h"
#include "../src/modelingTools/LinearEC.h"
#include "../src/modelingTools/LinearTIEC.h"
#include "../src/modelingTools/NonSmoothDynamicalSystem.h"


#include "../src/simulationTools/TimeDiscretisation.h"
#include "../src/simulationTools/OneStepNSProblem.h"
#include "../src/simulationTools/Relay.h"
#include "../src/simulationTools/LCP.h"
#include "../src/simulationTools/QP.h"
#include "../src/simulationTools/OneStepIntegrator.h"
#include "../src/simulationTools/Moreau.h"
#include "../src/simulationTools/Lsodar.h"
#include "../src/simulationTools/Adams.h"
#include "../src/simulationTools/Strategy.h"
#include "../src/simulationTools/EventDriven.h"
#include "../src/simulationTools/TimeStepping.h"


#include "../src/model/KernelDefaultConfig.h"
#include "../src/model/Model.h"
#include "../src/model/SiconosConst.h"
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





//namespace std {
//   %template(intVector) vector<int>;
//   %template(doubleVector) vector<double>;
//};


// --- Utils ---

%include "../src/utils/SiconosException/SiconosException.h"
%include "../src/utils/SiconosException/RuntimeException.h"
%include "../src/utils/SiconosException/XMLException.h"
%include "../src/utils/SiconosSharedLibrary/SiconosSharedLibraryException.h"
%include "../src/utils/SiconosSharedLibrary/SiconosSharedLibrary.h"
%include "../src/utils/SiconosAlgebra/SiconosMatrixException.h"
%include "../src/utils/SiconosAlgebra/SiconosVectorException.h"
%include "../src/utils/SiconosAlgebra/NewSiconosVector.h"
%include "../src/utils/SiconosAlgebra/SimpleVector.h"
%include "../src/utils/SiconosAlgebra/CompositeVector.h"
%include "../src/utils/SiconosAlgebra/SiconosMatrix.h"
%include "../src/utils/SiconosMemory/SiconosMemoryException.h"
%include "../src/utils/SiconosMemory/SiconosMemory.h"


// --- Xml ---

//--------
%include "../src/xml/SiconosDOMTreeTools.h"

// Xml - utils
%include "../src/xml/SiconosMemoryXML.h"

// Xml - formalisation
%include "../src/xml/RelationXML.h"
%include "../src/xml/LagrangianLinearRXML.h"
%include "../src/xml/LagrangianNonLinearRXML.h"
%include "../src/xml/LinearTIRXML.h"
// ---
%include "../src/xml/NonSmoothLawXML.h"
%include "../src/xml/NewtonImpactLawNSLXML.h"
%include "../src/xml/NewtonImpactFrictionNSLXML.h"
%include "../src/xml/RelayNSLXML.h"
%include "../src/xml/ComplementarityConditionNSLXML.h"
// ---
%include "../src/xml/InteractionXML.h"
// ---
%include "../src/xml/BoundaryConditionXML.h"
%include "../src/xml/LinearBCXML.h"
%include "../src/xml/PeriodicBCXML.h"
%include "../src/xml/NLinearBCXML.h"
// --- 
%include "../src/xml/DSXML.h"
%include "../src/xml/LagrangianNLDSXML.h"
%include "../src/xml/LagrangianTIDSXML.h"
%include "../src/xml/LinearSystemDSXML.h"
// ---
%include "../src/xml/NSDSXML.h"

// Xml - strategy
%include "../src/xml/TimeDiscretisationXML.h"
// ---
%include "../src/xml/OneStepNSProblemXML.h"
%include "../src/xml/QPXML.h"
%include "../src/xml/LCPXML.h"
// ---
%include "../src/xml/OneStepIntegratorXML.h"
%include "../src/xml/AdamsXML.h"
%include "../src/xml/MoreauXML.h"
%include "../src/xml/LsodarXML.h"
// ---
%include "../src/xml/StrategyXML.h"

// Xml - EC and DSIO
%include "../src/xml/DSInputOutputXML.h"
%include "../src/xml/EqualityConstraintXML.h"
%include "../src/xml/LagrangianDSIOXML.h"
%include "../src/xml/LagrangianECXML.h"
%include "../src/xml/LinearDSIOXML.h"
%include "../src/xml/LinearECXML.h"
%include "../src/xml/LinearTIECXML.h"

// Xml - model
%include "../src/xml/SiconosModelXML.h"
%include "../src/xml/XMLTagsName.h"


// --- ModelFormalisation ---

%include "../src/modelingTools/Relation.h"
%include "../src/modelingTools/LagrangianLinearR.h"
%include "../src/modelingTools/LagrangianNonLinearR.h"
%include "../src/modelingTools/LinearTIR.h"
// ---
%include "../src/modelingTools/NonSmoothLaw.h"
%include "../src/modelingTools/RelayNSL.h"
%include "../src/modelingTools/NewtonImpactLawNSL.h"
%include "../src/modelingTools/NewtonImpactFrictionNSL.h"
%include "../src/modelingTools/ComplementarityConditionNSL.h"
// ---
%include "../src/modelingTools/BoundaryCondition.h"
%include "../src/modelingTools/PeriodicBC.h"
%include "../src/modelingTools/NLinearBC.h"
%include "../src/modelingTools/LinearBC.h"
// ---
%include "../src/modelingTools/Interaction.h"
// ---
%include "../src/modelingTools/DynamicalSystem.h"
%include "../src/modelingTools/LinearSystemDS.h"
%include "../src/modelingTools/LagrangianNLDS.h"
%include "../src/modelingTools/LagrangianTIDS.h"
// ---
%include "../src/modelingTools/DSInputOutput.h"
%include "../src/modelingTools/LagrangianDSIO.h"
%include "../src/modelingTools/LinearDSIO.h"
// ---
%include "../src/modelingTools/EqualityConstraint.h"
%include "../src/modelingTools/LagrangianEC.h"
%include "../src/modelingTools/LinearEC.h"
%include "../src/modelingTools/LinearTIEC.h"
// ---
%include "../src/modelingTools/NonSmoothDynamicalSystem.h"



// --- simulationTools ---
%include "../src/simulationTools/TimeDiscretisation.h"
// ---
%include "../src/simulationTools/OneStepNSProblem.h"
%include "../src/simulationTools/Relay.h"
%include "../src/simulationTools/LCP.h"
%include "../src/simulationTools/QP.h"
// ---
%include "../src/simulationTools/OneStepIntegrator.h"
%include "../src/simulationTools/Moreau.h"
%include "../src/simulationTools/Lsodar.h"
%include "../src/simulationTools/Adams.h"
// ---
%include "../src/simulationTools/Strategy.h"
%include "../src/simulationTools/EventDriven.h"
%include "../src/simulationTools/TimeStepping.h"


// --- Model ---
%include "../src/model/KernelDefaultConfig.h"
%include "../src/model/Model.h"
%include "../src/model/SiconosConst.h"

namespace std {
   %template(intVector) vector<int>;
  %template(doubleVector) vector<double>;
  %template(dsioVector) vector<DSInputOutput*>;
   %template(dsVector) vector<DynamicalSystem*>;
   %template(ecVector) vector<EqualityConstraint*>;
   %template(osiVector) vector<OneStepIntegrator*>;
};
