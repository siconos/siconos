// pySiconos.i - SWIG interface
%module pySiconos
%{
#include "../src/utils/SiconosException/SiconosException.h"
#include "../src/utils/SiconosException/RuntimeException.h"
#include "../src/utils/SiconosException/XMLException.h"

#include "../src/utils/SiconosSharedLibrary/SiconosSharedLibraryException.h"
#include "../src/utils/SiconosSharedLibrary/SiconosSharedLibrary.h"

#include "../src/utils/NewSiconosVector/SiconosMatrixException.h"
#include "../src/utils/NewSiconosVector/SiconosVectorException.h"
#include "../src/utils/NewSiconosVector/CompositeVector.h"

#include "../src/utils/NewSiconosVector/NewSiconosVector.h"
#include "../src/utils/NewSiconosVector/SimpleVector.h"
#include "../src/utils/NewSiconosVector/SiconosMatrix.h"

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

#include "../src/modelformalisation/Relation.h"
#include "../src/modelformalisation/LagrangianLinearR.h"
#include "../src/modelformalisation/LagrangianNonLinearR.h"
#include "../src/modelformalisation/LinearTIR.h"
#include "../src/modelformalisation/NonSmoothLaw.h"
#include "../src/modelformalisation/RelayNSL.h"
#include "../src/modelformalisation/NewtonImpactLawNSL.h"
#include "../src/modelformalisation/NewtonImpactFrictionNSL.h"
#include "../src/modelformalisation/ComplementarityConditionNSL.h"
#include "../src/modelformalisation/BoundaryCondition.h"
#include "../src/modelformalisation/PeriodicBC.h"
#include "../src/modelformalisation/NLinearBC.h"
#include "../src/modelformalisation/LinearBC.h"
#include "../src/modelformalisation/Interaction.h"
#include "../src/modelformalisation/DynamicalSystem.h"
#include "../src/modelformalisation/LinearSystemDS.h"
#include "../src/modelformalisation/LagrangianNLDS.h"
#include "../src/modelformalisation/LagrangianTIDS.h"
#include "../src/modelformalisation/DSInputOutput.h"
#include "../src/modelformalisation/LagrangianDSIO.h"
#include "../src/modelformalisation/LinearDSIO.h"
#include "../src/modelformalisation/EqualityConstraint.h"
#include "../src/modelformalisation/LagrangianEC.h"
#include "../src/modelformalisation/LinearEC.h"
#include "../src/modelformalisation/LinearTIEC.h"
#include "../src/modelformalisation/NonSmoothDynamicalSystem.h"


#include "../src/modelstrategy/TimeDiscretisation.h"
#include "../src/modelstrategy/OneStepNSProblem.h"
#include "../src/modelstrategy/Relay.h"
#include "../src/modelstrategy/LCP.h"
#include "../src/modelstrategy/QP.h"
#include "../src/modelstrategy/OneStepIntegrator.h"
#include "../src/modelstrategy/Moreau.h"
#include "../src/modelstrategy/Lsodar.h"
#include "../src/modelstrategy/Adams.h"
#include "../src/modelstrategy/Strategy.h"
#include "../src/modelstrategy/EventDriven.h"
#include "../src/modelstrategy/TimeStepping.h"


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
%include "../src/utils/NewSiconosVector/SiconosMatrixException.h"
%include "../src/utils/NewSiconosVector/SiconosVectorException.h"
%include "../src/utils/NewSiconosVector/NewSiconosVector.h"
%include "../src/utils/NewSiconosVector/SimpleVector.h"
%include "../src/utils/NewSiconosVector/CompositeVector.h"
%include "../src/utils/NewSiconosVector/SiconosMatrix.h"
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

%include "../src/modelformalisation/Relation.h"
%include "../src/modelformalisation/LagrangianLinearR.h"
%include "../src/modelformalisation/LagrangianNonLinearR.h"
%include "../src/modelformalisation/LinearTIR.h"
// ---
%include "../src/modelformalisation/NonSmoothLaw.h"
%include "../src/modelformalisation/RelayNSL.h"
%include "../src/modelformalisation/NewtonImpactLawNSL.h"
%include "../src/modelformalisation/NewtonImpactFrictionNSL.h"
%include "../src/modelformalisation/ComplementarityConditionNSL.h"
// ---
%include "../src/modelformalisation/BoundaryCondition.h"
%include "../src/modelformalisation/PeriodicBC.h"
%include "../src/modelformalisation/NLinearBC.h"
%include "../src/modelformalisation/LinearBC.h"
// ---
%include "../src/modelformalisation/Interaction.h"
// ---
%include "../src/modelformalisation/DynamicalSystem.h"
%include "../src/modelformalisation/LinearSystemDS.h"
%include "../src/modelformalisation/LagrangianNLDS.h"
%include "../src/modelformalisation/LagrangianTIDS.h"
// ---
%include "../src/modelformalisation/DSInputOutput.h"
%include "../src/modelformalisation/LagrangianDSIO.h"
%include "../src/modelformalisation/LinearDSIO.h"
// ---
%include "../src/modelformalisation/EqualityConstraint.h"
%include "../src/modelformalisation/LagrangianEC.h"
%include "../src/modelformalisation/LinearEC.h"
%include "../src/modelformalisation/LinearTIEC.h"
// ---
%include "../src/modelformalisation/NonSmoothDynamicalSystem.h"



// --- ModelStrategy ---
%include "../src/modelstrategy/TimeDiscretisation.h"
// ---
%include "../src/modelstrategy/OneStepNSProblem.h"
%include "../src/modelstrategy/Relay.h"
%include "../src/modelstrategy/LCP.h"
%include "../src/modelstrategy/QP.h"
// ---
%include "../src/modelstrategy/OneStepIntegrator.h"
%include "../src/modelstrategy/Moreau.h"
%include "../src/modelstrategy/Lsodar.h"
%include "../src/modelstrategy/Adams.h"
// ---
%include "../src/modelstrategy/Strategy.h"
%include "../src/modelstrategy/EventDriven.h"
%include "../src/modelstrategy/TimeStepping.h"


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
