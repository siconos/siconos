// Siconos-Front-End version 2.1.1, Copyright INRIA 2005-2007.
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.	
// Siconos is a free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// Siconos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Siconos; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
// Contact: Vincent ACARY vincent.acary@inrialpes.fr 
//	


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

#include "SiconosVector.h"
#include "SimpleVector.h"
#include "BlockVector.h"
#include "SiconosMatrix.h"
#include "SimpleMatrix.h"
#include "BlockMatrix.h"
#include "ioObject.h"
#include "ioMatrix.h"
#include "ioVector.h"

#include "SiconosMemory.h"
#include "SiconosMemoryException.h"

#include "SiconosDOMTreeTools.h"
#include "SiconosMemoryXML.h"
#include "RelationXML.h"
#include "FirstOrderRXML.h"
#include "FirstOrderLinearRXML.h"
#include "LagrangianRXML.h"
#include "LagrangianLinearRXML.h"
#include "NonSmoothLawXML.h"
#include "NewtonImpactNSLXML.h"
#include "NewtonImpactFrictionNSLXML.h"
#include "RelayNSLXML.h"
#include "ComplementarityConditionNSLXML.h"
#include "InteractionXML.h"
#include "DynamicalSystemXML.h"
#include "FirstOrderNonLinearDSXML.h"
#include "LagrangianDSXML.h"
#include "LagrangianLinearTIDSXML.h"
#include "FirstOrderLinearDSXML.h"
#include "NonSmoothDynamicalSystemXML.h"
#include "TimeDiscretisationXML.h"
#include "SolverXML.h"
#include "OneStepNSProblemXML.h"
#include "FrictionContactXML.h"
#include "QPXML.h"
#include "LCPXML.h"
#include "OneStepIntegratorXML.h"
#include "MoreauXML.h"
#include "LsodarXML.h"
#include "SimulationXML.h"
#include "SiconosModelXML.h"
#include "XMLTagsName.h"


#include "NonSmoothLaw.h"
#include "ComplementarityConditionNSL.h"
#include "NewtonImpactNSL.h"
#include "NewtonImpactFrictionNSL.h"
#include "RelayNSL.h"
#include "DynamicalSystem.h"
#include "FirstOrderNonLinearDS.h"
#include "FirstOrderLinearDS.h"
#include "FirstOrderLinearTIDS.h"
#include "LagrangianDS.h"
#include "LagrangianLinearTIDS.h"
#include "Interaction.h"
#include "Relation.h"
#include "FirstOrderR.h"
#include "FirstOrderLinearR.h"
#include "FirstOrderLinearTIR.h"
#include "LagrangianR.h"
#include "LagrangianLinearR.h"
#include "LagrangianScleronomousR.h"
#include "LagrangianRheonomousR.h"
#include "LagrangianCompliantR.h"
#include "NonSmoothDynamicalSystem.h"
#include "Topology.h"
#include "DynamicalSystemsSet.h"
#include "UnitaryRelation.h"
#include "UnitaryRelationsSet.h"
#include "InteractionsSet.h"
#include "TimeDiscretisation.h"
#include "Solver.h"
#include "Event.h"
#include "NonSmoothEvent.h"
#include "TimeDiscretisationEvent.h"
#include "EventsManager.h"
#include "OneStepNSProblem.h"
#include "QP.h"
#include "Relay.h"
#include "LCP.h"
#include "FrictionContact.h"
#include "FrictionContact2D.h"
#include "FrictionContact3D.h"
#include "OneStepIntegrator.h"
#include "Moreau.h"
#include "Lsodar.h"
#include "Simulation.h"
#include "EventDriven.h"
#include "TimeStepping.h"

#include "Model.h"
#include "SiconosConst.h"

#include "Sensor.h"
#include "SensorPosition.h"
#include "SensorEvent.h"
#include "Actuator.h"
#include "ExampleActuator.h"
#include "ActuatorEvent.h"
#include "ControlManager.h"

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
%include "SiconosVector.h"
%include "SimpleVector.h"
%include "BlockVector.h"
%include "SiconosMatrix.h"
%include "SimpleMatrix.h"
%include "BlockMatrix.h"
%include "ioObject.h"
%include "ioMatrix.h"
%include "ioVector.h"
%include "SiconosMemoryException.h"
%include "SiconosMemory.h"


// --- Xml ---
%include "SiconosDOMTreeTools.h"

// Xml - utils
%include "SiconosMemoryXML.h"

// Xml - formalisation
%include "RelationXML.h"
%include "FirstOrderRXML.h"
%include "FirstOrderLinearRXML.h"
%include "LagrangianRXML.h"
%include "LagrangianLinearRXML.h"
// ---
%include "NonSmoothLawXML.h"
%include "NewtonImpactNSLXML.h"
%include "NewtonImpactFrictionNSLXML.h"
%include "RelayNSLXML.h"
%include "ComplementarityConditionNSLXML.h"
// ---
%include "InteractionXML.h"
// --- 
%include "DynamicalSystemXML.h"
%include "FirstOrderNonLinearDSXML.h"
%include "LagrangianDSXML.h"
%include "LagrangianLinearTIDSXML.h"
%include "FirstOrderLinearDSXML.h"
// ---
%include "NonSmoothDynamicalSystemXML.h"

// Xml - simulation
%include "TimeDiscretisationXML.h"
%include "SolverXML.h"
// ---
%include "OneStepNSProblemXML.h"
%include "QPXML.h"
%include "LCPXML.h"
%include "FrictionContactXML.h"
// ---
%include "OneStepIntegratorXML.h"
%include "MoreauXML.h"
%include "LsodarXML.h"
// ---
%include "SimulationXML.h"

// Xml - model
%include "SiconosModelXML.h"
%include "XMLTagsName.h"


// --- ModelFormalisation ---

%include "Relation.h"
%include "FirstOrderR.h"
%include "FirstOrderLinearR.h"
%include "FirstOrderLinearTIR.h"
%include "LagrangianR.h"
%include "LagrangianScleronomousR.h"
%include "LagrangianRheonomousR.h"
%include "LagrangianCompliantR.h"
%include "LagrangianLinearR.h"
// ---
%include "NonSmoothLaw.h"
%include "RelayNSL.h"
%include "NewtonImpactNSL.h"
%include "NewtonImpactFrictionNSL.h"
%include "ComplementarityConditionNSL.h"
%include "Interaction.h"
// ---
%include "DynamicalSystem.h"
%include "FirstOrderNonLinearDS.h"
%include "FirstOrderLinearDS.h"
%include "FirstOrderLinearTIDS.h"
%include "LagrangianDS.h"
%include "LagrangianLinearTIDS.h"
// ---
%include "Topology.h"
%include "DynamicalSystemsSet.h"
%include "NonSmoothDynamicalSystem.h"
%include "UnitaryRelation.h"
%include "UnitaryRelationsSet.h"
%include "InteractionsSet.h"


// --- simulationTools ---
%include "TimeDiscretisation.h"
%include "Solver.h"
%include "Event.h"
%include "NonSmoothEvent.h"
%include "TimeDiscretisationEvent.h"
%include "EventsManager.h"

// ---
%include "OneStepNSProblem.h"
%include "QP.h"
%include "Relay.h"
%include "LCP.h"
%include "FrictionContact.h"
%include "FrictionContact2D.h"
%include "FrictionContact3D.h"
// ---
%include "OneStepIntegrator.h"
%include "Moreau.h"
%include "Lsodar.h"
// ---
%include "Simulation.h"
%include "EventDriven.h"
%include "TimeStepping.h"


// --- Model ---
%include "Model.h"
%include "SiconosConst.h"

%include "Sensor.h"
%include "SensorPosition.h"
%include "SensorEvent.h"
%include "Actuator.h"
%include "ExampleActuator.h"
%include "ActuatorEvent.h"
%include "ControlManager.h"

namespace std {
  %template(intVector) vector<int>;
  %template(doubleVector) vector<double>;
  %template(dsVector) vector<DynamicalSystem*>;
  %template(osiVector) vector<OneStepIntegrator*>;
};
