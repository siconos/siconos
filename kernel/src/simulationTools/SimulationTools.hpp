/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
/*! \file SimulationTools.hpp
Include files related to simulationTools
Note that not all files from the current location are listed below, since some of them are already included inside the ones below.
*/

//#include "Interaction.hpp"
#include "Topology.hpp"
#include "EventDriven.hpp"
#include "EventsManager.hpp"
#include "EventFactory.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"
#include "TimeSteppingD1Minus.hpp"
#include "TimeSteppingDirectProjection.hpp"
#include "TimeSteppingCombinedProjection.hpp"
#include "InteractionManager.hpp"

#include "Equality.hpp"
#include "LCP.hpp"
#include "OSNSMultipleImpact.hpp"
#include "MLCP.hpp"
#include "MLCPProjectOnConstraints.hpp"
#include "AVI.hpp"
#include "GenericMechanical.hpp"
//#include "mlcpDefaultSolver.hpp"
#include "QP.hpp"
#include "Relay.hpp"
#include "FrictionContact.hpp"
#include "GlobalFrictionContact.hpp"

#include "LsodarOSI.hpp"
#include "Hem5OSI.hpp"
#include "MoreauJeanOSI.hpp"
#include "MoreauJeanBilbaoOSI.hpp"
#include "EulerMoreauOSI.hpp"
#include "NewMarkAlphaOSI.hpp"
#include "MoreauJeanDirectProjectionOSI.hpp"
#include "MoreauJeanCombinedProjectionOSI.hpp"
#include "D1MinusLinearOSI.hpp"
#include "SchatzmanPaoliOSI.hpp"
#include "ZeroOrderHoldOSI.hpp"

#include "MoreauJeanGOSI.hpp"

#include "NonSmoothEvent.hpp"
#include "TimeDiscretisationEvent.hpp"
#include "BlockCSRMatrix.hpp"
#include "MatrixIntegrator.hpp"
#include "ExtraAdditionalTerms.hpp"
