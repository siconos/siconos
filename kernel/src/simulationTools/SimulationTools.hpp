/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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

#include "Equality.hpp"
#include "LCP.hpp"
#include "OSNSMultipleImpact.hpp"
#include "MLCP.hpp"
#include "MLCPProjectOnConstraints.hpp"
#include "AVI.hpp"
//#include "MLCP2.hpp"
#include "GenericMechanical.hpp"
//#include "mlcpDefaultSolver.hpp"
#include "QP.hpp"
#include "Relay.hpp"
#include "FrictionContact.hpp"
#include "GlobalFrictionContact.hpp"

#include "LsodarOSI.hpp"
#include "Hem5OSI.hpp"
#include "MoreauJeanOSI.hpp"
#include "EulerMoreauOSI.hpp"
#include "NewMarkAlphaOSI.hpp"
#include "MoreauJeanDirectProjectionOSI.hpp"
#include "MoreauJeanCombinedProjectionOSI.hpp"
#include "MoreauJeanOSI2.hpp"
#include "D1MinusLinearOSI.hpp"
#include "SchatzmanPaoliOSI.hpp"
#include "ZeroOrderHoldOSI.hpp"

#include "MoreauJeanGOSI.hpp"

#include "NonSmoothEvent.hpp"
#include "TimeDiscretisationEvent.hpp"
#include "BlockCSRMatrix.hpp"
#include "MatrixIntegrator.hpp"
#include "ExtraAdditionalTerms.hpp"
