/* Siconos version 1.0, Copyright INRIA 2005.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#ifndef __SICONOSCONST__
#define __SICONOSCONST__

#include<string>

const std::string N_ASCII = "ascii";
const std::string N_BINARY = "binary";

// const found in the DynamicalSystem
const std::string LNLDS = "LagrangianDS";
const std::string LTIDS = "LagrangianLinearTIDS";
const std::string LDS = "LinearDS";
const std::string NLDS = "NonLinearDS";

// const found in the BoundaryCondition
const std::string LINEARBC = "LinearBC";
const std::string NLINEARBC = "NonLinearBC";
const std::string PERIODICBC = "PeriodicBC";

// const found in the Interaction
// Relations:
const std::string RELATION = "Relation";
const std::string LINEARTIRELATION = "LinearTIR";
const std::string LAGRANGIANRELATION = "LagrangianR";
const std::string LAGRANGIANLINEARRELATION = "LagrangianLinearR";
// Non smooth laws
const std::string COMPLEMENTARITYCONDITIONNSLAW = "ComplementarityConditionNSL";
const std::string RELAYNSLAW = "RelayNSL";
const std::string NEWTONIMPACTLAWNSLAW = "NewtonImpactLawNSL";
const std::string NEWTONIMPACTFRICTIONNSLAW = "NewtonImpactFrictionNSL";

// const for DSInputOutput
const std::string LINEARDSIO = "LinearDSIO";
const std::string NLINEARDSIO = "NLinearDSIO";
const std::string LAGRANGIANDSIO = "LagrangianDSIO";
const std::string LAGRANGIANLINEARDSIO = "LagrangianLinearDSIO";

// const for EqualityConstraint
const std::string LINEAREC = "LinearEC";
const std::string NLINEAREC = "NLinearEC";
const std::string LINEARTIEC = "LinearTIEC";
const std::string LAGRANGIANEC = "LagrangianEC";
const std::string LAGRANGIANLINEAREC = "LagrangianLinearEC";

// const found in the OneStepNSProblem
const std::string LCP_OSNSP = "LCP";
const std::string DFC_2D_OSNSP = "DFC_2D";
const std::string QP_OSNSP = "QP";
const std::string RELAY_OSNSP = "Relay";
const std::string  OSNSP_TOLERANCE = "tolerance";
const std::string  OSNSP_MAXITER = "maxIter";
const std::string  OSNSP_NORMTYPE = "normType";
const std::string  OSNSP_SEARCHDIRECTION = "searchDirection";

const std::string  OSNSP_LCPSOLVING = "LcpSolving";
const std::string  OSNSP_PRSOLVING = "PrimalRelaySolving";
const std::string  OSNSP_DRSOLVING = "DualRelaySolving";
const std::string  OSNSP_PFC_2DSOLVING = "PrimalFrictionContact2DSolving";
const std::string  OSNSP_DFC_2DSOLVING = "DualFrictionContact2DSolving";
const std::string  OSNSP_LEMKE = "Lemke";
const std::string  OSNSP_LEXICOLEMKE = "LexicoLemke";
const std::string  OSNSP_QP = "QP" ;
const std::string  OSNSP_NSQP = "NSQP" ;
const std::string  OSNSP_NLGS = "NLGS";
const std::string  OSNSP_CPG = "CPG";
const std::string  OSNSP_LATIN = "Latin";

// const found in the Strategy
const std::string EVENTDRIVEN_STRATEGY = "EventDriven";
const std::string TIMESTEPPING_STRATEGY = "TimeStepping";

// const found in the OneStepIntegrator
const std::string MOREAU_INTEGRATOR = "Moreau";
const std::string ADAMS_INTEGRATOR = "Adams";
const std::string LSODAR_INTEGRATOR = "LSODAR";

// To initialize pointers:
#ifndef NULL
const int NULL = 0;
#endif


#endif

