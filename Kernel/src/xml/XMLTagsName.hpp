/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#ifndef __XMLTAGSNAME__
#define __XMLTAGSNAME__

#include<string>

/*
 * the different kind of attribute we can encounter
 */
const std::string TYPE_ATTRIBUTE = "type";
const std::string ID_ATTRIBUTE = "Id";
const std::string NUMBER_ATTRIBUTE = "number";
const std::string PLUGIN_ATTRIBUTE = "plugin";
const std::string ALL_ATTRIBUTE = "all";
const std::string SIZE_ATTRIBUTE = "size";
// usefull when a list of objects, identified thanks to numbers, is required;
// for example a list of DS in an interaction:
const std::string INDEX_LIST = "indexList";
//Tags
const std::string STEPSINMEMORY = "StepsInMemory";
const std::string MATRIXPLUGIN = "matrixPlugin";
const std::string VECTORPLUGIN = "vectorPlugin";

/*
 * the different node names we can encounter
 */

// DynamicalSystem tags
const std::string LAGRANGIAN_TIDS_TAG = "LagrangianLinearTIDS";
const std::string LAGRANGIAN_NON_LINEARDS_TAG = "LagrangianDS";
const std::string LINEAR_DS_TAG = "FirstOrderLinearDS";
const std::string LINEAR_TIDS_TAG = "FirstOrderLinearTIDS";
const std::string NON_LINEAR_DS_TAG = "FirstOrderNonLinearDS";

const std::string INTERACTION_TAG = "Interaction";
const std::string INTERACTION_CONTENT_TAG = "Interaction_Content";
const std::string RELATION_CONTENT_TAG = "Relation_Content";

// NSDS
const std::string DYNAMICAL_SYSTEM_TAG = "DS";
const std::string NON_SMOOTH_DYNAMICAL_SYSTEM_TAG = "NSDS";

const std::string DYNAMICAL_SYSTEM_DEFINITION_TAG = "DS_Definition";
const std::string INTERACTION_DEFINITION_TAG = "Interaction_Definition";
const std::string LMGC90_NSDS_TAG = "DS_LMGC90";
const std::string LMGC90_SIMULATION_TAG = "OneStepIntegrator_LMGC90";


// Interaction
// - Relations -
const std::string RELATION_TAG = "Relation";
const std::string LINEAR_TIME_INVARIANT_RELATION_TAG = "LinearTimeInvariantRelation";
const std::string LAGRANGIAN_RELATION_TAG = "LagrangianRelation";
// - Non-smooth laws -
const std::string COMPLEMENTARITY_CONDITION_NSLAW_TAG = "ComplementarityCondition";
const std::string RELAY_NSLAW_TAG = "Relay";
const std::string NEWTON_IMPACT_NSLAW_TAG = "NewtonImpactLaw";
const std::string NEWTON_IMPACT_FRICTION_NSLAW_TAG = "NewtonImpactFrictionLaw";

const std::string DS_CONCERNED = "DS_Concerned";

//===========================================

// SiconosModel
const std::string MODEL_TAG = "SiconosModel";
const std::string NSDS_TAG = "NSDS";
const std::string SIMULATION_TAG = "Simulation";

// Simulation
const std::string ONESTEPINTEGRATOR_DEFINITION_TAG = "OneStepIntegrator_Definition";
const std::string ONESTEPINTEGRATOR_TAG = "OneStepIntegrator";
const std::string ONESTEPNSPROBLEM_TAG = "OneStepNSProblems_List";

const std::string TIMEDISCRETISATION_TAG = "TimeDiscretisation";
const std::string TIMESTEPPING_TAG = "TimeStepping";
const std::string EVENTDRIVEN_TAG = "EventDriven";

/*Types of OneStepIntegrator defined*/
const std::string MOREAU_TAG = "Moreau";
const std::string LSODAR_TAG = "Lsodar";
const std::string ADAMS_TAG = "Adams";

/*Types of OneStepNSProblem defined*/
const std::string LCP_TAG = "LCP";
const std::string FRICTIONCONTACT_TAG = "FrictionContact";
const std::string QP_TAG = "QP";
const std::string RELAY_TAG = "Relay";


#endif

