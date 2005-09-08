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

// common tags
const std::string COMPUTE_INPUT_TAG = "computeInput";
const std::string COMPUTE_OUTPUT_TAG = "computeOutput";

/*
 * the different node names we can encounter
 */
// BoundaryCondition tags
const std::string LINEARBC_TAG = "Linear";
const std::string NON_LINEARBC_TAG = "NLinear";
const std::string PERIODICBC_TAG = "Periodic";

// DynamicalSystem tags
const std::string LAGRANGIAN_TIME_INVARIANTDS_TAG = "LagrangianLinearTIDS"; // LL
const std::string LAGRANGIAN_NON_LINEARDS_TAG = "LagrangianDS"; // LNL
const std::string LINEAR_SYSTEMDS_TAG = "LinearDS"; // LTI
const std::string NON_LINEAR_SYSTEMDS_TAG = "NonLinearDS";

const std::string BOUNDARYCONDITION_TAG = "BoundaryCondition";
const std::string INTERACTION_TAG = "Interaction";
const std::string RELATION_TAG = "Relation";
const std::string INTERACTION_CONTENT_TAG = "Interaction_Content";
const std::string RELATION_CONTENT_TAG = "Relation_Content";

// NSDS
const std::string DYNAMICAL_SYSTEM_TAG = "DS";
const std::string NON_SMOOTH_DYNAMICAL_SYSTEM_TAG = "NSDS";
const std::string DSINPUTOUTPUT_TAG = "DSInputOutput";
const std::string EQUALITYCONSTRAINT_TAG = "EqualityConstraint";

const std::string EQUALITYCONSTRAINT_DEFINITION_TAG = "EqualityConstraint_Definition";
const std::string DSINPUTOUTPUT_DEFINITION_TAG = "DSInputOutput_Definition";
const std::string DYNAMICAL_SYSTEM_DEFINITION_TAG = "DS_Definition";
const std::string INTERACTION_DEFINITION_TAG = "Interaction_Definition";
const std::string LMGC90_NSDS_TAG = "DS_LMGC90";
const std::string LMGC90_STRATEGY_TAG = "OneStepIntegrator_LMGC90";


// Interaction
const std::string LAGRANGIAN_LINEAR_RELATION_TAG = "LagrangianLinear";
const std::string LAGRANGIAN_NON_LINEAR_RELATION_TAG = "LagrangianNonLinear";
const std::string LINEAR_TIME_INVARIANT_RELATION_TAG = "LinearTimeInvariant";
const std::string NON_LINEAR_RELATION_TAG = "NonLinear";

const std::string COMPLEMENTARITY_CONDITION_NSLAW_TAG = "ComplementarityCondition";
const std::string RELAY_NSLAW_TAG = "Relay";
const std::string NEWTON_IMPACT_LAW_NSLAW_TAG = "NewtonImpactLaw";
const std::string NEWTON_IMPACT_FRICTION_NSLAW_TAG = "NewtonImpactFrictionLaw";


// DSIO
const std::string LINEAR_DSIO_TAG = "LinearDSIO";
const std::string NON_LINEAR_DSIO_TAG = "NonLinearDSIO";
const std::string LAGRANGIAN_DSIO_TAG = "LagrangianDSIO";
const std::string LAGRANGIAN_LINEAR_DSIO_TAG = "LagrangianLinearDSIO";

// EqualityConstraint
const std::string LINEAR_EC_TAG = "LinearEC";
const std::string NON_LINEAR_EC_TAG = "NonLinearEC";
const std::string LINEAR_TIME_INVARIANT_EC_TAG = "LinearTIEC";
const std::string LAGRANGIAN_EC_TAG = "LagrangianEC";
const std::string LAGRANGIAN_LINEAR_EC_TAG = "LagrangianLinearEC";


const std::string DSIO_CONCERNED = "DSInputOutput_Concerned";
const std::string DS_CONCERNED = "DS_Concerned";

//===========================================

// SiconosModel
const std::string MODEL_TAG = "SiconosModel";
const std::string NSDS_TAG = "NSDS";
const std::string STRATEGY_TAG = "Strategy";

// Strategy
const std::string ONESTEPINTEGRATOR_DEFINITION_TAG = "OneStepIntegrator_Definition";
const std::string ONESTEPINTEGRATOR_TAG = "OneStepIntegrator";
const std::string ONESTEPNSPROBLEM_TAG = "OneStepNSProblem";

const std::string TIMEDISCRETISATION_TAG = "TimeDiscretisation";
const std::string TIMESTEPPING_TAG = "TimeStepping";
const std::string EVENTDRIVEN_TAG = "EventDriven";

/*Types of OneStepIntegrator defined*/
const std::string MOREAU_TAG = "Moreau";
const std::string LSODAR_TAG = "LSODAR";
const std::string ADAMS_TAG = "Adams";

/*Types of OneStepNSProblem defined*/
const std::string LCP_TAG = "LCP";
const std::string CFD_TAG = "CFD";
const std::string QP_TAG = "QP";
const std::string RELAY_TAG = "Relay";

//===========================================
//===========================================

//== BoundaryConditionXML ==
//const std::string BC_TYPE = "type";
//const std::string BC_NLINEAR = "NLinear";
//const std::string BC_LINEAR = "Linear";
//const std::string BC_PERIODIC = "Periodic";
//==========================


//== DSInputOutputXML ==
//const std::string DSIO_TYPE = "type";
//const std::string DSIO_INPUT = "computeInput";
//const std::string DSIO_OUTPUT = "computeOutput";
//const std::string DSIO_PLUGIN = "plugin";
//
//const std::string DSIO_LL = "LL";
//const std::string DSIO_LNL = "LNL";
//const std::string DSIO_LTI = "LTI";
//==========================


////== InteractionXML ==
////const std::string INTERACTION_NODE = "Interaction";
////const std::string INTERACTION_CONTENT = "Interaction_Content";
////const std::string INTERACTION_ID = "Id";
//const std::string INTERACTION_Y = "y";
//const std::string INTERACTION_LAMBDA = "lambda";
//const std::string INTERACTION_NINTER = "nInter";
////const std::string INTERACTION_ISACTIVE = "isActive";
//
//const std::string INTERACTION_DS_CONCERNED = "DS_Concerned";
////const std::string INTERACTION_NS_LAW = "NS_Law";
////const std::string INTERACTION_RELATION = "Relation";
////const std::string INTERACTION_DS = "DS";
//const std::string INTERACTION_INTERACTWITHDS_NUMBER = "interactsWithDS_Number";
////const std::string INTERACTION_NUMBER="number";
//
////const std::string INTERACTION_TYPE="type";
////const std::string INTERACTION_SIZE="size";
////const std::string INTERACTION_ALL="all";
//
////const std::string INTERACTION_LL = "LL";
////const std::string INTERACTION_LNL = "LNL";
////const std::string INTERACTION_LTI = "LTI";
////
////const std::string INTERACTION_COMPLEMENTARITYCONDITIONNSLAW = "ComplementarityCondition";
////const std::string INTERACTION_RELAYNSLAW = "Relay";
////const std::string INTERACTION_NEWTONIMPACTLAWNSLAW = "NewtonImpactLaw";
//==========================


////== DSXML ==
////const std::string DS_ID = "Id";
//const std::string DS_N = "n";
//const std::string DS_X0 = "x0";
//const std::string DS_X = "x";
//const std::string DS_XDOT = "xDot";
//const std::string DS_R = "R";
//const std::string DS_XMEMORY = "xMemory";
//const std::string DS_XDOTMEMORY = "xDotMemory";
//const std::string DS_RMEMORY = "rMemory";
//const std::string DS_STEPSINMEMORY = "StepsInMemory";
////const std::string DS_BOUNDARYCONDITION = "BoundaryCondition";
//const std::string DS_VECTORFIELD = "vectorField";
//const std::string DS_COMPUTEJACOBIANX = "computeJacobianX";
//
////const std::string DS_PLUGIN = "plugin";

//attributes
//const std::string DS_NUMBER = "number";
//const std::string DS_TYPE = "type";

//const std::string DS_NLINEAR = "NLinear";
//const std::string DS_LINEAR = "Linear";
//const std::string DS_PERIODIC = "Periodic";
//==========================


////== LagrangianLinearRXML ==
//const std::string LLR_H = "H";
//const std::string LLR_B = "b";
////==========================
//
////== LagrangianNonLinearDS ==
//const std::string LNLDS_Q = "q";
//const std::string LNLDS_Q0 = "q0";
//const std::string LNLDS_QMEMORY = "qMemory";
//
//const std::string LNLDS_VELOCITY = "Velocity";
//const std::string LNLDS_VELOCITY0 = "Velocity0";
//const std::string LNLDS_VELOCITYMEMORY = "VelocityMemory";
//
//const std::string LNLDS_QNLINERTIA = "QNLInertia";
//const std::string LNLDS_FINT = "Fint";
//const std::string LNLDS_FEXT = "Fext";
//
//const std::string LNLDS_JACOBIANQFINT = "JacobianQFint";
//const std::string LNLDS_JACOBIANVELOCITYFINT = "JacobianVelocityFint";
//const std::string LNLDS_JACOBIANQQNLINERTIA = "JacobianQQNLInertia";
//const std::string LNLDS_JACOBIANVELOCITYQNLINERTIA = "JacobianVelocityQNLInertia";
//
//const std::string LNLDS_M = "M";
//const std::string LNLDS_NDOF = "ndof";
//const std::string LNLDS_MATRIXPLUGIN = "matrixPlugin";
//const std::string LNLDS_VECTORPLUGIN = "vectorPlugin";
//==========================


//== LagrangianNonLinearDS ==
//==========================



// NSDS
//const std::string NSDS_NODE = "NSDS";
//const std::string NSDS_EQUALITYCONSTRAINT_DEFINITION = "EqualityConstraint_Definition";
//const std::string NSDS_DSINPUTOUTPUT_DEFINITION = "DSInputOutput_Definition";
//const std::string NSDS_DS_DEFINITION = "DS_Definition";
//const std::string NSDS_INTERACTION_DEFINITION = "Interaction_Definition";

//const std::string NSDS_INTERACTION = "Interaction";

//const std::string NSDS_DS = "DS";
//const std::string NSDS_NUMBER = "number";
//const std::string NSDS_TYPE = "type";
//const std::string NSDS_BVP = "bvp";


/*Types of DS defined*/
//const std::string NSDS_LAGRANGIANNLDS = "LagrangianDS";
//const std::string NSDS_LAGRANGIANTIDS = "LagrangianLinearTIDS";
//const std::string NSDS_LINEARSYSTEM = "LinearDS";
//const std::string NSDS_NONLINEARSYSTEM = "NonLinearDS";



#endif

