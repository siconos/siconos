#ifndef __XMLTAGSNAME__
#define __XMLTAGSNAME__

/*
 * the different kind of attribute we can encounter
 */
const string TYPE_ATTRIBUTE = "type";
const string ID_ATTRIBUTE = "Id";
const string NUMBER_ATTRIBUTE = "number";
const string PLUGIN_ATTRIBUTE = "plugin";
const string ALL_ATTRIBUTE = "all";
const string SIZE_ATTRIBUTE = "size";

// common tags
const string COMPUTE_INPUT_TAG = "computeInput";
const string COMPUTE_OUTPUT_TAG = "computeOutput";



/*
 * the different node names we can encounter
 */
// BoundaryCondition tags
const string LINEARBC_TAG = "Linear";
const string NON_LINEARBC_TAG = "NLinear";
const string PERIODICBC_TAG = "Periodic";

// DynamicalSystem tags
const string LAGRANGIAN_TIME_INVARIANTDS_TAG = "LagrangianLinearTIDS"; // LL
const string LAGRANGIAN_NON_LINEARDS_TAG = "LagrangianDS"; // LNL
const string LINEAR_SYSTEMDS_TAG = "LinearSystemDS"; // LTI
const string NON_LINEAR_SYSTEMDS_TAG = "NonLinearSystemDS";

const string BOUNDARYCONDITION_TAG = "BoundaryCondition";
const string INTERACTION_TAG = "Interaction";
const string RELATION_TAG = "Relation";
const string INTERACTION_CONTENT_TAG = "Interaction_Content";
const string RELATION_CONTENT_TAG = "Relation_Content";

// NSDS
const string DYNAMICAL_SYSTEM_TAG = "DS";
const string NON_SMOOTH_DYNAMICAL_SYSTEM_TAG = "NSDS";
const string DSINPUTOUTPUT_TAG = "DSInputOutput";
const string EQUALITYCONSTRAINT_TAG = "EqualityConstraint";

const string EQUALITYCONSTRAINT_DEFINITION_TAG = "EqualityConstraint_Definition";
const string DSINPUTOUTPUT_DEFINITION_TAG = "DSInputOutput_Definition";
const string DYNAMICAL_SYSTEM_DEFINITION_TAG = "DS_Definition";
const string INTERACTION_DEFINITION_TAG = "Interaction_Definition";
const string LMGC90_NSDS_TAG = "DS_LMGC90";
const string LMGC90_STRATEGY_TAG = "OneStepIntegrator_LMGC90";


// Interaction
const string LAGRANGIAN_LINEAR_RELATION_TAG = "LL";
const string LAGRANGIAN_NON_LINEAR_RELATION_TAG = "LNL";
const string LINEAR_TIME_INVARIANT_RELATION_TAG = "LTI";

const string COMPLEMENTARITY_CONDITION_NSLAW_TAG = "ComplementarityCondition";
const string RELAY_NSLAW_TAG = "Relay";
const string NEWTON_IMPACT_LAW_NSLAW_TAG = "NewtonImpactLaw";
const string NEWTON_IMPACT_FRICTION_NSLAW_TAG = "NewtonImpactFrictionLaw";


// DSIO
const string LINEAR_DSIO_TAG = "LinearDSIO";
const string NON_LINEAR_DSIO_TAG = "NonLinearDSIO";
const string LAGRANGIAN_DSIO_TAG = "LagrangianDSIO";

// EqualityConstraint
const string LINEAR_EC_TAG = "LinearEC";
const string NON_LINEAR_EC_TAG = "NonLinearEC";
const string LINEAR_TIME_INVARIANT_EC_TAG = "LinearTIEC";
const string LAGRANGIAN_EC_TAG = "LagrangianEC";


const string DSIO_CONCERNED = "DSInputOutput_Concerned";
const string DS_CONCERNED = "DS_Concerned";

//===========================================

// SiconosModel
const string MODEL_TAG = "SiconosModel";
const string NSDS_TAG = "NSDS";
const string STRATEGY_TAG = "Strategy";

// Strategy
const string ONESTEPINTEGRATOR_DEFINITION_TAG = "OneStepIntegrator_Definition";
const string ONESTEPINTEGRATOR_TAG = "OneStepIntegrator";
const string ONESTEPNSPROBLEM_TAG = "OneStepNSProblem";

const string TIMEDISCRETISATION_TAG = "TimeDiscretisation";
const string TIMESTEPPING_TAG = "TimeStepping";
const string EVENTDRIVEN_TAG = "EventDriven";

/*Types of OneStepIntegrator defined*/
const string MOREAU_TAG = "Moreau";
const string LSODAR_TAG = "LSODAR";
const string ADAMS_TAG = "Adams";

/*Types of OneStepNSProblem defined*/
const string LCP_TAG = "LCP";
const string QP_TAG = "QP";
const string RELAY_TAG = "Relay";

//===========================================
//===========================================

//== BoundaryConditionXML ==
//const string BC_TYPE = "type";
//const string BC_NLINEAR = "NLinear";
//const string BC_LINEAR = "Linear";
//const string BC_PERIODIC = "Periodic";
//==========================


//== DSInputOutputXML ==
//const string DSIO_TYPE = "type";
//const string DSIO_INPUT = "computeInput";
//const string DSIO_OUTPUT = "computeOutput";
//const string DSIO_PLUGIN = "plugin";
//
//const string DSIO_LL = "LL";
//const string DSIO_LNL = "LNL";
//const string DSIO_LTI = "LTI";
//==========================


////== InteractionXML ==
////const string INTERACTION_NODE = "Interaction";
////const string INTERACTION_CONTENT = "Interaction_Content";
////const string INTERACTION_ID = "Id";
//const string INTERACTION_STATUS = "Status";
//const string INTERACTION_Y = "y";
//const string INTERACTION_LAMBDA = "lambda";
//const string INTERACTION_NINTER = "nInter";
////const string INTERACTION_ISACTIVE = "isActive";
//
//const string INTERACTION_DS_CONCERNED = "DS_Concerned";
////const string INTERACTION_NS_LAW = "NS_Law";
////const string INTERACTION_RELATION = "Relation";
////const string INTERACTION_DS = "DS";
//const string INTERACTION_INTERACTWITHDS_NUMBER = "interactsWithDS_Number";
////const string INTERACTION_NUMBER="number";
//
////const string INTERACTION_TYPE="type";
////const string INTERACTION_SIZE="size";
////const string INTERACTION_ALL="all";
//
////const string INTERACTION_LL = "LL";
////const string INTERACTION_LNL = "LNL";
////const string INTERACTION_LTI = "LTI";
////
////const string INTERACTION_COMPLEMENTARITYCONDITIONNSLAW = "ComplementarityCondition";
////const string INTERACTION_RELAYNSLAW = "Relay";
////const string INTERACTION_NEWTONIMPACTLAWNSLAW = "NewtonImpactLaw";
//==========================


////== DSXML ==
////const string DS_ID = "Id";
//const string DS_N = "n";
//const string DS_X0 = "x0";
//const string DS_X = "x";
//const string DS_XDOT = "xDot";
//const string DS_R = "R";
//const string DS_XMEMORY = "xMemory";
//const string DS_XDOTMEMORY = "xDotMemory";
//const string DS_RMEMORY = "rMemory";
//const string DS_STEPSINMEMORY = "StepsInMemory";
////const string DS_BOUNDARYCONDITION = "BoundaryCondition";
//const string DS_VECTORFIELD = "vectorField";
//const string DS_COMPUTEJACOBIANX = "computeJacobianX";
//
////const string DS_PLUGIN = "plugin";

//attributes
//const string DS_NUMBER = "number";
//const string DS_TYPE = "type";

//const string DS_NLINEAR = "NLinear";
//const string DS_LINEAR = "Linear";
//const string DS_PERIODIC = "Periodic";
//==========================


////== LagrangianLinearRXML ==
//const string LLR_H = "H";
//const string LLR_B = "b";
////==========================
//
////== LagrangianNonLinearDS ==
//const string LNLDS_Q = "q";
//const string LNLDS_Q0 = "q0";
//const string LNLDS_QMEMORY = "qMemory";
//
//const string LNLDS_VELOCITY = "Velocity";
//const string LNLDS_VELOCITY0 = "Velocity0";
//const string LNLDS_VELOCITYMEMORY = "VelocityMemory";
//
//const string LNLDS_QNLINERTIA = "QNLInertia";
//const string LNLDS_FINT = "Fint";
//const string LNLDS_FEXT = "Fext";
//
//const string LNLDS_JACOBIANQFINT = "JacobianQFint";
//const string LNLDS_JACOBIANVELOCITYFINT = "JacobianVelocityFint";
//const string LNLDS_JACOBIANQQNLINERTIA = "JacobianQQNLInertia";
//const string LNLDS_JACOBIANVELOCITYQNLINERTIA = "JacobianVelocityQNLInertia";
//
//const string LNLDS_M = "M";
//const string LNLDS_NDOF = "ndof";
//const string LNLDS_MATRIXPLUGIN = "matrixPlugin";
//const string LNLDS_VECTORPLUGIN = "vectorPlugin";
//==========================


//== LagrangianNonLinearDS ==
//==========================



// NSDS
//const string NSDS_NODE = "NSDS";
//const string NSDS_EQUALITYCONSTRAINT_DEFINITION = "EqualityConstraint_Definition";
//const string NSDS_DSINPUTOUTPUT_DEFINITION = "DSInputOutput_Definition";
//const string NSDS_DS_DEFINITION = "DS_Definition";
//const string NSDS_INTERACTION_DEFINITION = "Interaction_Definition";

//const string NSDS_INTERACTION = "Interaction";
//const string NSDS_STATUS = "status";

//const string NSDS_DS = "DS";
//const string NSDS_NUMBER = "number";
//const string NSDS_TYPE = "type";
//const string NSDS_BVP = "bvp";


/*Types of DS defined*/
//const string NSDS_LAGRANGIANNLDS = "LagrangianDS";
//const string NSDS_LAGRANGIANTIDS = "LagrangianLinearTIDS";
//const string NSDS_LINEARSYSTEM = "LinearSystemDS";
//const string NSDS_NONLINEARSYSTEM = "NonLinearSystemDS";



#endif

