//$Id: SiconosConst.h,v 1.5 2005/03/22 15:55:04 jbarbier Exp $
#ifndef __SICONOSCONST__
#define __SICONOSCONST__


// const found in the DynamicalSystem
const string LNLDS = "LagrangianNLDS";
const string LTIDS = "LagrangianTIDS";
const string LSDS = "LinearSystemDS";
const string NLSDS = "NonLinearSystemDS";

// const found in the BoundaryCondition
const string LINEARBC = "LinearBC";
const string NLINEARBC = "NonLinearBC";
const string PERIODICBC = "PeriodicBC";

// const found in the Interaction
const string LINEARTIRELATION = "LinearTIR";
const string LAGRANGIANLINEARRELATION = "LagrangianLinearR";
const string LAGRANGIANNONLINEARRELATION = "LagrangianNonLinearR";
const string COMPLEMENTARITYCONDITIONNSLAW = "ComplementarityConditionNSL";
const string RELAYNSLAW = "RelayNSL";
const string NEWTONIMPACTLAWNSLAW = "NewtonImpactLawNSL";
const string NEWTONIMPACTFRICTIONNSLAW = "NewtonImpactFrictionNSL";

// const for DSInputOutput
const string LINEARDSIO = "LinearDSIO";
const string NLINEARDSIO = "NLinearDSIO";
const string LAGRANGIANDSIO = "LagrangianDSIO";

// const for EqualitySonctraint
const string LINEAREC = "LinearEC";
const string NLINEAREC = "NLinearEC";
const string LINEARTIEC = "LinearTIEC";
const string LAGRANGIANEC = "LagrangianEC";

// const found in the OneStepNSProblem
const string LCP_OSNSP = "LCP";
const string QP_OSNSP = "QP";
const string RELAY_OSNSP = "Relay";
//const string NSPB_LCPSOLVING = "LcpSolving";
//const string NSPB_RPSOLVING = "RelayPrimalSolving";
//const string NSPB_RDSOLVING = "RelayDualSolving";
//const string NSPB_CFPSOLVING = "ContactFrictionPrimalSolving";
//const string NSPB_CFDSOLVING = "ContactFrictionDualSolving";
//const string NSPB_LEMKE = "Lemke";
//const string NSPB_GSNL = "Gsnl";
//const string NSPB_GCP = "Gcp";
//const string NSPB_LATIN = "Latin";
//const string NSPB_lemke = "lemke";
//const string NSPB_gsnl = "gsnl";
//const string NSPB_gcp = "gcp";
//const string NSPB_latin = "latin";

// const found in the Strategy
const string EVENTDRIVEN_STRATEGY = "EventDriven";
const string TIMESTEPPING_STRATEGY = "TimeStepping";

// const found in the OneStepIntegrator
const string MOREAU_INTEGRATOR = "Moreau";
const string ADAMS_INTEGRATOR = "Adams";
const string LSODAR_INTEGRATOR = "LSODAR";


#endif

//$Log: SiconosConst.h,v $
//Revision 1.5  2005/03/22 15:55:04  jbarbier
//- class NewtonImpactFriction non smooth law added to the kernel
//
//- xml schema modified for this new class
//- xml schema modified to accept a "joker" for further use of a LMGC90 mechanical plugin
//
//- new test added for the loading/saving of a NewtonImpactFrictionNSL
//
//Revision 1.4  2005/03/08 14:23:42  jbarbier
//- modification of constant variables :
//in the XML module, main tags of the XML objects of the strategy are in XMLTagsName.h
//
//in simualtion tools, some constants have been moved to SiconosConst.h
//
//Revision 1.3  2005/03/08 12:41:35  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.2  2005/03/07 13:17:18  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.1  2005/02/10 10:35:18  jbarbier
//- new file regrouping all the const values of the model, modelingTools and numericalStrategy
//
//- new function in the LagrangianLinearR to get the H matrix corresponding to one of the 2 dynamical systems linked to the relation
//
//- new atribute of the OneStepNSProblem. A visibility table of the Interaction.
//
