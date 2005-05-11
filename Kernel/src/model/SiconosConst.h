#ifndef __SICONOSCONST__
#define __SICONOSCONST__


// const found in the DynamicalSystem
const string LNLDS = "LagrangianDS";
const string LTIDS = "LagrangianLinearTIDS";
const string LDS = "LinearDS";
const string NLDS = "NonLinearDS";

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
const string LAGRANGIANLINEARDSIO = "LagrangianLinearDSIO";

// const for EqualitySonctraint
const string LINEAREC = "LinearEC";
const string NLINEAREC = "NLinearEC";
const string LINEARTIEC = "LinearTIEC";
const string LAGRANGIANEC = "LagrangianEC";
const string LAGRANGIANLINEAREC = "LagrangianLinearEC";

// const found in the OneStepNSProblem
const string LCP_OSNSP = "LCP";
const string CFD_OSNSP = "CFD";
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

