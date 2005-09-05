#ifndef __SICONOSCONST__
#define __SICONOSCONST__

#include<string>

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
const std::string LINEARTIRELATION = "LinearTIR";
const std::string LAGRANGIANLINEARRELATION = "LagrangianLinearR";
const std::string LAGRANGIANNONLINEARRELATION = "LagrangianNonLinearR";
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
const std::string CFD_OSNSP = "CFD";
const std::string QP_OSNSP = "QP";
const std::string RELAY_OSNSP = "Relay";
const std::string  OSNSP_TOLERANCE = "tolerance";
const std::string  OSNSP_MAXITER = "maxIter";
const std::string  OSNSP_NORMTYPE = "normType";
const std::string  OSNSP_SEARCHDIRECTION = "searchDirection";

const std::string  OSNSP_LCPSOLVING = "LcpSolving";
const std::string  OSNSP_RPSOLVING = "RelayPrimalSolving";
const std::string  OSNSP_RDSOLVING = "RelayDualSolving";
const std::string  OSNSP_CFPSOLVING = "ContactFrictionPrimalSolving";
const std::string  OSNSP_CFDSOLVING = "ContactFrictionDualSolving";
const std::string  OSNSP_LEMKE = "Lemke";
const std::string  OSNSP_LEXICOLEMKE = "LexicoLemke";
const std::string  OSNSP_QP = "Qp" ;
const std::string  OSNSP_QPNONSYM = "Qpnonsym" ;
const std::string  OSNSP_GSNL = "Gsnl";
const std::string  OSNSP_GCP = "Gcp";
const std::string  OSNSP_LATIN = "Latin";





//const std::string NSPB_LCPSOLVING = "LcpSolving";
//const std::string NSPB_RPSOLVING = "RelayPrimalSolving";
//const std::string NSPB_RDSOLVING = "RelayDualSolving";
//const std::string NSPB_CFPSOLVING = "ContactFrictionPrimalSolving";
//const std::string NSPB_CFDSOLVING = "ContactFrictionDualSolving";
//const std::string NSPB_LEMKE = "Lemke";
//const std::string NSPB_GSNL = "Gsnl";
//const std::string NSPB_GCP = "Gcp";
//const std::string NSPB_LATIN = "Latin";
//const std::string NSPB_lemke = "lemke";
//const std::string NSPB_gsnl = "gsnl";
//const std::string NSPB_gcp = "gcp";
//const std::string NSPB_latin = "latin";

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

