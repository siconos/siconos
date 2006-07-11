#include "stack-c.h"
#include <string.h>
#include "machine.h"
#include "sun/link.h"
//#include "siconos.h"

#ifdef WINDOWS
#define extern __declspec (dllexport)
#endif


/***************************************************
 * Declarations for siconos gateway                *
 ***************************************************/
typedef int (*gate_function)(char *);

extern int C2F(SiconosGateway)();

extern int sci_gateway(char *name, gate_function f);

/***************************************************
 * Declarations for Siconos function interface     *
 * with scilab                                     *
 ***************************************************/


/* Simulation control loop */
extern int sicInitSimulationInterface(char *fname);

extern int sicTimeGetHInterface(char *fname);

extern int sicTimeGetNInterface(char *fname);

extern int sicTimeGetKInterface(char *fname);

extern int sicSTNextStepInterface(char *fname);

extern int sicSTComputeFreeStateInterface(char *fname);

extern int sicSTcomputePbInterface(char *fname);

extern int sicSTupdateStateInterface(char *fname);

extern int sicSTnewtonSolveInterface(char *fname);

/* Pick Datas */
extern int sicModelgetQInterface(char *fname);

/* DS Construction */

extern int sicLoadModelInterface(char *fname);

extern int sicLagrangianLinearTIDSInterface(char *fname);

extern int sicLagrangianDSInterface(char *fname);

extern int sicSetMassInterface(char *fname);
extern int sicSetNNLInterface(char *fname);
extern int sicSetJacQNNLInterface(char *fname);
extern int sicSetJacVelNNLInterface(char *fname);
extern int sicSetFIntInterface(char *fname);
extern int sicSetJacQFInt(char *fname);
extern int sicSetJacVelFIntInterface(char *fname);
extern int sicSetFExtInterface(char *fname);

/* Interaction Construction */

extern int sicInteractionInterface(char *fname);

extern int sicLagrangianLinearRInterface(char *fname);

extern int sicLagrangianRInterface(char *fname);

extern int sicNewtonImpactLawNSLInterface(char *fname);

extern int sicNonSmoothDynamicalSystemInterface(char *fname);


/* Main global Objects */

extern int sicModelInterface(char *fname);

/* Simulation Simulation */

extern int sicSimulationTimeSteppingInterface(char *fname);

extern int sicOneStepIntegratorMoreauInterface(char *fname);

extern int sicOneStepNSProblemLCPInterface(char *fname);

/* Clean global Objects */

extern int sicCleanInterface(char *fname);
