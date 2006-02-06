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

extern int sicLoadModelInterface(char *fname);

extern int sicInitStrategyInterface(char *fname);

extern int sicTimeGetHInterface(char *fname);

extern int sicTimeGetNInterface(char *fname);

extern int sicTimeGetKInterface(char *fname);

extern int sicSTNextStepInterface(char *fname);

extern int sicSTComputeFreeStateInterface(char *fname);

extern int sicSTcomputePbInterface(char *fname);

extern int sicSTupdateStateInterface(char *fname);

extern int sicModelgetQInterface(char *fname);

extern int sicLagrangianLinearTIDSInterface(char *fname);

extern int sicInteractionInterface(char *fname);

extern int sicLagrangianLinearRInterface(char *fname);

extern int sicNewtonImpactLawNSLInterface(char *fname);

extern int sicNonSmoothDynamicalSystemInterface(char *fname);

extern int sicModelInterface(char *fname);

extern int sicStrategyTimeSteppingInterface(char *fname);

extern int sicOneStepIntegratorMoreauInterface(char *fname);

extern int sicOneStepNSProblemLCPInterface(char *fname);

extern int sicCleanInterface(char *fname);
