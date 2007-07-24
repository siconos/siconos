#define SIC_OK    0
#define SIC_ERROR -1

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif

/*! \file SiconosApiC.h
\brief API C for Scilab Interface.
*/


EXTERN int sicLoadModel(char ModelXmlFile[]);

EXTERN int sicInitSimulation();

EXTERN int sicTimeGetH(double *H);

EXTERN int sicTimeGetN(int *N);

EXTERN int sicSTNextStep();

EXTERN int sicSTSaveInMemory();

EXTERN int sicSTComputeOneStep();

EXTERN int sicSTnewtonSolve(double criterion, int maxIter);

EXTERN int sicSTupdateState();

EXTERN int sicAdvanceToEvent();

EXTERN int sicProcessEvents();

EXTERN int sicHasNextEvent(int *hasnextevent);

EXTERN int sicGetTypeEvent(char *type);

EXTERN void sicDebug(int *ret);

EXTERN int sicModelgetQ(double *value, int indexDS, int indexVector);

EXTERN int sicLagrangianLinearTIDS(int nDof, double *Q0, double *Vel0, double *Mass);

EXTERN int sicLagrangianDS(int nDof, double *Q0, double *Vel0);

EXTERN int sicSetComputeMassFunction(int nIdDs, char *libname, char *func);

EXTERN int sicSetComputeNNLFunction(int nIdDs, char *libname, char *func);

EXTERN int sicSetComputeJacobianQNNLFunction(int nIdDs, char *libname, char *func);

EXTERN int sicSetComputeJacobianVelocityNNLFunction(int nIdDs, char *libname, char *func);

EXTERN int sicSetComputeFIntFunction(int nIdDs, char *libname, char *func);

EXTERN int sicSetComputeJacobianQFIntFunction(int nIdDs, char *libname, char *func);

EXTERN int sicSetComputeJacobianVelocityFIntFunction(int nIdDs, char *libname, char *func);

EXTERN int sicSetComputeFExtFunction(int nIdDs, char *libname, char *func);

EXTERN int sicSetFExt(int nIdDs, double *Fext);

EXTERN int sicInteraction(char *name, int nbDS, int *DS, int idLaw, int idRelation, int nSize);

//EXTERN int sicLagrangianLinearR(int nIdInteraction, double *H,double *b);

EXTERN int sicLagrangianLinearR(int nDof, int nRel, double *H, double *b);

// EXTERN int sicLagrangianR(int nIdInteraction, char *relationType, char *funcH, char *funcG);

EXTERN int sicNewtonImpactNSL(double e);

EXTERN int sicNonSmoothDynamicalSystem(int isBVP);

EXTERN int sicModel(double t0, double T);

EXTERN int sicTimeDiscretisation(double h);

EXTERN int sicSimulationTimeStepping(int idTime);

EXTERN int sicOneStepIntegratorMoreau(double *theta);

EXTERN int sicOneStepNSProblemLCP(char* solverName, int maxiter, double tolerance);

EXTERN int sicClean();

