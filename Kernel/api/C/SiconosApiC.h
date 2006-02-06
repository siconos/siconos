#define SIC_OK    0
#define SIC_ERROR -1

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicLoadModel(char ModelXmlFile[]);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicInitStrategy();

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicTimeGetH(double *H);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicTimeGetN(int *N);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicTimeGetK(int *K);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicSTNextStep();

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicSTComputeFreeState();

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicSTComputePb();

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicSTupdateState();

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicModelgetQ(double *value, int indexDS, int indexVector);


extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
void sicDebug(int *ret);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicLagrangianLinearTIDS(int nDof, double *Q0, double *Vel0, double *Mass, double *K, double *C, char *libname, char * fctname);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicLagrangianDS(int nDof, double *Q0, double *Vel0, char *libname,
                    char *fctMass, char *fctNNL, char *fctJacNNL,
                    char *fctFInt, char *fctJacFInt, char *fctFext);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicInteraction(char *name, int nbDS, int *DS, int nbRel);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicLagrangianLinearR(int nIdInteraction, double *H, double *b);


extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicNewtonImpactLawNSL(int nIdInteraction, char *lawtype, double e);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicNonSmoothDynamicalSystem(int isBVP);


extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicModel(double t0, double T);


extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicStrategyTimeStepping(double h);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicOneStepIntegratorMoreau(double *theta);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicOneStepNSProblemLCP(double maxiter, double tolerance);

extern
#ifdef __cplusplus
"C"
#endif /* __cplusplus */
int sicClean();

