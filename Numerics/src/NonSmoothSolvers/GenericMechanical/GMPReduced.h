#ifndef NUMERICS_GMPREDUCED
#define NUMERICS_GMPREDUCED
void buildReducedGMP(GenericMechanicalProblem* pInProblem, GenericMechanicalProblem* pOutProblem);
void GMPReducedsolve(GenericMechanicalProblem* pInProblem, double *reaction , double *velocity, SolverOptions* options);
#endif
