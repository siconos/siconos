#ifndef NUMERICS_GMPREDUCED
#define NUMERICS_GMPREDUCED
void buildReducedGMP(GenericMechanicalProblem* pInProblem, GenericMechanicalProblem* pOutProblem);
void GMPReducedsolve(GenericMechanicalProblem* pInProblem, double *reaction , double *velocity, int* info, SolverOptions* options);
void GMPReducedSolToSol(GenericMechanicalProblem* pInProblem, double * reaction, double * velocity,
                        double * Re, double * Rreduced, double * Vreduced);
#endif
