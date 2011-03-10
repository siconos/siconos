#ifndef NUMERICS_GMPREDUCED
#define NUMERICS_GMPREDUCED
void buildReducedGMP(GenericMechanicalProblem* pInProblem, GenericMechanicalProblem* pOutProblem);
void GMPReducedSolve(GenericMechanicalProblem* pInProblem, double *reaction , double *velocity, int* info, SolverOptions* options);
void GMPReducedEqualitySolve(GenericMechanicalProblem* pInProblem, double *reaction , double *velocity, int* info, SolverOptions* options);
void GMPReducedSolToSol(GenericMechanicalProblem* pInProblem, double * reaction, double * velocity,
                        double * Re, double * Rreduced, double * Vreduced);
void printDenseMatrice(char* name, FILE * titi, double * m, int N, int M);
#endif
