#ifndef MLCP_DIRECT_H
#define MLCP_DIRECT_H

/*
 * who use the direct solver:
 * First, you have to call the mlcp_direct_init().
 * Use a lot the mlcp_direct function.
 * add configuration with mlcp_direct_addConfigFromWSolution to add configuration.
 * mlcp_direct_reset
 *
 */


void mlcp_direct_addConfig(MixedLinearComplementarityProblem* problem, int * zw);
void mlcp_direct_addConfigFromWSolution(MixedLinearComplementarityProblem* problem, double * wSol);
void mlcp_direct_init(MixedLinearComplementarityProblem* problem, SolverOptions* options);
void mlcp_direct_reset();

int mlcp_direct_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options);
int mlcp_direct_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options);

#endif //MLCP_DIRECT_H
