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


void mlcp_direct_addConfig(MixedLinearComplementarity_Problem* problem, int * zw);
void mlcp_direct_addConfigFromWSolution(MixedLinearComplementarity_Problem* problem, double * wSol);
void mlcp_direct_init(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
void mlcp_direct_reset();

int mlcp_direct_getNbIWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
int mlcp_direct_getNbDWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);

#endif //MLCP_DIRECT_H
