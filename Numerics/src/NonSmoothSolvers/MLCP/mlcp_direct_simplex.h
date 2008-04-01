#ifndef MLCP_DIRECT_SIMPLEX_H
#define MLCP_DIRECT_SIMPLEX_H



/*
 * who use the direct solver:
 * First, you have to call the mlcp_direct_simplex_init().
 * Use a lot the mlcp_direct_simplex function
 * mlcp_direct_simplex_reset at the end.
 *
 *
 */


int mlcp_direct_simplex_getNbIWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
int mlcp_direct_simplex_getNbDWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);

void mlcp_direct_simplex_init(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
void mlcp_direct_simplex_reset();

#endif //MLCP_DIRECT_SIMPLEX_H
