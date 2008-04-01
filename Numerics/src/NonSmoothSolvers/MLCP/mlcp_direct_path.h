#ifndef MLCP_DIRECT_PATH_H
#define MLCP_DIRECT_PATH_H



/*
 * who use the direct solver:
 * First, you have to call the mlcp_direct_path_init().
 * Use a lot the mlcp_direct_path function
 * mlcp_direct_path_reset at the end.
 *
 *
 */


int mlcp_direct_path_getNbIWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
int mlcp_direct_path_getNbDWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);

void mlcp_direct_path_init(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
void mlcp_direct_path_reset();

#endif //MLCP_DIRECT_PATH_H
