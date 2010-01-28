#ifndef MLCP_PATH_ENUM_H
#define MLCP_PATH_ENUM_H



/*
 * who use the direct solver:
 * First, you have to call the mlcp_direct_enum_init().
 * Use a lot the mlcp_direct_enum function
 * mlcp_direct_enum_reset at the end.
 *
 *
 */


int mlcp_path_enum_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options);
int mlcp_path_enum_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options);

void mlcp_path_enum_init(MixedLinearComplementarityProblem* problem, SolverOptions* options);
void mlcp_path_enum_reset();
void mlcp_path_enum(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);
#endif //MLCP_PATH_ENUM_H
