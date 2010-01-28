#ifndef MLCP_DIRECT_PATH_ENUM_H
#define MLCP_DIRECT_PATH_ENUM_H


int mlcp_direct_path_enum_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options);
int mlcp_direct_path_enum_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options);

void mlcp_direct_path_enum(MixedLinearComplementarityProblem* problem, double *z, double *w, int *info, SolverOptions* options);
void mlcp_direct_path_enum_reset();
void mlcp_direct_path_enum_init(MixedLinearComplementarityProblem* problem, SolverOptions* options);

#endif //MLCP_DIRECT_PATH_ENUM_H
