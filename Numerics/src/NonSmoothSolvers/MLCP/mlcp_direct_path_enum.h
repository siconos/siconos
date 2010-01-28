#ifndef MLCP_DIRECT_PATH_ENUM_H
#define MLCP_DIRECT_PATH_ENUM_H


int mlcp_direct_path_enum_getNbIWork(MixedLinearComplementarity_Problem* problem, SolverOptions* options);
int mlcp_direct_path_enum_getNbDWork(MixedLinearComplementarity_Problem* problem, SolverOptions* options);

void mlcp_direct_path_enum(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, SolverOptions* options);
void mlcp_direct_path_enum_reset();
void mlcp_direct_path_enum_init(MixedLinearComplementarity_Problem* problem, SolverOptions* options);

#endif //MLCP_DIRECT_PATH_ENUM_H
