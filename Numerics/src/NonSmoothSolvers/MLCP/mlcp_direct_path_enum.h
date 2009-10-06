#ifndef MLCP_DIRECT_PATH_ENUM_H
#define MLCP_DIRECT_PATH_ENUM_H

int mlcp_direct_path_enum_getNbIWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
int mlcp_direct_path_enum_getNbDWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);

void mlcp_direct_path_enum(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options);
void mlcp_direct_path_enum_reset();
void mlcp_direct_path_enum_init(MixedLinearComplementarity_Problem* problem, Solver_Options* options);

#endif //MLCP_DIRECT_PATH_ENUM_H
