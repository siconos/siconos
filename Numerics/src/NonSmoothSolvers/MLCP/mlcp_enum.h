#ifndef MLCP_ENUM_H
#define MLCP_ENUM_H

int mlcp_enum_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options);
int mlcp_enum_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options);
/*Alloc memory iff options->iWork options->dWork and are  null.
 Return 0 if the memory is not allocated. else return 1.*/
int mlcp_enum_alloc_working_memory(MixedLinearComplementarityProblem* problem, SolverOptions* options);
/*free the Work memory, and set pointer to zero.*/
void mlcp_enum_free_working_memory(MixedLinearComplementarityProblem* problem, SolverOptions* options);
#endif //MLCP_ENUM_H
