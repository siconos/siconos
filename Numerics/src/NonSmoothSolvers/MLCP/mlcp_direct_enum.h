#ifndef MLCP_DIRECT_ENUM_H
#define MLCP_DIRECT_ENUM_H



/*
 * who use the direct solver:
 * First, you have to call the mlcp_direct_enum_init().
 * Use a lot the mlcp_direct_enum function
 * mlcp_direct_enum_reset at the end.
 *
 *
 */


int mlcp_direct_enum_getNbIWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
int mlcp_direct_enum_getNbDWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);

void mlcp_direct_enum_init(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
void mlcp_direct_enum_reset();

#endif //MLCP_DIRECT_ENUM_H
