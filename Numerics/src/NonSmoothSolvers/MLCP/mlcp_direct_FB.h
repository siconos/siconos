#ifndef MLCP_DIRECT_FB_H
#define MLCP_DIRECT_FB_H



/*
 * who use the direct solver:
 * First, you have to call the mlcp_direct_path_init().
 * Use a lot the mlcp_direct_path function
 * mlcp_direct_FB_reset at the end.
 *
 *
 */
void mlcp_direct_FB_init(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
void mlcp_direct_FB_reset();

int mlcp_direct_FB_getNbIWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
int mlcp_direct_FB_getNbDWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);


#endif //MLCP_DIRECT_FB_H
