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

void mlcp_direct_FB_init(MixedLinearComplementarityProblem* problem, SolverOptions* options);
void mlcp_direct_FB_reset();

int mlcp_direct_FB_getNbIWork(MixedLinearComplementarityProblem* problem, SolverOptions* options);
int mlcp_direct_FB_getNbDWork(MixedLinearComplementarityProblem* problem, SolverOptions* options);


#endif //MLCP_DIRECT_FB_H
