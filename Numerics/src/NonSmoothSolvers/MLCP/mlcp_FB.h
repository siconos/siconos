#ifndef MLCP_FB_H
#define MLCP_FB_H

void mlcp_FB_init(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
void mlcp_FB_reset();


int mlcp_FB_getNbIWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);
int mlcp_FB_getNbDWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options);

#endif //MLCP_FB_H
