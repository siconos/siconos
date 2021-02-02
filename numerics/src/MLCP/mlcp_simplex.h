#ifndef MLCP_SIMPLEX_H
#define MLCP_SIMPLEX_H

/*
 *
 */


#include "NumericsFwd.h"  // for MixedLinearComplementarityProblem, SolverOp...
void mlcp_simplex_init(MixedLinearComplementarityProblem* problem, SolverOptions* options);
void mlcp_simplex_reset(void);

#endif //MLCP_SIMPLEX_H
