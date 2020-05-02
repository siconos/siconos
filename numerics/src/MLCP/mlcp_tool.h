#ifndef MLCP_TOOL
#define MLCP_TOOL
/*
  |Z1       |W1
M*|   + Q = |
  |Z2       |W2


  M size = (Nblines) X ( n+m)

  W2 and Z2 size = m
  Z1 size = n
  W1 size = NbLines - m

*/

#include "NumericsFwd.h"  // for MixedLinearComplementarityProblem
void mlcp_enum_build_M(int * zw, double * M, double * Mref, int n, int m, int NbLines);
void mlcp_enum_build_M_Block(int * zw, double * M, double * Mref, int n, int m, int NbLines, int *indexInBlock);
void mlcp_enum_fill_solution(double * z1, double * z2, double * w1, double * w2, int n, int m, int NbLines, int* zw, double * Q);
void mlcp_enum_fill_solution_Block(double * z, double * w, int n, int m, int NbLines, int* zw, double * Q, int *indexInBlock);


void mlcp_enum_display_solution(double * z1, double * z2, double * w1, double * w2, int n, int m, int NbLines);
void mlcp_enum_display_solution_Block(double * z, double * w, int n, int m, int Nblines, int *indexInBlock);
void mlcp_enum_build_indexInBlock(MixedLinearComplementarityProblem* problem, int *indexInBlock);

#endif
