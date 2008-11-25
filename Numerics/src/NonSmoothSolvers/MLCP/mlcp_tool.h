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


void mlcp_buildM(int * zw, double * M, double * Mref, int n, int m, int NbLines);
void mlcp_fillSolution(double * z1, double * z2, double * w1, double * w2, int n, int m, int NbLines, int* zw, double * Q);
void mlcp_DisplaySolution(double * z1, double * z2, double * w1, double * w2, int n, int m, int NbLines);

#endif
