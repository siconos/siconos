#ifndef MLCP_TOOL
#define MLCP_TOOL

void mlcp_buildM(int * zw, double * M, double * Mref, int n, int m);
void mlcp_fillSolution(double * z1, double * z2, double * w1, double * w2, int n, int m, int* zw, double * Q);
void   mlcp_DisplaySolution(double * z1, double * z2, double * w1, double * w2, int n, int m);

#endif
