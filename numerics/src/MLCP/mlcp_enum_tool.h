// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "MixedLinearComplementarityProblem.h"

#ifndef MLCP_ENUM_TOOL_H
#define MLCP_ENUM_TOOL_H

void mlcp_enum_build_M(int* zw, double* M, double* Mref, int n, int m, int NbLines);
void mlcp_enum_build_M_Block(int* zw, double* M, double* Mref, int n, int m, int NbLines,
                             int* indexInBlock);
void mlcp_enum_fill_solution(double* z1, double* z2, double* w1, double* w2, int n, int m,
                             int NbLines, int* zw, double* Q);
void mlcp_enum_fill_solution_Block(double* z, double* w, int n, int m, int NbLines, int* zw,
                                   double* Q, int* indexInBlock);

void mlcp_enum_display_solution(double* z1, double* z2, double* w1, double* w2, int n, int m,
                                int NbLines);
void mlcp_enum_display_solution_Block(double* z, double* w, int n, int m, int Nblines,
                                      int* indexInBlock);
void mlcp_enum_build_indexInBlock(MixedLinearComplementarityProblem* problem,
                                  int* indexInBlock);

#endif  // MLCP_ENUM_H
