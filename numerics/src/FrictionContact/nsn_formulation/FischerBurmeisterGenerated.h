// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef FischerBurmeisterGenerated_h
#define FischerBurmeisterGenerated_h

#if defined(__cplusplus)
extern "C" {
#endif

int fc3d_FischerBurmeisterFABGenerated(double rn, double rt1, double rt2, double un,
                                       double ut1, double ut2, double mu, double rhon,
                                       double rhot1, double rhot2, double* result);

int fc3d_FischerBurmeisterFGenerated(double rn, double rt1, double rt2, double un, double ut1,
                                     double ut2, double mu, double rhon, double rhot1,
                                     double rhot2, double* result);

int fc3d_FischerBurmeisterABGenerated(double rn, double rt1, double rt2, double un, double ut1,
                                      double ut2, double mu, double rhon, double rhot1,
                                      double rhot2, double* result);

int fc3d_FischerBurmeisterFunctionGenerated(double* reaction, double* velocity, double mu,
                                            double* rho, double* f, double* A, double* B);

int fc3d_FischerBurmeisterFMeritGenerated(double rn, double rt1, double rt2, double un,
                                          double ut1, double ut2, double mu, double rhon,
                                          double rhot1, double rhot2, double* result);

int fc3d_FischerBurmeisterGradFMeritGenerated(double rn, double rt1, double rt2, double un,
                                              double ut1, double ut2, double mu, double rhon,
                                              double rhot1, double rhot2, double* result);

int fc3d_FischerBurmeisterGradMeritFunctionGenerated(double* reaction, double* velocity,
                                                     double mu, double* rho, double* gf,
                                                     double*, double*);

#if defined(__cplusplus)
}
#endif
#endif
