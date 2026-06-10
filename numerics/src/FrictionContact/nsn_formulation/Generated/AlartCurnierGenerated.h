/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef AlartCurnierGenerated_h
#define AlartCurnierGenerated_h

#if defined(__cplusplus)
extern "C" {
#endif

void fc3d_AlartCurnierFunctionGenerated(double *reaction, double *velocity, double mu,
                                        double *rho, double *f, double *A, double *B);

void fc3d_AlartCurnierJeanMoreauFunctionGenerated(double *reaction, double *velocity,
                                                  double mu, double *rho, double *f, double *A,
                                                  double *B);

#if defined(__cplusplus)
}
#endif
#endif
