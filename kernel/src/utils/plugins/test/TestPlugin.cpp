/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#include <cstdio>

#if defined(_MSC_VER)
#define DLLEXPORT __declspec(dllexport)
#else
#define DLLEXPORT
#endif

// ==== Dynamical System ====

extern "C" DLLEXPORT void computef(double time, unsigned int sizeOfX, double* x, double* f,
                                   unsigned int sizeZ, double* z);
extern "C" DLLEXPORT void computef(double time, unsigned int sizeOfX, double* x, double* f,
                                   unsigned int sizeZ, double* z) {
  for (unsigned int i = 0; i < sizeOfX; i++) f[i] = time * x[i];
}

extern "C" DLLEXPORT void computeJacobianfx(double time, unsigned int sizeOfX, double* x,
                                            double* jacob, unsigned int sizeZ, double* z);
extern "C" DLLEXPORT void computeJacobianfx(double time, unsigned int sizeOfX, double* x,
                                            double* jacob, unsigned int sizeZ, double* z) {
  for (unsigned int i = 0; i < sizeOfX * sizeOfX; i++) jacob[i] = i + 1;
}

extern "C" DLLEXPORT void computeM(double time, unsigned int sizeOfX, double* x, double* M,
                                   unsigned int sizeZ, double* z);
extern "C" DLLEXPORT void computeM(double time, unsigned int sizeOfX, double* x, double* M,
                                   unsigned int sizeZ, double* z) {
  for (unsigned int i = 0; i < sizeOfX * sizeOfX; ++i) M[i] = 0;
  M[0] = 1 * time;
  M[4] = 2 * time;
  M[8] = 3 * time;
}

//==================  FirstOrderLinearDS ==================

extern "C" DLLEXPORT void computeb(double time, unsigned int sizeOfB, double* b,
                                   unsigned int sizeOfZ, double* z);
extern "C" DLLEXPORT void computeb(double time, unsigned int sizeOfB, double* b,
                                   unsigned int sizeOfZ, double* z) {
  for (unsigned int i = 0; i < sizeOfB; i++) b[i] = time * i;
}
extern "C" DLLEXPORT void computeA(double time, unsigned int rowA, unsigned int colA,
                                   double* A, unsigned int sizeOfZ, double* z);
extern "C" DLLEXPORT void computeA(double time, unsigned int rowA, unsigned int colA,
                                   double* A, unsigned int sizeOfZ, double* z) {
  for (unsigned int j = 0; j < rowA; j++) {
    for (unsigned int i = 0; i < rowA; i++) A[i + j * rowA] = 4 * (i + 1);
  }
}

//==================  LagrangianScleronomousR ==================
extern "C" DLLEXPORT void hSclero(unsigned int sizeDS, double* q, unsigned int sizeY,
                                  double* y, unsigned int sizeZ, double* z);
extern "C" DLLEXPORT void hSclero(unsigned int sizeDS, double* q, unsigned int sizeY,
                                  double* y, unsigned int sizeZ, double* z) {
  printf("Call of the function 'hSclero' of the test plugin.\n");
}

extern "C" DLLEXPORT void G0Sclero(unsigned int sizeDS, double* q, unsigned int sizeY,
                                   double* G0, unsigned int sizeZ, double* z);
extern "C" DLLEXPORT void G0Sclero(unsigned int sizeDS, double* q, unsigned int sizeY,
                                   double* G0, unsigned int sizeZ, double* z) {
  printf("Call of the function 'G0' of the test plugin.\n");
}

//==================  LagrangianRheonomousR ==================

extern "C" DLLEXPORT void hRheo(unsigned int, double*, double, unsigned int, double*,
                                unsigned int, double*);
extern "C" DLLEXPORT void hRheo(unsigned int, double*, double, unsigned int, double*,
                                unsigned int, double*) {
  printf("Call of the function 'hRheo' of the test plugin.\n");
}

extern "C" DLLEXPORT void G0Rheo(unsigned int, double*, double, unsigned int, double*,
                                 unsigned int, double*);
extern "C" DLLEXPORT void G0Rheo(unsigned int, double*, double, unsigned int, double*,
                                 unsigned int, double*) {
  printf("Call of the function 'G0Rheo' of the test plugin.\n");
}

extern "C" DLLEXPORT void hDot(unsigned int, double*, double, unsigned int, double*,
                               unsigned int, double*);
extern "C" DLLEXPORT void hDot(unsigned int, double*, double, unsigned int, double*,
                               unsigned int, double*) {
  printf("Call of the function 'hDot' of the test plugin.\n");
}

//==================  LagrangianCompliantR ==================

extern "C" DLLEXPORT void hCompl(unsigned int, double*, unsigned int, double*, double*,
                                 unsigned int, double*);
extern "C" DLLEXPORT void hCompl(unsigned int, double*, unsigned int, double*, double*,
                                 unsigned int, double*) {
  printf("Call of the function 'hCompl' of the test plugin.\n");
}


//==================  FirstOrderLinearR ==================

extern "C" DLLEXPORT void C(double, unsigned int, unsigned int, double*, unsigned int,
                            double*);
extern "C" DLLEXPORT void C(double, unsigned int, unsigned int, double*, unsigned int,
                            double*) {
  printf("Warning: call of the function 'C' of the test plugin.\n");
}

extern "C" DLLEXPORT void D(double, unsigned int, unsigned int, double*, unsigned int,
                            double*);
extern "C" DLLEXPORT void D(double, unsigned int, unsigned int, double*, unsigned int,
                            double*) {
  printf("Warning: call of the function 'D' of the test plugin.\n");
}

extern "C" DLLEXPORT void F(double, unsigned int, unsigned int, double*, unsigned int,
                            double*);
extern "C" DLLEXPORT void F(double, unsigned int, unsigned int, double*, unsigned int,
                            double*) {
  printf("Warning: call of the function 'F' of the test plugin.\n");
}

extern "C" DLLEXPORT void e(double, unsigned int, double*, unsigned int, double*);
extern "C" DLLEXPORT void e(double, unsigned int, double*, unsigned int, double*) {
  printf("Warning: call of the function 'e' of the test plugin.\n");
}

extern "C" DLLEXPORT void B(double, unsigned int, unsigned int, double*, unsigned int,
                            double*);
extern "C" DLLEXPORT void B(double, unsigned int, unsigned int, double*, unsigned int,
                            double*) {
  printf("Warning: call of the function 'B' of the test plugin.\n");
}
