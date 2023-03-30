/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2023 INRIA.
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

/*! \file
  C++ API (iso_c_binding) to fortran solvers (netlib/odepack, hairer/hem5)

  All fortran functions are supposed to be defined in the module siconos_fortran2c
  (siconos_fortran2c.f90)

*/

#ifndef SICONOSFORTRAN_H
#define SICONOSFORTRAN_H

namespace siconos::netlib {

typedef void (*fpointer)(int*, double*, double*, double*);
typedef void (*gpointer)(int*, double*, double*, int*, double*);
typedef void (*jacopointer)(int*, double*, double*, int*, int*, double*, int*);

#if defined(HAS_FORTRAN)
extern "C" void lsodar(fpointer, int* NEQ, double* Y, double* T, double* TOUT, int* ITOL,
                       double* RTOL, double* ATOL, int* ITASK, int* ISTATE, int* IOPT,
                       double* RWORK, int* LRW, int* IWORK, int* LIW, jacopointer C_JAC,
                       int* JT, gpointer C_G, int* NG, int* JROOT);

#else
extern "C" inline void lsodar(fpointer, int* NEQ, double* Y, double* T, double* TOUT,
                              int* ITOL, double* RTOL, double* ATOL, int* ITASK, int* ISTATE,
                              int* IOPT, double* RWORK, int* LRW, int* IWORK, int* LIW,
                              jacopointer C_JAC, int* JT, gpointer C_G, int* NG, int* JROOT) {
  printf("Siconos Fortran API is off. This function (lsodar) has no effects.\n");
}
#endif

}  // namespace siconos::netlib

namespace siconos::hairer {

typedef void (*fprobpointer)(int* IFCN, int* NQ, int* NV, int* NU, int* NL, int* LDG, int* LDF,
                             int* LDA, int* NBLK, int* NMRC, int* NPGP, int* NPFL, int* INDGR,
                             int* INDGC, int* INDFLR, int* INDFLC, double* time, double* q,
                             double* v, double* u, double* xl, double* G, double* GQ,
                             double* F, double* GQQ, double* GT, double* FL, double* QDOT,
                             double* UDOT, double* AM);

typedef void (*soloutpointer)(int* MODE, int* NSTEP, int* NQ, int* NV, int* NU, int* NL,
                              int* LDG, int* LDF, int* LDA, int* LRDO, int* LIDO,
                              fprobpointer FPROB, double* q, double* v, double* u,
                              double* DOWK, int* IDOWK);

#if defined(HAS_FORTRAN)
extern "C" void hem5(int* NQ, int* NV, int* NU, int* NL, fprobpointer FPROB, double* T,
                     double* Q, double* V, double* U, double* A, double* RLAM, double* TEND,
                     double* H, double* RTOL, double* ATOL, int* ITOL, soloutpointer SOLOUT,
                     int* IOUT, double* WK, int* LWK, int* IWK, int* LIWK, int* IDID);
#else
// #if !defined(HAS_FORTRAN)
extern "C" inline void hem5(int* NQ, int* NV, int* NU, int* NL, fprobpointer FPROB, double* T,
                            double* Q, double* V, double* U, double* A, double* RLAM,
                            double* TEND, double* H, double* RTOL, double* ATOL, int* ITOL,
                            soloutpointer SOLOUT, int* IOUT, double* WK, int* LWK, int* IWK,
                            int* LIWK, int* IDID) {
  printf("Siconos Fortran API is off. This function (hem5) has no effects.\n");
}

#endif
}  // namespace siconos::hairer

#endif
