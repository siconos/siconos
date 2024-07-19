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
                       double* RTOL, double* ATOL, int* ISTATE, double* RWORK, int* LRW,
                       int* IWORK, int* LIW, jacopointer C_JAC, int* JT, gpointer C_G, int* NG,
                       int* JROOT);

#else
extern "C" inline void lsodar(fpointer, int* NEQ, double* Y, double* T, double* TOUT,
                              int* ITOL, double* RTOL, double* ATOL, int* ISTATE,
                              double* RWORK, long* LRW, int* IWORK, int* LIW,
                              jacopointer C_JAC, int* JT, gpointer C_G, int* NG, int* JROOT) {
  printf("Siconos Fortran API is off. This function (lsodar) has no effects.\n");
}
#endif

}  // namespace siconos::netlib

namespace siconos::fortran {
namespace hairer {

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
}  // namespace hairer

namespace optim {
#if defined(HAS_FORTRAN)

/** Computes the minimum of a constrained function

    see
   https://who.rocq.inria.fr/Jean-Charles.Gilbert/modulopt/optimization-routines/m2qn1/m2qn1.pdf

    \param[in] n number of variables on which f depends
    \param[in, out] x starting/final point, size n
    \param[in, out] function value at x
    \param[in, out] g gradient of, evaluated at x
    \param[in], dxmin, vector of size n, used to control precision among other things
    \param[in,out] df1 see paper
    \param[in, out] epsabs convergence criteria
    \param[in, out] mode control the way the algo is initialized (in) and give details
     about how it stops (out)
    \param[in] binf constraint on x (lower bound)
    \param[in] bsup constraint on x (upper bound)
    \param[in, out] iz work vector (size = 2n+1)
    \param[in, out] rz work vector (size = 1/2*n*(n+9))
*/
extern "C" void n2qn1(const int* n, double* x, double* f, double* g, double* dxmin,
                      double* df1, double* epsabs, int* mode, const double* binf,
                      const double* bsup, int* iz, double* rz);

#else
extern "C" inline void n2qn1(int* n, double* x, double* f, double* g, double* dxmin,
                             double* df1, double* epsabs, int* mode, const double* binf,
                             const double* bsup, int* iz, double* rz) {
  printf("Siconos Fortran API is off. This function (n2qn1) has no effects.\n");
}
#endif

}  // namespace optim
}  // namespace siconos::fortran
#endif
