/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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
#ifndef ODEPACK_H
#define ODEPACK_H


#include "SiconosConfig.h"
#ifdef HAS_FORTRAN
#include "SiconosFortran.h"

typedef void (*fpointer)(integer *, doublereal *, doublereal *, doublereal *);
typedef void (*gpointer)(integer *, doublereal *, doublereal*, integer *, doublereal *);
typedef void (*jacopointer)(integer *, doublereal *, doublereal *, integer* , integer *,  doublereal *, integer *);

#ifdef __cplusplus
extern "C" {
#endif

  void CNAME(dlsode)(fpointer, integer * neq, doublereal * y, doublereal *t, doublereal *tout, integer * itol, doublereal * rtol, doublereal *atol, integer * itask, integer *istate, integer * iopt, doublereal * rwork, integer * lrw, integer * iwork, integer * liw, jacopointer, integer * mf);

  void CNAME(dlsodar)(fpointer, integer * neq, doublereal * y, doublereal *t, doublereal *tout, integer * itol, doublereal * rtol, doublereal *atol, integer * itask, integer *istate, integer * iopt, doublereal * rwork, integer * lrw, integer * iwork, integer * liw, jacopointer, integer * jt, gpointer, integer* ng, integer * jroot);

#ifdef __cplusplus
}
#endif

#endif

#endif
