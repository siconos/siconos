/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include "SiconosConfig.h"
#ifdef HAS_FORTRAN
#include "SiconosFortran.h"

typedef void fprobfunction (integer* IFCN,
                             integer* NQ,
                             integer* NV,
                             integer* NU,
                             integer* NL,
                             integer* LDG, integer* LDF, integer* LDA,
                             integer* NBLK, integer* NMRC,
                             integer* NPGP, integer* NPFL,
                             integer* INDGR, integer* INDGC, integer * INDFLR, integer * INDFLC,
                             doublereal* time,
                             doublereal* q, doublereal* v, doublereal* u,  doublereal* xl,
                             doublereal* G, doublereal* GQ, doublereal * F,
                             doublereal* GQQ, doublereal* GT, doublereal * FL,
                             doublereal* QDOT, doublereal* UDOT, doublereal * AM);

typedef fprobfunction* fprobpointer;

typedef void soloutfunction(integer* MODE,
                              integer* NSTEP,
                              integer* NQ,
                              integer* NV,
                              integer* NU,
                              integer* NL,
                              integer* LDG, integer* LDF, integer* LDA,
                              integer* LRDO, integer* LIDO,
                              fprobpointer FPROB,
                              doublereal* q, doublereal* v, doublereal* u,
                              doublereal *DOWK, integer* IDOWK);

typedef soloutfunction* soloutpointer;

#ifdef __cplusplus
extern "C" {
#endif
  void CNAME(hem5)(integer* NQ, integer* NV, integer* NU, integer* NL,
                     fprobpointer FPROB,
                     doublereal* T,
                     doublereal* Q, doublereal* V, doublereal* U,  doublereal* A,  doublereal* RLAM,
                     doublereal* TEND, doublereal * H,
                     doublereal* RTOL, doublereal* ATOL, integer* ITOL,
                     soloutpointer SOLOUT, integer * IOUT,
                     doublereal* WK, integer * LWK, integer* IWK, integer* LIWK,
                     integer* IDID);


#ifdef __cplusplus
}
#endif

#endif
