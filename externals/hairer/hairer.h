// Siconos-Numerics, Copyright INRIA 2005-2011.
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.
// Siconos is a free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// Siconos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Siconos; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//
// Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
//

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
