/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <assert.h>

#include "gfc3d_Solvers.h"
#include "NonSmoothDrivers.h"
#include "fc3d_Solvers.h"
#include "fc3d_solvers_wr.h"
#include <string.h>

#include "NumericsSparseMatrix.h"
#include "NumericsVector.h"

#include "sanitizer.h"
#include "numerics_verbose.h"
//#define TEST_COND
/* #define OUTPUT_DEBUG */


/* #define DEBUG_NOCOLOR */
/* #define DEBUG_MESSAGES */
/* #define DEBUG_STDOUT */
#include "debug.h"

//#define USE_LAPACK_DGETRS

#pragma GCC diagnostic ignored "-Wmissing-prototypes"

/* Global Variable for the reformulation of the problem */

GlobalFrictionContactProblem*  fc3d_reformulation_global_problem(FrictionContactProblem* problem)
{

  int dimension  =  problem->dimension;
  int nc  = problem->numberOfContacts;
  int m = dimension * nc;



  GlobalFrictionContactProblem*  globalproblem = (GlobalFrictionContactProblem*) malloc(sizeof(GlobalFrictionContactProblem));
  globalFrictionContact_null(globalproblem);


  globalproblem->numberOfContacts = problem->numberOfContacts;
  globalproblem->dimension =  problem->dimension;

  globalproblem->mu = (double *)malloc(m*sizeof(double));
  cblas_dcopy(nc,problem->mu,1,globalproblem->mu,1);
  globalproblem->q = (double *)malloc(m*sizeof(double));
  cblas_dcopy(m,problem->q,1,globalproblem->q,1);
  globalproblem->b = (double *)calloc(m,sizeof(double));

  NumericsMatrix *M = problem->M;


  NumericsMatrix *Mglobal = NM_create(NM_SPARSE,m,m);
  CSparseMatrix * M_triplet = NM_triplet(M);
  Mglobal->matrix2->triplet = M_triplet;
  Mglobal->matrix2->origin = NSM_TRIPLET;
  globalproblem->M=Mglobal;

  NumericsMatrix *Hglobal = NM_eye(m);
  NM_display(Hglobal);
  globalproblem->H=M=Hglobal;



   return globalproblem;
}
