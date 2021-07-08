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
#include <assert.h>                  // for assert
#include <stdio.h>                   // for printf, NULL, fprintf, stderr
#include <stdlib.h>                  // for exit, malloc, EXIT_FAILURE
#include "FrictionContactProblem.h"  // for FrictionContactProblem
#include "Friction_cst.h"            // for SICONOS_FRICTION_2D_NSGS, SICONO...
#include "NonSmoothDrivers.h"        // for fc2d_driver
#include "NumericsFwd.h"             // for FrictionContactProblem, SolverOp...
#include "NumericsMatrix.h"          // for NumericsMatrix, RawNumericsMatrix
#include "SolverOptions.h"           // for SolverOptions, solver_options_id...
#include "SparseBlockMatrix.h"       // for SparseBlockStructuredMatrix, SBM...
#include "fc2d_Solvers.h"            // for fc2d_cpg, fc2d_enum
#include "numerics_verbose.h"        // for numerics_error, verbose, numeric...

const char* const   SICONOS_FRICTION_2D_NSGS_STR  = "FC2D_NSGS";
const char* const   SICONOS_FRICTION_2D_CPG_STR  = "FC2D_CPG";
const char* const   SICONOS_FRICTION_2D_LEMKE_STR  = "FC2D_LEMKE";
const char* const   SICONOS_FRICTION_2D_ENUM_STR  = "FC2D_ENUM";
//#define DUMP_PROBLEM
#ifdef DUMP_PROBLEM
static int fccounter = 0;
#endif
//#define DUMP_PROBLEM_IF_INFO
#ifdef DUMP_PROBLEM_IF_INFO
static int fccounter = 0;
#endif


int fc2d_driver(FrictionContactProblem* problem, double *reaction, double *velocity, SolverOptions* options)
{

#ifdef DUMP_PROBLEM
  char fname[256];
  int ncc = problem->numberOfContacts;
  sprintf(fname, "fc2d_granularflowonwall_%.5d_%.5d.dat", ncc, fccounter++);
  printf("Dump %s file\n", fname);

  FILE * foutput  =  fopen(fname, "w");
  frictionContact_printInFile(problem, foutput);
  fclose(foutput);
#endif
  
  if(options == NULL)
    numerics_error("fc2d_driver", "null input for solver options");

  /* Checks inputs */
  if(problem == NULL || reaction == NULL || velocity == NULL)
    numerics_error("fc2d_driver", "null input for FrictionContactProblem and/or unknowns (reaction,velocity)");

  assert(options->isSet);

  if(verbose > 0)
    solver_options_print(options);


  /* Solver name */
  /*const char* const  name = options->solverName;*/


  int info = -1 ;

  if(problem->dimension != 2)
    numerics_error("fc2d_driver", "Dimension of the problem : problem-> dimension is not compatible or is not set");


  /* Non Smooth Gauss Seidel (NSGS) */

  if(problem->M->storageType == NM_SPARSE_BLOCK || problem->M->storageType == NM_SPARSE )
  {

    if(options->solverId == SICONOS_FRICTION_2D_NSGS)
    {
      numerics_printf(" ======================= Call Sparse NSGS solver for Friction-Contact 2D problem ======================");
      fc2d_nsgs(problem, reaction, velocity, &info, options);
    }
    else
    {

      NumericsMatrix* M_dense = NM_create(NM_DENSE, problem->M->size0, problem->M->size1);
      NM_to_dense(problem->M, M_dense);
      NumericsMatrix* M_original = problem->M;
      problem->M = M_dense;

      info = fc2d_driver(problem, reaction, velocity, options);

      NM_free(M_dense);
      problem->M = M_original;
    }

  }
  else if(problem->M->storageType == NM_DENSE)
  {

    switch(options->solverId)
    {
    /****** NLGS algorithm ******/
    case SICONOS_FRICTION_2D_NSGS:
    {
      if(verbose)
        printf(" ========================== Call NLGS solver for Friction-Contact 2D problem ==========================\n");
      fc2d_nsgs_dense(problem, reaction, velocity, &info, options);
      break;
    }
    /****** CPG algorithm ******/
    case SICONOS_FRICTION_2D_CPG:
    {
      if(verbose)
        printf(" ========================== Call CPG solver for Friction-Contact 2D problem ==========================\n");
      fc2d_cpg(problem, reaction, velocity, &info, options);
      break;
    }
    /****** Lexicolemke algorithm ******/
    case SICONOS_FRICTION_2D_LEMKE:
    {
      if(verbose)
        printf(" ========================== Call Lemke solver for Friction-Contact 2D problem ==========================\n");
      fc2d_lexicolemke(problem, reaction, velocity, &info, options);
      break;
    }
    /****** Enum algorithm ******/
    case SICONOS_FRICTION_2D_ENUM:
    {
      if(verbose)
        printf(" ========================== Call Enumerative solver for Friction-Contact 2D problem ==========================\n");
      fc2d_enum(problem, reaction, velocity, &info, options);
      break;
    }
    /*error */
    default:
    {
      fprintf(stderr, "fc2d_driver error: unknown solver named: %s\n", solver_options_id_to_name(options->solverId));
      exit(EXIT_FAILURE);
    }
    }
#ifdef DUMP_PROBLEM_IF_INFO
    if(info)
    {
      char fname[256];
      sprintf(fname, "FrictionContactProblem%.5d.dat", fccounter++);
      printf("Dump of FrictionContactProblem%.5d.dat\n", fccounter);

      FILE * foutput  =  fopen(fname, "w");
      frictionContact_printInFile(problem, foutput);
      fclose(foutput);
    }
#endif
  }
  else
  {
    numerics_error("fc2d_driver",
                   " error: unknown storagetype named");
    exit(EXIT_FAILURE);
  }

  return info;

}
