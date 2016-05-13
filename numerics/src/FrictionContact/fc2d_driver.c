/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#include "NumericsOptions.h"
#include "fc2d_Solvers.h"
#include "NonSmoothDrivers.h"
char *  SICONOS_FRICTION_2D_NSGS_STR  = "F2D_NSGS";
char *  SICONOS_FRICTION_2D_PGS_STR  = "F2D_PGS";
char *  SICONOS_FRICTION_2D_CPG_STR  = "F2D_CPG";
char *  SICONOS_FRICTION_2D_LATIN_STR  = "F2D_LATIN";
char *  SICONOS_FRICTION_2D_LEMKE_STR  = "F2D_LEMKE";
char *  SICONOS_FRICTION_2D_ENUM_STR  = "F2D_ENUM";
//#define DUMP_PROBLEM
#ifdef DUMP_PROBLEM
static int fccounter = 0;
#endif
//#define DUMP_PROBLEM_IF_INFO
#ifdef DUMP_PROBLEM_IF_INFO
static int fccounter = 0;
#endif


int fc2d_driver(FrictionContactProblem* problem, double *reaction , double *velocity, SolverOptions* options, NumericsOptions* global_options)
{

#ifdef DUMP_PROBLEM
  char fname[256];
  sprintf(fname, "FrictionContactProblem%.5d.dat", fccounter++);
  printf("Dump of FrictionContactProblem%.5d.dat", fccounter);
  
  FILE * foutput  =  fopen(fname, "w");
  frictionContact_printInFile(problem, foutput);
  fclose(foutput);
#endif

  if (options == NULL || global_options == NULL)
    numericsError("fc2d_driver", "null input for solver and/or global options");

  /* Set global options */
  setNumericsOptions(global_options);

  /* Checks inputs */
  if (problem == NULL || reaction == NULL || velocity == NULL)
    numericsError("fc2d_driver", "null input for FrictionContactProblem and/or unknowns (reaction,velocity)");

  /* If the options for solver have not been set, read default values in .opt file */
  int NoDefaultOptions = options->isSet; /* true(1) if the SolverOptions structure has been filled in else false(0) */

  if (!NoDefaultOptions)
    readSolverOptions(3, options);

  if (verbose > 0)
    printSolverOptions(options);





  /* Solver name */
  /*char * name = options->solverName;*/


  int info = -1 ;

  if (problem->dimension != 2)
    numericsError("fc2d_driver", "Dimension of the problem : problem-> dimension is not compatible or is not set");


  /* Non Smooth Gauss Seidel (NSGS) */

  if (problem->M->storageType == 1)
  {

    if (options->solverId == SICONOS_FRICTION_2D_NSGS)
    {
      if (verbose)
        printf(" ======================= Call Sparse NSGS solver for Friction-Contact 2D problem ======================\n");
      fc2d_sparse_nsgs(problem, reaction , velocity , &info , options);
    }
    else
    {
      fprintf(stderr, "fc2d_driver error: unknown solver named: %s\n", idToName(options->solverId));
      exit(EXIT_FAILURE);
    }
  }
  else if (problem->M->storageType == 0)
  {

    switch (options->solverId)
    {
      /****** NLGS algorithm ******/
    case SICONOS_FRICTION_2D_PGS:
    case SICONOS_FRICTION_2D_NSGS:
    {
      if (verbose)
        printf(" ========================== Call NLGS solver for Friction-Contact 2D problem ==========================\n");
      fc2d_nsgs(problem, reaction, velocity, &info, options);
      break;
    }
    /****** CPG algorithm ******/
    case SICONOS_FRICTION_2D_CPG:
    {
      if (verbose)
        printf(" ========================== Call CPG solver for Friction-Contact 2D problem ==========================\n");
      fc2d_cpg(problem, reaction, velocity, &info, options);
      break;
    }
    /****** Latin algorithm ******/
    case SICONOS_FRICTION_2D_LATIN:
    {
      if (verbose)
        printf(" ========================== Call Latin solver for Friction-Contact 2D problem ==========================\n");
      fc2d_latin(problem, reaction, velocity, &info, options);
      break;
    }
    /****** Lexicolemke algorithm ******/
    case SICONOS_FRICTION_2D_LEMKE:
    {
      if (verbose)
        printf(" ========================== Call Lemke solver for Friction-Contact 2D problem ==========================\n");
      fc2d_lexicolemke(problem, reaction, velocity, &info, options, global_options);
      break;
    }
    /****** Enum algorithm ******/
    case SICONOS_FRICTION_2D_ENUM:
    {
      if (verbose)
        printf(" ========================== Call Enumerative solver for Friction-Contact 2D problem ==========================\n");
      fc2d_enum(problem, reaction, velocity, &info, options, global_options);
      break;
    }
    /*error */
    default:
    {
      fprintf(stderr, "fc2d_driver error: unknown solver named: %s\n", idToName(options->solverId));
      exit(EXIT_FAILURE);
    }
    }
#ifdef DUMP_PROBLEM_IF_INFO
    if (info)
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
    numericsError("fc2d_driver", 
                  " error: unknown storagetype named");
    exit(EXIT_FAILURE);
  }

  return info;

}
