/* Siconos-Numerics version 2.1.1, Copyright INRIA 2005-2007.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

#include "FrictionContact2D_Solvers.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef MEXFLAG
#include "NonSmoothDrivers.h"
#endif

int dfc_2D_driver(FrictionContact_Problem* problem, double *reaction , double *velocity, Solver_Options* options, Numerics_Options* global_options, int* iparamDFC, double* J1)
{
  if (options == NULL || global_options == NULL)
    numericsError("dfc_2D_driver", "null input for solver and/or global options");

  /* Set global options */
  setNumericsOptions(global_options);

  /* Checks inputs */
  if (problem == NULL || reaction == NULL || velocity == NULL)
    numericsError("dfc_2D_driver", "null input for FrictionContact_Problem and/or unknowns (reaction,velocity)");

  /* Output info. : 0: ok -  >0: problem (depends on solver) */
  int info = -1;

  int storageType = problem->M->storageType;
  if (storageType == 1)
    numericsError("dfc_2D_driver", "not yet implemented for Sparse Block Storage");

  /* If the options for solver have not been set, read default values in .opt file */
  int NoDefaultOptions = options->isSet; /* true(1) if the Solver_Options structure has been filled in else false(0) */

  if (NoDefaultOptions == 0)
  {
    readSolverOptions(1, options);
    options->filterOn = 1;
  }

  if (verbose > 0)
    printSolverOptions(options);

  /*************************************************
   *  1 - Call specific solver
   *************************************************/

  /* Solver name */
  char * name = options->solverName;

  if (verbose == 1)
    printf(" ========================== Call %s solver for dual Friction Contact 2D problem ==========================\n", name);

  int i;
  int n = problem->numberOfContacts * 2;
  /****** CFD Latin algorithm ******/
  if (strcmp(name , "Cfd_latin") == 0)
  {
    FrictionContact_Problem newProblem;
    newProblem.numberOfContacts = iparamDFC[3];
    int dim_q  = 2 * iparamDFC[3];
    int dim_MM = dim_q * dim_q;
    newProblem.M = malloc(sizeof(*newProblem.M));
    newProblem.M->matrix0 = malloc(dim_MM * sizeof(double));
    newProblem.q = malloc(dim_q * sizeof(double));
    newProblem.mu = problem->mu;
    double * z = malloc(dim_q * sizeof(double));
    double * w = malloc(dim_q * sizeof(double));
    for (i = 0; i < dim_q; i++)
    {
      problem->q[i] = 0.0;
      z[i] = 0.0;
      w[i] = 0.0;
    }

    dfc_2D2cond_2D(&n , problem->mu , problem->M->matrix0 , problem->q, &iparamDFC[0] , &iparamDFC[1] , &iparamDFC[3] ,
                   &iparamDFC[2] , &iparamDFC[4] ,  J1 , newProblem.M->matrix0 , newProblem.q);

    dfc_2D_latin(&newProblem, z, w, &info, options);

    cond_2D2dfc_2D(&n , z , w ,  problem->M->matrix0 , problem->q , J1 ,  &iparamDFC[0] , &iparamDFC[1] , &iparamDFC[3] ,
                   &iparamDFC[2] , &iparamDFC[4], reaction , velocity);

    free(z);
    free(w);
    newProblem.mu = NULL;
    free(newProblem.q);
    free(newProblem.M->matrix0);
    free(newProblem.M);
  }


  /****** Lemke algorithm ******/
  else if (strcmp(name , "Lemke")  == 0 || (strcmp(name , "PGS")) == 0 || (strcmp(name , "CPG")) == 0)
  {
    /* Convert fc2D to LCP */
    LinearComplementarity_Problem LCP;
    LCP.size = 3 * iparamDFC[3];
    int dim_q  =  LCP.size ;
    int dim_MM = dim_q * dim_q;
    LCP.M = malloc(sizeof(*LCP.M));
    LCP.M->matrix0 = malloc(dim_MM * sizeof(double));
    LCP.q = malloc(dim_q * sizeof(double));

    double *z      = (double *)malloc(dim_q * sizeof(double));
    double* w      = (double *)malloc(dim_q * sizeof(double));

    for (int i = 0; i < dim_q; i++)
    {
      LCP.q[i] = 0.0;
      z[i] = 0.0;
      w[i] = 0.0;
    }

    dfc_2D2lcp(&n , problem->mu , problem->M->matrix0 , problem->q, &iparamDFC[0] , &iparamDFC[1] , &iparamDFC[3] ,
               &iparamDFC[2] , &iparamDFC[4] ,  J1 , LCP.M->matrix0 , LCP.q);

    /* Options for the LCP solver */
    Solver_Options optionsLCP;
    strcpy(optionsLCP.solverName, name);
    int iparam[2] = {options->iparam[0], 0};
    double dparam[2] = {options->dparam[0], 1.0};
    optionsLCP.iSize = 2;
    optionsLCP.dSize = 2;
    optionsLCP.iparam = iparam;
    optionsLCP.dparam = dparam ;
    optionsLCP.isSet = 1;
    optionsLCP.filterOn = 1;
    /* Call LCP solver */
    if (strcmp(name , "Lemke") == 0)
    {
      lcp_lexicolemke(&LCP, z, w, &info, &optionsLCP);
    }
    else if (strcmp(name , "PGS") == 0)
    {
      lcp_pgs(&LCP, z, w, &info, &optionsLCP);
      options->dparam[1] = optionsLCP.dparam[1];
      options->iparam[1] = optionsLCP.iparam[1];
    }
    else if (strcmp(name , "CPG") == 0)
    {
      lcp_cpg(&LCP, z, w, &info, &optionsLCP);
      options->dparam[1] = optionsLCP.dparam[1];
      options->iparam[1] = optionsLCP.iparam[1];
    }

    lcp2dfc_2D(&n , z , w ,  problem->M->matrix0 , problem->q , J1 ,  &iparamDFC[0] , &iparamDFC[1] , &iparamDFC[3] ,
               &iparamDFC[2] , &iparamDFC[4],  reaction , velocity);

    free(z);
    free(w);
    free(LCP.q);
    free(LCP.M->matrix0);
    free(LCP.M);
  }

  /*error */
  else
  {
    fprintf(stderr, "dfc_2D_driver error: unknown solver named: %s\n", name);
    exit(EXIT_FAILURE);
  }

  return info;

}

