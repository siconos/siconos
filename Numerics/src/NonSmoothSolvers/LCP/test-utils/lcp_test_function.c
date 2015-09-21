/* Siconos-Numerics, Copyright INRIA 2005-2012.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include <stdio.h>
#include <stdlib.h>
#include "NonSmoothDrivers.h"
#include "lcp_test_function.h"
#include "GAMSlink.h"

void fillParamWithRespectToSolver_SBM(SolverOptions *options, int solverId, LinearComplementarityProblem* problem)
{
  int maxIter = 1001;
  double tolerance = 1e-8;
  double lighttolerance = 1e-8;

  switch (solverId)
  {
  case SICONOS_LCP_PGS:
  case SICONOS_LCP_CPG:
  case SICONOS_LCP_LEMKE:
  case SICONOS_LCP_PIVOT:
  case SICONOS_LCP_BARD:
  case SICONOS_LCP_MURTY:
  case SICONOS_LCP_AVI_CAOFERRIS:
  case SICONOS_LCP_NEWTONMIN:
  case SICONOS_LCP_NEWTON_FBLSA:
  case SICONOS_LCP_NEWTON_MINFBLSA:
  {
    options->iparam[0] = maxIter;
    options->dparam[0] = tolerance;
    break;
  }
  case SICONOS_LCP_RPGS:
  {
    options->iparam[0] = maxIter;
    options->dparam[0] = tolerance;
    options->dparam[2] = 1.0;
    break;
  }
  case SICONOS_LCP_LATIN :
  {
    options->iparam[0] = maxIter;
    options->dparam[0] = lighttolerance;
    options->dparam[2] = 1.0;
    break;
  }
  case SICONOS_LCP_LATIN_W:
  {
    options->iparam[0] = maxIter;
    options->dparam[0] = lighttolerance;
    options->dparam[2] = 0.3;
    options->dparam[3] = 1.0;
    break;
  }
  case SICONOS_LCP_PATH:
  case SICONOS_LCP_QP:
  case SICONOS_LCP_NSQP:
  {
    options->dparam[0] = tolerance;
    break;
  }
  case SICONOS_LCP_ENUM:
  {
    /*       options->dparam[0]=tolerance; */
    /*       options->dWork=(double*) malloc((3*problem->size +problem->size*problem->size)*sizeof(double)); */
    /*       options->iWork=(int*) malloc(2*problem->size*sizeof(int)); */

    if (options->iparam != NULL)
      free(options->iparam);
    if (options->dparam != NULL)
      free(options->dparam);

    linearComplementarity_enum_setDefaultSolverOptions(problem,  options);

    break;
  }
  default:
    ;
  }


}



int lcp_test_function(FILE * f, int solverId, char* filename)
{

  int i, info = 0 ;
  LinearComplementarityProblem* problem = (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));

  info = linearComplementarity_newFromFile(problem, f);

  FILE * foutput  =  fopen("./lcp_mmc.verif", "w");
  info = linearComplementarity_printInFile(problem, foutput);
  fclose(foutput);

  NumericsOptions global_options;
  setDefaultNumericsOptions(&global_options);
  global_options.verboseMode = 1;
  SolverOptions options;
  set_SolverOptions(&options, solverId);

#ifdef HAVE_GAMS_C_API
  if (solverId == SICONOS_LCP_GAMS)
  {
    SN_GAMSparams* GP = (SN_GAMSparams*)options.solverParameters;
    assert(GP);
    GP->model_dir = GAMS_MODELS_SOURCE_DIR;
    assert(filename);
    GP->filename = filename;
  }
#endif

  double * z = (double *)calloc(problem->size, sizeof(double));
  double * w = (double *)calloc(problem->size, sizeof(double));

  info = linearComplementarity_driver(problem, z , w, &options, &global_options);

  for (i = 0 ; i < problem->size ; i++)
  {
    printf("z[%i] = %12.8e\t,w[%i] = %12.8e\n", i, z[i], i, w[i]);
  }

  if (!info)
  {
    printf("test succeeded err = %e \n", options.dparam[1]);
  }
  else
  {
    printf("test unsuccessful err =%e  \n", options.dparam[1]);
  }
  free(z);
  free(w);

  deleteSolverOptions(&options);

  if (solverId == SICONOS_LCP_GAMS)
  {
    free(options.solverParameters);
    options.solverParameters = NULL;
  }

  freeLinearComplementarityProblem(problem);
  printf("End of test.\n");


  return info;
}

int lcp_test_function_SBM(FILE * f, int solverId)
{

  int i, info = 0 ;
  LinearComplementarityProblem* problem = (LinearComplementarityProblem *)malloc(sizeof(LinearComplementarityProblem));

  info = linearComplementarity_newFromFile(problem, f);

  FILE * foutput  =  fopen("./lcp_mmc.verif", "w");
  info = linearComplementarity_printInFile(problem, foutput);


  NumericsOptions global_options;
  setDefaultNumericsOptions(&global_options);
  global_options.verboseMode = 1;

  SolverOptions * options = (SolverOptions *)malloc(sizeof(SolverOptions));



  info = linearComplementarity_setDefaultSolverOptions(problem, options, SICONOS_LCP_NSGS_SBM);

  options->internalSolvers->solverId = solverId;

  fillParamWithRespectToSolver_SBM(options->internalSolvers, solverId, problem);

#ifdef HAVE_GAMS_C_API
  if (solverId == SICONOS_LCP_GAMS)
  {
    // no testing for now
    deleteSolverOptions(options);
    free(options);
    freeLinearComplementarityProblem(problem);
    fclose(foutput);
    return 0;
/*    SN_GAMSparams* GP = (SN_GAMSparams*)options->internalSolvers->solverParameters;
    assert(GP);
    GP->model_dir = GAMS_MODELS_SOURCE_DIR;*/
  }
#endif



  double * z = (double *)calloc(problem->size, sizeof(double));
  double * w = (double *)calloc(problem->size, sizeof(double));

  info = linearComplementarity_driver(problem, z , w, options, &global_options);

  for (i = 0 ; i < problem->size ; i++)
  {
    printf("z[%i] = %12.8e\t,w[%i] = %12.8e\n", i, z[i], i, w[i]);
  }

  if (!info)
  {
    printf("test succeeded err=%e \n", options->dparam[1]);
  }
  else
  {
    printf("test unsuccessful err =%e \n", options->dparam[1]);
  }
  free(z);
  free(w);
  // info = linearComplementarity_deleteDefaultSolverOptions(&options,solvername);

  deleteSolverOptions(options);
  free(options);

  freeLinearComplementarityProblem(problem);
  fclose(foutput);

  return info;


}


