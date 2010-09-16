/* Siconos-Numerics, Copyright INRIA 2005-2010.
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
#include "SolverOptions.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mlcp_cst.h"
#include "lcp_cst.h"
#include "relay_cst.h"
#include "Friction_cst.h"

char SICONOS_NUMERICS_PROBLEM_LCP_STR[] = "LCP";
char SICONOS_NUMERICS_PROBLEM_MLCP_STR[] = "MLCP";
char SICONOS_NUMERICS_PROBLEM_EQUALITY_STR[] = "EQUALITY";
char SICONOS_NUMERICS_PROBLEM_FC2D_STR[] = "FC2D";
char SICONOS_NUMERICS_PROBLEM_FC3D_STR[] = "FC3D";


char * idProblemToChar(int id)
{
  switch (id)
  {
  case (SICONOS_NUMERICS_PROBLEM_LCP):
  {
    return SICONOS_NUMERICS_PROBLEM_LCP_STR;
    break;
  }
  case (SICONOS_NUMERICS_PROBLEM_MLCP):
  {
    return SICONOS_NUMERICS_PROBLEM_MLCP_STR;
    break;
  }
  case (SICONOS_NUMERICS_PROBLEM_EQUALITY):
  {
    return SICONOS_NUMERICS_PROBLEM_EQUALITY_STR;
    break;
  }
  case (SICONOS_NUMERICS_PROBLEM_FC2D):
  {
    return SICONOS_NUMERICS_PROBLEM_FC2D_STR;
    break;
  }
  case (SICONOS_NUMERICS_PROBLEM_FC3D):
  {
    return SICONOS_NUMERICS_PROBLEM_FC3D_STR;
    break;
  }
  default:
    printf("Numerics:idProblemToChar, id unknown : %d \n", id);
    return NULL;
  }

}

void readSolverOptions(int driverType, SolverOptions* options)
{
  /* To each problem, corresponds a XXX_parameters.opt file where default parameters can be read, XXX being the problem name (LCP, FrictionContact3D ...) */

  if (verbose > 0)
    printf("\n ========== Numerics Non Smooth Solver - Read default parameters for the solver.\n ==========");

  // Checks if NUMERICSSPATH is set.
  if (getenv("SICONOSPATH") == NULL)
  {
    fprintf(stderr, "Numerics, readSolverOptions error, SICONOSPATH environment variable not set. Can not find default solver options file.\n");
    exit(EXIT_FAILURE);
  }

  FILE * ficin;
  /* Name of the default parameters file */
  char name[64];

  strcpy(name, getenv("SICONOSPATH"));
  strcat(name, "/include/Siconos/Numerics/");

  char buffer[64];
  char bufferName[64];
  /* Return value for reading */
  int nval;

  // set default size to 4 ...
  if (options->iparam == NULL)
    options->iparam = (int*)malloc(4 * sizeof(*options->iparam));
  if (options->dparam == NULL)
    options->dparam = (double*)malloc(4 * sizeof(*options->dparam));

  switch (driverType)
  {

  case 0:
    strcat(name, "LCP_parameters.opt");
  case 1:
    strcat(name, "dfc2D_parameters.opt");
  case 2:
    strcat(name, "FrictionContact2D_parameters.opt");
  case 3:
    strcat(name, "FrictionContact3D_parameters.opt");
    ficin = fopen(name, "r");
    if (verbose > 0)
      printf("The default-parameters file is: %s\n", name);
    if (!ficin)
    {
      printf("Numerics, readSolverOptions error. Can not open file %60s", name);
      exit(-1);
    }
    //nval = fscanf(ficin, "%c", &(options->solverName));
    fgets(buffer, 64, ficin);
    fgets(buffer, 64, ficin);
    fgets(buffer, 64, ficin);
    /* Solver name */
    fgets(bufferName , 64, ficin);
    options->solverId = nameToId(bufferName);
    fgets(buffer, 64, ficin);
    /* iparam */
    nval = fscanf(ficin, "%d%d", &(options->iparam[0]), &(options->iparam[1]));
    if (nval != 4)
    {
      fprintf(stderr, "Numerics, readSolverOptions error, wrong number of parameters for iparam.\n");
      exit(EXIT_FAILURE);

    }
    /* dparam */
    nval = fscanf(ficin, "%lf%lf%lf", &(options->dparam[0]), &(options->dparam[1]), &(options->dparam[2]));
    if (nval != 3)
    {
      fprintf(stderr, "Numerics, readSolverOptions error, wrong number of parameters for dparam.\n");
      exit(EXIT_FAILURE);

    }
    fclose(ficin);
    break;
  default:
    fprintf(stderr, "Numerics, readSolverOptions error, unknown problem type.\n");
    exit(EXIT_FAILURE);
  }
}

void recursive_printSolverOptions(SolverOptions* options, int level)
{
  char* marge;
  int i;
  marge = (char*) malloc((level + 1) * sizeof(char));
  for (i = 0; i < level; i++)
    marge[i] = ' ';
  marge[level] = '\0';

  printf("%s\n ========== Numerics Non Smooth Solver parameters: \n", marge);
  if (options->isSet == 0)
    printf("%sThe solver parameters have not been set. \t options->isSet = %i \n", marge, options->isSet);
  else
  {
    printf("%sThe solver parameters below have  been set \t options->isSet = %i\n", marge, options->isSet);
    printf("%sId of the solver\t\t\t\t options->solverId = %d \n", marge, options->solverId);
    printf("%sName of the solver\t\t\t\t  %s \n", marge, idToName(options->solverId));
    if (options->iparam != NULL)
    {
      printf("%sint parameters \t\t\t\t\t options->iparam\n", marge);
      printf("%ssize of the int parameters\t\t\t options->iSize = %i\n", marge, options->iSize);
      for (int i = 0; i < options->iSize; ++i)
        printf("%s\t\t\t\t\t\t options->iparam[%i] = %d\n", marge, i, options->iparam[i]);
    }
    if (options->dparam != NULL)
    {
      printf("%sdouble parameters \t\t\t\t options->dparam\n", marge);
      printf("%ssize of the double parameters\t\t\t options->dSize = %i\n", marge, options->dSize);
      for (int i = 0; i < options->iSize; ++i)
        printf("%s\t\t\t\t\t\t options->dparam[%i] = %.6le\n", marge, i, options->dparam[i]);
    }
  }
  if (options->iWork == NULL)
  {
    printf("%sinteger work array have not been allocated. \t options->iWork = NULL \n", marge);
  }
  else
  {
    printf("%sinteger work array have been allocated. \t options->iWork = %p \n", marge, options->iWork);
    printf("%sinteger work array size \t\t\t options->iSize = %i \n", marge, options->iSize);
  }
  if (options->dWork == NULL)
  {
    printf("%sdouble work array have not been allocated. \t options->dWork = NULL \n", marge);
  }
  else
  {
    printf("%sdouble work array have been allocated. \t options->dWork = %p \n", marge, options->dWork);
    printf("%sdouble work array size \t\t\t options->dSize = %i \n", marge, options->dSize);
  }




  printf("%sSee %s documentation for parameters definition)\n", marge, idToName(options->solverId));

  printf("\n");

  printf("%snumber of internal (or local) solvers \t\t options->numberOfInternalSolvers = %i\n", marge, options->numberOfInternalSolvers);
  for (i = 0; i < options->numberOfInternalSolvers; i++)
  {
    recursive_printSolverOptions(options->internalSolvers + i, level + 1);
  }
  free(marge);

}
void printSolverOptions(SolverOptions* options)
{
  recursive_printSolverOptions(options, 0);
}
void recursive_deleteSolverOptions(SolverOptions* op)
{

  for (int i = 0; i < op->numberOfInternalSolvers; i++)
    recursive_deleteSolverOptions(&(op->internalSolvers[i]));

  if (op->numberOfInternalSolvers && op->internalSolvers)
    free(op->internalSolvers);
  op->internalSolvers = 0;
  if (op->iparam)
    free(op->iparam);
  op->iparam = 0;
  if (op->dparam)
    free(op->dparam);
  op->dparam = 0;
  if (op->iWork)
    free(op->iWork);
  op->iWork = 0;
  if (op->dWork)
    free(op->dWork);
  op->dWork = 0;
}


void deleteSolverOptions(SolverOptions* op)
{
  for (int i = 0; i < op->numberOfInternalSolvers; i++)
    recursive_deleteSolverOptions(&(op->internalSolvers[i]));
  if (op->numberOfInternalSolvers && op->internalSolvers)
    free(op->internalSolvers);
  op->internalSolvers = 0;
  if (op->iparam)
    free(op->iparam);
  op->iparam = 0;
  if (op->dparam)
    free(op->dparam);
  op->dparam = 0;
  if (op->iWork)
    free(op->iWork);
  op->iWork = 0;
  if (op->dWork)
    free(op->dWork);
  op->dWork = 0;



}




char * idToName(int Id)
{
  switch (Id)
  {

    /*MCLP*/
  case    SICONOS_MLCP_PGS:
    return SICONOS_MLCP_PGS_STR;
  case   SICONOS_MLCP_RPGS:
    return SICONOS_MLCP_RPGS_STR;
  case   SICONOS_MLCP_PSOR :
    return SICONOS_MLCP_PSOR_STR;
  case   SICONOS_MLCP_RPSOR :
    return SICONOS_MLCP_RPSOR_STR;
  case   SICONOS_MLCP_PATH :
    return SICONOS_MLCP_PATH_STR;
  case   SICONOS_MLCP_ENUM :
    return SICONOS_MLCP_ENUM_STR;
  case   SICONOS_MLCP_SIMPLEX :
    return SICONOS_MLCP_SIMPLEX_STR;
  case   SICONOS_MLCP_DIRECT_ENUM :
    return SICONOS_MLCP_DIRECT_ENUM_STR;
  case   SICONOS_MLCP_PATH_ENUM :
    return SICONOS_MLCP_PATH_ENUM_STR;
  case   SICONOS_MLCP_DIRECT_SIMPLEX :
    return SICONOS_MLCP_DIRECT_SIMPLEX_STR;
  case   SICONOS_MLCP_DIRECT_PATH :
    return SICONOS_MLCP_DIRECT_PATH_STR;
  case   SICONOS_MLCP_DIRECT_PATH_ENUM :
    return SICONOS_MLCP_DIRECT_PATH_ENUM_STR;
  case   SICONOS_MLCP_FB :
    return SICONOS_MLCP_FB_STR;
  case   SICONOS_MLCP_DIRECT_FB :
    return SICONOS_MLCP_DIRECT_FB_STR;
    /*LCP*/
  case   SICONOS_LCP_LEMKE:
    return SICONOS_LCP_LEMKE_STR;
  case    SICONOS_LCP_NSGS_SBM :
    return SICONOS_LCP_NSGS_SBM_STR;
  case    SICONOS_LCP_PGS:
    return SICONOS_LCP_PGS_STR;
  case    SICONOS_LCP_CPG :
    return SICONOS_LCP_CPG_STR;
  case    SICONOS_LCP_LATIN :
    return SICONOS_LCP_LATIN_STR;
  case    SICONOS_LCP_LATIN_W:
    return SICONOS_LCP_LATIN_W_STR;
  case    SICONOS_LCP_QP:
    return SICONOS_LCP_QP_STR;
  case    SICONOS_LCP_NSQP:
    return SICONOS_LCP_NSQP_STR;
  case    SICONOS_LCP_NEWTONMIN:
    return SICONOS_LCP_NEWTONMIN_STR ;
  case    SICONOS_LCP_NEWTONFB :
    return SICONOS_LCP_NEWTONFB_STR;
  case    SICONOS_LCP_PSOR :
    return SICONOS_LCP_PSOR_STR;
  case    SICONOS_LCP_RPGS :
    return SICONOS_LCP_RPGS_STR;
  case    SICONOS_LCP_PATH:
    return SICONOS_LCP_PATH_STR;
  case     SICONOS_LCP_ENUM :
    return SICONOS_LCP_ENUM_STR;
    /*RELAY*/
  case SICONOS_RELAY_PGS:
    return SICONOS_RELAY_PGS_STR;
  case SICONOS_RELAY_ENUM:
    return SICONOS_RELAY_ENUM_STR;
  case SICONOS_RELAY_PATH:
    return SICONOS_RELAY_PATH_STR;
  case SICONOS_RELAY_LEMKE:
    return SICONOS_RELAY_LEMKE_STR;
  case SICONOS_RELAY_NLGS:
    return SICONOS_RELAY_NLGS_STR;
  case SICONOS_RELAY_LATIN:
    return SICONOS_RELAY_LATIN_STR;
    /*FRICTION_2D*/
  case SICONOS_FRICTION_2D_NSGS:
    return SICONOS_FRICTION_2D_NSGS_STR;
  case SICONOS_FRICTION_2D_NLGS:
    return SICONOS_FRICTION_2D_NLGS_STR;
  case SICONOS_FRICTION_2D_PGS:
    return SICONOS_FRICTION_2D_PGS_STR;
  case SICONOS_FRICTION_2D_CPG:
    return SICONOS_FRICTION_2D_CPG_STR;
  case SICONOS_FRICTION_2D_LATIN:
    return SICONOS_FRICTION_2D_LATIN_STR;
    /*FRICTION_3D*/
  case SICONOS_FRICTION_3D_NSGS:
    return SICONOS_FRICTION_3D_NSGS_STR;
  case SICONOS_FRICTION_3D_NSGSV:
    return SICONOS_FRICTION_3D_NSGSV_STR;
  case SICONOS_FRICTION_3D_PROX:
    return SICONOS_FRICTION_3D_PROX_STR;
  case SICONOS_FRICTION_3D_TFP:
    return SICONOS_FRICTION_3D_TFP_STR;
  case SICONOS_FRICTION_3D_GLOBALAC:
    return SICONOS_FRICTION_3D_GLOBALAC_STR;
  case SICONOS_FRICTION_3D_DSFP:
    return SICONOS_FRICTION_3D_DSFP_STR;
  case SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint:
    return SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint_STR;
  case SICONOS_FRICTION_3D_AlartCurnierNewton:
    return SICONOS_FRICTION_3D_AlartCurnierNewton_STR;
  case SICONOS_FRICTION_3D_DampedAlartCurnierNewton:
    return SICONOS_FRICTION_3D_DampedAlartCurnierNewton_STR;
  case SICONOS_FRICTION_3D_NCPGlockerFBNewton:
    return SICONOS_FRICTION_3D_NCPGlockerFBNewton_STR;
  case SICONOS_FRICTION_3D_ProjectionOnConeWithDiagonalization:
    return SICONOS_FRICTION_3D_ProjectionOnConeWithDiagonalization_STR;
  case SICONOS_FRICTION_3D_ProjectionOnCone:
    return SICONOS_FRICTION_3D_ProjectionOnCone_STR;
  case SICONOS_FRICTION_3D_ProjectionOnConeWithLocalIteration:
    return SICONOS_FRICTION_3D_ProjectionOnConeWithLocalIteration_STR;
  case SICONOS_FRICTION_3D_projectionOnConeWithRegularization:
    return SICONOS_FRICTION_3D_projectionOnConeWithRegularization_STR;
  case SICONOS_FRICTION_3D_NCPGlockerFBPATH:
    return SICONOS_FRICTION_3D_NCPGlockerFBPATH_STR;
  case SICONOS_FRICTION_3D_projectionOnCylinder:
    return SICONOS_FRICTION_3D_projectionOnCylinder_STR;
  case SICONOS_FRICTION_3D_ProjectionOnCone_velocity:
    return SICONOS_FRICTION_3D_ProjectionOnCone_velocity_STR;
  case SICONOS_FRICTION_3D_PGoC:
    return SICONOS_FRICTION_3D_PGoC_STR;
  case SICONOS_FRICTION_3D_DeSaxceFixedPoint:
    return SICONOS_FRICTION_3D_DeSaxceFixedPoint_STR;
  case SICONOS_FRICTION_3D_QUARTIC:
    return SICONOS_FRICTION_3D_QUARTIC_STR;
  case SICONOS_FRICTION_3D_QUARTIC_NU:
    return SICONOS_FRICTION_3D_QUARTIC_NU_STR;
    /*3D_PRIMAL*/
  case SICONOS_FRICTION_3D_PRIMAL_NSGS_WR:
    return SICONOS_FRICTION_3D_PRIMAL_NSGS_WR_STR;
  case SICONOS_FRICTION_3D_PRIMAL_NSGSV_WR:
    return SICONOS_FRICTION_3D_PRIMAL_NSGSV_WR_STR;
  case SICONOS_FRICTION_3D_PRIMAL_PROX_WR:
    return SICONOS_FRICTION_3D_PRIMAL_PROX_WR_STR;
  case SICONOS_FRICTION_3D_PRIMAL_DSFP_WR:
    return SICONOS_FRICTION_3D_PRIMAL_DSFP_WR_STR;
  case SICONOS_FRICTION_3D_PRIMAL_TFP_WR:
    return SICONOS_FRICTION_3D_PRIMAL_TFP_WR_STR;
  case SICONOS_FRICTION_3D_PRIMAL_GLOBALAC_WR:
    return SICONOS_FRICTION_3D_PRIMAL_GLOBALAC_WR_STR;
  case SICONOS_FRICTION_3D_PRIMAL_NSGS:
    return SICONOS_FRICTION_3D_PRIMAL_NSGS_STR;
    /*DEFAULT*/
  default:
    return SICONOS_NONAME_STR;


  }
}
int nameToId(char * pName)
{
  /*MLCP*/
  if (strcmp(SICONOS_MLCP_PGS_STR, pName) == 0)
    return SICONOS_MLCP_PGS;
  else if (strcmp(SICONOS_MLCP_RPGS_STR, pName) == 0)
    return SICONOS_MLCP_RPGS;
  else if (strcmp(SICONOS_MLCP_PSOR_STR, pName) == 0)
    return SICONOS_MLCP_PSOR;
  else if (strcmp(SICONOS_MLCP_RPSOR_STR, pName) == 0)
    return SICONOS_MLCP_RPSOR;
  else if (strcmp(SICONOS_MLCP_PATH_STR, pName) == 0)
    return SICONOS_MLCP_PATH;
  else if (strcmp(SICONOS_MLCP_ENUM_STR, pName) == 0)
    return SICONOS_MLCP_ENUM;
  else if (strcmp(SICONOS_MLCP_SIMPLEX_STR, pName) == 0)
    return SICONOS_MLCP_SIMPLEX;
  else if (strcmp(SICONOS_MLCP_DIRECT_ENUM_STR, pName) == 0)
    return SICONOS_MLCP_DIRECT_ENUM;
  else if (strcmp(SICONOS_MLCP_PATH_ENUM_STR, pName) == 0)
    return SICONOS_MLCP_PATH_ENUM ;
  else if (strcmp(SICONOS_MLCP_DIRECT_SIMPLEX_STR, pName) == 0)
    return SICONOS_MLCP_DIRECT_SIMPLEX;
  else if (strcmp(SICONOS_MLCP_DIRECT_PATH_STR, pName) == 0)
    return SICONOS_MLCP_DIRECT_PATH;
  else if (strcmp(SICONOS_MLCP_FB_STR, pName) == 0)
    return SICONOS_MLCP_FB;
  else if (strcmp(SICONOS_MLCP_DIRECT_FB_STR, pName) == 0)
    return SICONOS_MLCP_DIRECT_FB;
  /*LCP*/
  else if (strcmp(SICONOS_LCP_LEMKE_STR, pName) == 0)
    return SICONOS_LCP_LEMKE;
  else if (strcmp(SICONOS_LCP_NSGS_SBM_STR, pName) == 0)
    return SICONOS_LCP_NSGS_SBM;
  else if (strcmp(SICONOS_LCP_PGS_STR, pName) == 0)
    return SICONOS_LCP_PGS;
  else if (strcmp(SICONOS_LCP_CPG_STR, pName) == 0)
    return SICONOS_LCP_CPG;
  else if (strcmp(SICONOS_LCP_LATIN_STR, pName) == 0)
    return SICONOS_LCP_LATIN;
  else if (strcmp(SICONOS_LCP_LATIN_W_STR, pName) == 0)
    return SICONOS_LCP_LATIN_W;
  else if (strcmp(SICONOS_LCP_QP_STR, pName) == 0)
    return SICONOS_LCP_QP;
  else if (strcmp(SICONOS_LCP_NSQP_STR, pName) == 0)
    return SICONOS_LCP_NSQP;
  else if (strcmp(SICONOS_LCP_NEWTONMIN_STR, pName) == 0)
    return SICONOS_LCP_NEWTONMIN;
  else if (strcmp(SICONOS_LCP_NEWTONFB_STR, pName) == 0)
    return SICONOS_LCP_NEWTONFB;
  else if (strcmp(SICONOS_LCP_PSOR_STR, pName) == 0)
    return SICONOS_LCP_PSOR;
  else if (strcmp(SICONOS_LCP_RPGS_STR, pName) == 0)
    return SICONOS_LCP_RPGS;
  else if (strcmp(SICONOS_LCP_PATH_STR, pName) == 0)
    return SICONOS_LCP_PATH;
  else if (strcmp(SICONOS_LCP_ENUM_STR, pName) == 0)
    return SICONOS_LCP_ENUM;
  /*RELAY*/
  else if (strcmp(SICONOS_RELAY_PGS_STR, pName) == 0)
    return SICONOS_RELAY_PGS;
  else if (strcmp(SICONOS_RELAY_ENUM_STR, pName) == 0)
    return SICONOS_RELAY_ENUM ;
  else if (strcmp(SICONOS_RELAY_PATH_STR, pName) == 0)
    return SICONOS_RELAY_PATH ;
  else if (strcmp(SICONOS_RELAY_LEMKE_STR, pName) == 0)
    return SICONOS_RELAY_LEMKE;
  else if (strcmp(SICONOS_RELAY_NLGS_STR, pName) == 0)
    return SICONOS_RELAY_NLGS;
  /*FRICTION_2D*/
  else if (strcmp(SICONOS_FRICTION_2D_NSGS_STR, pName) == 0)
    return SICONOS_FRICTION_2D_NSGS;
  else if (strcmp(SICONOS_FRICTION_2D_NLGS_STR, pName) == 0)
    return SICONOS_FRICTION_2D_NLGS;
  else if (strcmp(SICONOS_FRICTION_2D_PGS_STR, pName) == 0)
    return SICONOS_FRICTION_2D_PGS;
  else if (strcmp(SICONOS_FRICTION_2D_CPG_STR, pName) == 0)
    return SICONOS_FRICTION_2D_CPG;
  else if (strcmp(SICONOS_FRICTION_2D_LATIN_STR, pName) == 0)
    return SICONOS_FRICTION_2D_LATIN;
  else if (strcmp(SICONOS_RELAY_LATIN_STR, pName) == 0)
    return SICONOS_RELAY_LATIN;
  /*FRICTION_3D*/
  else if (strcmp(SICONOS_FRICTION_3D_NSGS_STR, pName) == 0)
    return SICONOS_FRICTION_3D_NSGS;
  else if (strcmp(SICONOS_FRICTION_3D_NSGSV_STR, pName) == 0)
    return SICONOS_FRICTION_3D_NSGSV;
  else if (strcmp(SICONOS_FRICTION_3D_PROX_STR, pName) == 0)
    return SICONOS_FRICTION_3D_PROX;
  else if (strcmp(SICONOS_FRICTION_3D_TFP_STR, pName) == 0)
    return SICONOS_FRICTION_3D_TFP;
  else if (strcmp(SICONOS_FRICTION_3D_GLOBALAC_STR, pName) == 0)
    return SICONOS_FRICTION_3D_GLOBALAC;
  else if (strcmp(SICONOS_FRICTION_3D_DSFP_STR, pName) == 0)
    return SICONOS_FRICTION_3D_DSFP;
  else if (strcmp(SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint_STR, pName) == 0)
    return SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint;
  else if (strcmp(SICONOS_FRICTION_3D_AlartCurnierNewton_STR, pName) == 0)
    return SICONOS_FRICTION_3D_AlartCurnierNewton;
  else if (strcmp(SICONOS_FRICTION_3D_DampedAlartCurnierNewton_STR, pName) == 0)
    return SICONOS_FRICTION_3D_DampedAlartCurnierNewton;
  else if (strcmp(SICONOS_FRICTION_3D_NCPGlockerFBNewton_STR, pName) == 0)
    return SICONOS_FRICTION_3D_NCPGlockerFBNewton;
  else if (strcmp(SICONOS_FRICTION_3D_ProjectionOnConeWithDiagonalization_STR, pName) == 0)
    return SICONOS_FRICTION_3D_ProjectionOnConeWithDiagonalization;
  else if (strcmp(SICONOS_FRICTION_3D_ProjectionOnCone_STR, pName) == 0)
    return SICONOS_FRICTION_3D_ProjectionOnCone;
  else if (strcmp(SICONOS_FRICTION_3D_ProjectionOnConeWithLocalIteration_STR, pName) == 0)
    return SICONOS_FRICTION_3D_ProjectionOnConeWithLocalIteration;
  else if (strcmp(SICONOS_FRICTION_3D_projectionOnConeWithRegularization_STR, pName) == 0)
    return SICONOS_FRICTION_3D_projectionOnConeWithRegularization;
  else if (strcmp(SICONOS_FRICTION_3D_NCPGlockerFBPATH_STR, pName) == 0)
    return SICONOS_FRICTION_3D_NCPGlockerFBPATH;
  else if (strcmp(SICONOS_FRICTION_3D_projectionOnCylinder_STR, pName) == 0)
    return SICONOS_FRICTION_3D_projectionOnCylinder;
  else if (strcmp(SICONOS_FRICTION_3D_ProjectionOnCone_velocity_STR, pName) == 0)
    return SICONOS_FRICTION_3D_ProjectionOnCone_velocity;
  else if (strcmp(SICONOS_FRICTION_3D_PGoC_STR, pName) == 0)
    return SICONOS_FRICTION_3D_PGoC;
  else if (strcmp(SICONOS_FRICTION_3D_DeSaxceFixedPoint_STR, pName) == 0)
    return SICONOS_FRICTION_3D_DeSaxceFixedPoint;
  else if (strcmp(SICONOS_FRICTION_3D_QUARTIC_STR, pName) == 0)
    return SICONOS_FRICTION_3D_QUARTIC;
  else if (strcmp(SICONOS_FRICTION_3D_QUARTIC_NU_STR, pName) == 0)
    return SICONOS_FRICTION_3D_QUARTIC_NU;
  /*FRICTION_3D_PRIMAL**/
  else if (strcmp(SICONOS_FRICTION_3D_PRIMAL_NSGS_WR_STR, pName) == 0)
    return SICONOS_FRICTION_3D_PRIMAL_NSGS_WR;
  else if (strcmp(SICONOS_FRICTION_3D_PRIMAL_NSGSV_WR_STR, pName) == 0)
    return SICONOS_FRICTION_3D_PRIMAL_NSGSV_WR;
  else if (strcmp(SICONOS_FRICTION_3D_PRIMAL_PROX_WR_STR, pName) == 0)
    return SICONOS_FRICTION_3D_PRIMAL_PROX_WR;
  else if (strcmp(SICONOS_FRICTION_3D_PRIMAL_DSFP_WR_STR, pName) == 0)
    return SICONOS_FRICTION_3D_PRIMAL_DSFP_WR;
  else if (strcmp(SICONOS_FRICTION_3D_PRIMAL_TFP_WR_STR, pName) == 0)
    return SICONOS_FRICTION_3D_PRIMAL_TFP_WR;
  else if (strcmp(SICONOS_FRICTION_3D_PRIMAL_GLOBALAC_WR_STR, pName) == 0)
    return SICONOS_FRICTION_3D_PRIMAL_GLOBALAC_WR;
  else if (strcmp(SICONOS_FRICTION_3D_PRIMAL_NSGS_STR, pName) == 0)
    return SICONOS_FRICTION_3D_PRIMAL_NSGS;
  return 0;

}
