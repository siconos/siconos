/* Siconos-Numerics version 3.0.0, Copyright INRIA 2005-2008.
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

|A C| |u| |a| |0|
|   |*| |+| |=| |
|D B| |v| |b| |w|
0<z*v>0
dim(u)=mm
dim(v)=nn

**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LA.h"
#include "MLCP_Solvers.h"
#include "mlcp_tool.h"
#include <math.h>
#include "mlcp_enum.h"
#include "mlcp_enum_tool.h"



#ifdef MLCP_DEBUG
static int *sLastIWork;
static double *sLastDWork;
#endif

static double * sQ = 0;
static double * sColNul = 0;
static double * sM = 0;
static double * sMref = 0;
static double * sQref = 0;
static int sVerbose = 0;
static int sNn = 0;
static int sMm = 0;
static int* sW2V = 0;

/*OUTPUT */
/*sW2 is a pointer on the output w*/
static double* sW2;
/*sW1 is a pointer on the output w*/
static double* sW1;
/*sV is a pointer on the output v*/
static double* sV;
/*sU is a pointer on the output u*/
static double* sU;


/*case defined with sCurrentEnum
 *if sW2V[i]==0
 *  v[i] not null w2[i] null
 *else
 *  v[i] null and w2[i] not null
 */
/* void affectW2V(){ */
/*   unsigned long  int aux = sCurrentEnum; */
/*   for (unsigned int i =0; i <sMm; i++){ */
/*     sW2V[i]=aux & 1; */
/*     aux = aux >> 1; */
/*   } */

/* } */
/* void initEnum(){ */
/*   int cmp; */
/*   sCmpEnum=0;   */
/*   sNbCase = 1; */
/*   for (cmp=0;cmp<sMm;cmp++) */
/*     sNbCase = sNbCase << 1; */
/*   sProgress=0; */
/* } */
/* int nextEnum(){ */
/*   if (sCmpEnum == sNbCase) */
/*     return 0; */
/*   if (sCurrentEnum >= sNbCase){ */
/*     sCurrentEnum=0; */
/*   } */
/*   if (sVerbose) */
/*     printf("try enum :%d\n",(int)sCurrentEnum); */
/*   affectW2V(); */
/*   sCurrentEnum++; */
/*   sCmpEnum++; */
/*   if (sVerbose && sCmpEnum > sProgress*sNbCase){ */
/*     sProgress+=0.001; */
/*     printf(" progress %f %d \n",sProgress,(int) sCurrentEnum); */
/*   } */

/*   return 1; */
/* } */


void buildQ()
{
  memcpy(sQ, sQref, (sMm + sNn)*sizeof(double));
}
void printCurrentSystem()
{
  int npm = sNn + sMm;
  printf("printCurrentSystemM:\n");
  displayMat(sM, npm, npm, 0);
  printf("printCurrentSystemQ (ie -Q from mlcp beause of linear system MZ=Q):\n");
  displayMat(sQ, sMm + sNn, 1, 0);
}
void printRefSystem()
{
  int npm = sNn + sMm;
  printf("ref M n %d  m %d :\n", sNn, sMm);
  displayMat(sMref, npm, npm, 0);
  printf("ref Q (ie -Q from mlcp beause of linear system MZ=Q):\n");
  displayMat(sQref, sMm + sNn, 1, 0);
}
int mlcp_enum_getNbIWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options)
{
  return 2 * (problem->n + problem->m);
}
int mlcp_enum_getNbDWork(MixedLinearComplementarity_Problem* problem, Solver_Options* options)
{
  return 3 * (problem->n + problem->m) + (problem->n + problem->m) * (problem->n + problem->m);
}

/*
 * The are no memory allocation in mlcp_enum, all necessary memory must be allocated by the user.
 *
 *options:
 * iparam[0] : (in) verbose.
 * dparam[0] : (in) a positive value, tolerane about the sign.
 * dWork : working float zone size : (nn+mm)*(nn+mm) + 3*(nn+mm). MUST BE ALLOCATED BY THE USER.
 * iWork : working int zone size : 2(nn+mm). MUST BE ALLOCATED BY THE USER.
 * double *z : size n+m
 * double *w : size n+m
 */

void mlcp_enum(MixedLinearComplementarity_Problem* problem, double *z, double *w, int *info, Solver_Options* options)
{
  double tol ;
  double * workingFloat = options->dWork;
  int * workingInt = options->iWork;
  int lin;
  int npm = (problem->n) + (problem->m);
  int npm2 = npm * npm;
  int NRHS = 1;
  int one = 1;
  int * ipiv;
  int check;
  int DGESVinfo;

  /*OUTPUT param*/
  sW1 = w;
  sW2 = w + problem->n;
  sU = z;
  sV = z + problem->n;
  tol = options->dparam[0];

  sMref = problem->M->matrix0;

  sNn = problem->n;
  sMm = problem->m;
  sVerbose = options->iparam[0];
  if (sVerbose)
    printf("mlcp_enum begin, n %d m %d tol %lf\n", sNn, sMm, tol);

  sM = workingFloat;
  sQ = sM + npm2;
  sColNul = sQ + sMm + sNn;
  sQref = sColNul + sMm + sNn;
  for (lin = 0; lin < sMm + sNn; lin++)
    sQref[lin] =  - problem->q[lin];
  for (lin = 0; lin < npm; lin++)
    sColNul[lin] = 0;
  /*  printf("sColNul\n");
      displayMat(sColNul,npm,1);*/
  if (sVerbose)
    printRefSystem();
  sW2V = workingInt;
  ipiv = sW2V + sMm;
  *info = 0;
  initEnum(problem->m);
  while (nextEnum(sW2V))
  {
    mlcp_buildM(sW2V, sM, sMref, sNn, sMm);
    buildQ();
    if (sVerbose)
      printCurrentSystem();
    DGESV(npm, NRHS, sM, npm, ipiv, sQ, npm, DGESVinfo);

    if (!DGESVinfo)
    {
      if (sVerbose)
      {
        printf("LU foctorization success, solution in cone?\n");
        displayMat(sQ, npm, 1, 0);
      }

      check = 1;
      for (lin = 0 ; lin < sMm; lin++)
      {
        if (sQ[sNn + lin] < - tol)
        {
          check = 0;
          break;/*out of the cone!*/
        }
      }
      if (!check)
        continue;
      else
      {
        mlcp_fillSolution(sU, sV, sW1, sW2, sNn, sMm, sW2V, sQ);
        if (sVerbose)
        {
          printf("mlcp_enum find a solution!\n");
          mlcp_DisplaySolution(sU, sV, sW1, sW2, sNn, sMm);
        }
        //  options->iparam[1]=sCurrentEnum-1;
        return;
      }
    }
    else
    {
      if (sVerbose)
      {
        printf("LU foctorization failed:\n");
      }
    }
  }
  *info = 1;
  if (sVerbose)
    printf("mlcp_enum failed!\n");
}
