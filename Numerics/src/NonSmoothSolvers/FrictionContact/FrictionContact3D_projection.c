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

#include "LA.h"
#include "FrictionContact3D_Solvers.h"
#include "projectionOnCone.h"
#include "projectionOnCylinder.h"


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

/* Static variables */

/* The global problem of size n= 3*nc, nc being the number of contacts, is locally saved in MGlobal and qGlobal */
/* mu corresponds to the vector of friction coefficients */
/* note that either MGlobal or MBGlobal is used, depending on the chosen storage */
/* static int n=0; */
/* static const NumericsMatrix* MGlobal = NULL; */
/* static const double* qGlobal = NULL; */
/* static const double* mu = NULL; */

/* Local problem operators */
/* static const int nLocal = 3; */
/* static double* MLocal; */
/* static int isMAllocatedIn = 0; /\* True if a malloc is done for MLocal, else false *\/ */
/* static double qLocal[3]; */
/* static double mu_i = 0.0; */

#define VERBOSE_DEBUG

void frictionContact3D_projection_initialize(FrictionContactProblem * problem, FrictionContactProblem * localproblem)

{

}

void frictionContact3D_projection_update(int contact, FrictionContactProblem* problem, FrictionContactProblem* localproblem, double* reaction, SolverOptions* options)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */

  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  frictionContact3D_nsgs_fillMLocal(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/
  frictionContact3D_nsgs_computeqLocal(problem, localproblem, reaction, contact);

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];

}

void frictionContact3D_projectionWithDiagonalization_update(int contact, FrictionContactProblem* problem, FrictionContactProblem* localproblem,  double* reaction, SolverOptions* options)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */

  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  frictionContact3D_nsgs_fillMLocal(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/

  NumericsMatrix * MGlobal = problem->M;
  double * MLocal =  localproblem->M->matrix0;


  double *qLocal = localproblem->q;
  double * qGlobal = problem->q;
  int n = 3 * problem->numberOfContacts;


  int in = 3 * contact, it = in + 1, is = it + 1;
  /* reaction current block set to zero, to exclude current contact block */
  /*   double rin= reaction[in] ; double rit= reaction[it] ; double ris= reaction[is] ;  */
  /* qLocal computation*/
  qLocal[0] = qGlobal[in];
  qLocal[1] = qGlobal[it];
  qLocal[2] = qGlobal[is];

  if (MGlobal->storageType == 0)
  {
    double * MM = MGlobal->matrix0;
    int incx = n, incy = 1;
    qLocal[0] += DDOT(n , &MM[in] , incx , reaction , incy);
    qLocal[1] += DDOT(n , &MM[it] , incx , reaction , incy);
    qLocal[2] += DDOT(n , &MM[is] , incx , reaction , incy);
    // Substract diagonal term
    qLocal[0] -= MM[in + n * in] * reaction[in];
    qLocal[1] -= MM[it + n * it] * reaction[it];
    qLocal[2] -= MM[is + n * is] * reaction[is];
  }
  else if (MGlobal->storageType == 1)
  {
    /* qLocal += rowMB * reaction
       with rowMB the row of blocks of MGlobal which corresponds to the current contact
    */
    subRowProdSBM(n, 3, contact, MGlobal->matrix1, reaction, qLocal, 0);
    // Substract diagonal term
    qLocal[0] -= MLocal[0] * reaction[in];
    qLocal[1] -= MLocal[4] * reaction[it];
    qLocal[2] -= MLocal[8] * reaction[is];

  }
  /*   reaction[in] = rin; reaction[it] = rit; reaction[is] = ris; */

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];
}


void frictionContact3D_projection_initialize_with_regularization(FrictionContactProblem * problem, FrictionContactProblem * localproblem)
{
  if (!localproblem->M->matrix0)
    localproblem->M->matrix0 = (double*)malloc(9 * sizeof(double));
}

void frictionContact3D_projection_update_with_regularization(int contact, FrictionContactProblem * problem, FrictionContactProblem * localproblem, double* reaction, SolverOptions* options)
{


  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */

  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */

  NumericsMatrix * MGlobal = problem->M;

  int n = 3 * problem->numberOfContacts;


  int storageType = MGlobal->storageType;
  if (storageType == 0)
    // Dense storage
  {
    int in = 3 * contact, it = in + 1, is = it + 1;
    int inc = n * in;
    double * MM = MGlobal->matrix0;
    double * MLocal =  localproblem->M->matrix0;

    /* The part of MM which corresponds to the current block is copied into MLocal */
    MLocal[0] = MM[inc + in];
    MLocal[1] = MM[inc + it];
    MLocal[2] = MM[inc + is];
    inc += n;
    MLocal[3] = MM[inc + in];
    MLocal[4] = MM[inc + it];
    MLocal[5] = MM[inc + is];
    inc += n;
    MLocal[6] = MM[inc + in];
    MLocal[7] = MM[inc + it];
    MLocal[8] = MM[inc + is];
  }
  else if (storageType == 1)
  {
    int diagPos = getDiagonalBlockPos(MGlobal->matrix1, contact);
    /*     for (int i =0 ; i< 3*3 ; i++) localproblem->M->matrix0[i] = MGlobal->matrix1->block[diagPos][i] ; */
    DCOPY(9, MGlobal->matrix1->block[diagPos], 1, localproblem->M->matrix0 , 1);

  }
  else
    numericsError("FrictionContact3D_projection -", "unknown storage type for matrix M");

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/
  frictionContact3D_nsgs_computeqLocal(problem, localproblem, reaction, contact);

  double rho = options->dparam[3];
  for (int i = 0 ; i < 3 ; i++) localproblem->M->matrix0[i + 3 * i] += rho ;

  double *qLocal = localproblem->q;
  int in = 3 * contact, it = in + 1, is = it + 1;

  /* qLocal computation*/
  qLocal[0] -= rho * reaction[in];
  qLocal[1] -= rho * reaction[it];
  qLocal[2] -= rho * reaction[is];

  /* Friction coefficient for current block*/
  localproblem->mu[0] = problem->mu[contact];


}

int frictionContact3D_projectionWithDiagonalization_solve(FrictionContactProblem* localproblem, double* reaction, SolverOptions * options)
{



  /* Current block position */

  /* Builds local problem for the current contact */
  /*  frictionContact3D_projection_update(contact, reaction); */
  /*  frictionContact3D_projectionWithDiagonalization_update(contact, reaction);  */


  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  double mu_i = localproblem->mu[0];
  int nLocal = 3;

  double mrn, num, mu2 = mu_i * mu_i;


  /* projection */
  if (qLocal[0] > 0.)
  {
    reaction[0] = 0.;
    reaction[1] = 0.;
    reaction[2] = 0.;
  }
  else
  {
    if (MLocal[0] < DBL_EPSILON || MLocal[nLocal + 1] < DBL_EPSILON || MLocal[2 * nLocal + 2] < DBL_EPSILON)
    {
      fprintf(stderr, "FrictionContact3D_projection error: null term on MLocal diagonal.\n");
      exit(EXIT_FAILURE);
    }

    reaction[0] = -qLocal[0] / MLocal[0];
    reaction[1] = -qLocal[1] / MLocal[nLocal + 1];
    reaction[2] = -qLocal[2] / MLocal[2 * nLocal + 2];

    mrn = reaction[1] * reaction[1] + reaction[2] * reaction[2];

    if (mrn > mu2 * reaction[0]*reaction[0])
    {
      num = mu_i * reaction[0] / sqrt(mrn);
      reaction[1] = reaction[1] * num;
      reaction[2] = reaction[2] * num;
    }
  }
  return 0;
}


int frictionContact3D_projectionOnConeWithLocalIteration_solve(FrictionContactProblem* localproblem, double* reaction, SolverOptions* options)
{
  /* int and double parameters */
  int* iparam = options->iparam;
  double* dparam = options->dparam;

  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  double mu_i = localproblem->mu[0];
  int nLocal = 3;


  /*   /\* Builds local problem for the current contact *\/ */
  /*   frictionContact3D_projection_update(localproblem, reaction); */


  /*double an = 1./(MLocal[0]);*/
  /*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; */
  /*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] + MLocal[1*nLocal+2]; */
  /*   double beta = alpha*alpha - 4*det; */
  /*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */

  double an = 1. / (MLocal[0]);

  double at = 1.0 / (MLocal[4] + mu_i);
  double as = 1.0 / (MLocal[8] + mu_i);
  at = an;
  as = an;
  int incx = 1, incy = 1;



  double worktmp[3];
  double normUT;
  double localerror = 1.0;
  //printf ("localerror = %14.7e\n",localerror );
  int localiter = 0;
  double localtolerance = dparam[0];

  /*     printf ("localtolerance = %14.7e\n",localtolerance ); */
  while ((localerror > localtolerance) && (localiter < iparam[0]))
  {
    localiter ++;

    /*    printf ("reaction[0] = %14.7e\n",reaction[0]); */
    /*    printf ("reaction[1] = %14.7e\n",reaction[1]); */
    /*    printf ("reaction[2] = %14.7e\n",reaction[2]); */
    DCOPY(nLocal , qLocal, incx , worktmp , incy);
    DGEMV(LA_NOTRANS, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, worktmp, incy);
    normUT = sqrt(worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2]);
    /*    printf ("qLocal[0] = %14.7e\n",qLocal[0]); */
    /*    printf ("qLocal[1] = %14.7e\n",qLocal[1]); */
    /*    printf ("qLocal[2] = %14.7e\n",qLocal[2]); */
    /*    printf ("MLocal[0] = %14.7e\n",MLocal[0]); */
    /*    printf ("MLocal[1] = %14.7e\n",MLocal[1]); */
    /*    printf ("MLocal[2] = %14.7e\n",MLocal[2]); */
    /*    printf ("MLocal[3] = %14.7e\n",MLocal[3]); */
    /*    printf ("MLocal[4] = %14.7e\n",MLocal[4]); */
    /*    printf ("MLocal[5] = %14.7e\n",MLocal[5]); */
    /*    printf ("MLocal[6] = %14.7e\n",MLocal[6]); */
    /*    printf ("MLocal[7] = %14.7e\n",MLocal[7]); */
    /*    printf ("MLocal[8] = %14.7e\n",MLocal[8]); */
    /*    printf ("velocity[0] = %14.7e\n",worktmp[0]); */
    /*    printf ("velocity[1] = %14.7e\n",worktmp[1]); */
    /*    printf ("velocity[2] = %14.7e\n",worktmp[2]); */
    /*    printf ("normUT = %14.7e\n",normUT); */
    /*    printf ("Modified velocity[0] = %14.7e\n",worktmp[0]+mu_i*normUT); */
    /*    printf ("Modified velocity[1] = %14.7e\n",worktmp[1]); */
    /*    printf ("Modified velocity[2] = %14.7e\n",worktmp[2]); */



    reaction[0] -= an * (worktmp[0] + mu_i * normUT);
    reaction[1] -= at * worktmp[1];
    reaction[2] -= as * worktmp[2];
    worktmp[0] = reaction[0] + an * (worktmp[0] + mu_i * normUT);
    worktmp[1] = reaction[1] + at * worktmp[1];
    worktmp[2] = reaction[2] + as * worktmp[2];


    projectionOnCone(&reaction[0], mu_i);

    localerror =  sqrt((worktmp[0] - reaction[0]) * (worktmp[0] - reaction[0]) +
                       (worktmp[1] - reaction[1]) * (worktmp[1] - reaction[1]) +
                       (worktmp[2] - reaction[2]) * (worktmp[2] - reaction[2]));
    /*     printf ("localerror = %14.7e\n",localerror );  */
    /*    printf ("worktmp[0] = %14.7e\n",worktmp[0]); */
    /*    printf ("worktmp[1] = %14.7e\n",worktmp[1]); */
    /*    printf ("worktmp[2] = %14.7e\n",worktmp[2]); */
    /*           if (verbose>1) */
    /*             { */
    /*               printf ("reaction[0] = %14.7e\n",reaction[0]);  */
    /*               printf ("reaction[1] = %14.7e\n",reaction[1]);  */
    /*               printf ("reaction[2] = %14.7e\n",reaction[2]);  */
    /*             } */

    if (verbose > 1)
    {
      printf("----------------------  localiter = %i\t error = %.10e \n", localiter, localerror);

    }
    /*    DCOPY( nLocal , qLocal, incx , worktmp , incy ); */
    /*    DGEMV(LA_NOTRANS, nLocal, nLocal, 1.0, MLocal, 3, &reaction[0], incx, 1.0, worktmp, incy ); */
    /*    normUT = sqrt(worktmp[1]*worktmp[1]+worktmp[2]*worktmp[2]); */
    /*    worktmp[0] = reaction[0] - (worktmp[0]+mu_i*normUT); */
    /*    worktmp[1] = reaction[1] - worktmp[1]; */
    /*    worktmp[2] = reaction[2] - worktmp[2]; */

    /*    projectionOnCone(worktmp, mu_i);   */
    /*    localerror2 =  sqrt( (worktmp[0]-reaction[0])*(worktmp[0]-reaction[0]) + */
    /*          (worktmp[1]-reaction[1])*(worktmp[1]-reaction[1])+ */
    /*          (worktmp[2]-reaction[2])*(worktmp[2]-reaction[2]) ); */
    /*     printf ("localerror2 = %14.7e\n",localerror2 );  */


    /* /\*    printf ("velocity[0] = %14.7e\n",worktmp[0]); *\/ */
    /* /\*    printf ("velocity[1] = %14.7e\n",worktmp[1]); *\/ */
    /* /\*    printf ("velocity[2] = %14.7e\n",worktmp[2]); *\/ */
    /* /\*    printf ("normUT = %14.7e\n",normUT); *\/ */
    /* /\*    printf ("Modified velocity[0] = %14.7e\n",worktmp[0]+mu_i*normUT); *\/ */
    /* /\*    printf ("Modified velocity[1] = %14.7e\n",worktmp[1]); *\/ */
    /* /\*    printf ("Modified velocity[2] = %14.7e\n",worktmp[2]); *\/ */
  }
  if (localerror > localtolerance)
    return 1;
  return 0;

}
void frictionContact3D_projectionOnCylinder_update(int contact, FrictionContactProblem* problem, FrictionContactProblem* localproblem, double* reaction, SolverOptions* options)
{
  /* Build a local problem for a specific contact
     reaction corresponds to the global vector (size n) of the global problem.
  */

  /* Call the update function which depends on the storage for MGlobal/MBGlobal */
  /* Build a local problem for a specific contact
   reaction corresponds to the global vector (size n) of the global problem.
  */

  /* The part of MGlobal which corresponds to the current block is copied into MLocal */
  frictionContact3D_nsgs_fillMLocal(problem, localproblem, contact);

  /****  Computation of qLocal = qBlock + sum over a row of blocks in MGlobal of the products MLocal.reactionBlock,
     excluding the block corresponding to the current contact. ****/
  frictionContact3D_nsgs_computeqLocal(problem, localproblem, reaction, contact);

  /* Friction coefficient for current block*/
  localproblem->mu[0] = (options->dWork[contact]);

}


int frictionContact3D_projectionOnCone_solve(FrictionContactProblem* localproblem, double* reaction, SolverOptions * options)
{


  /*  /\* Builds local problem for the current contact *\/ */
  /*   frictionContact3D_projection_update(contact, reaction); */



  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  double mu_i = localproblem->mu[0];
  int nLocal = 3;

  /*double an = 1./(MLocal[0]);*/
  /*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; */
  /*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] + MLocal[1*nLocal+2]; */
  /*   double beta = alpha*alpha - 4*det; */
  /*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */

  //double an = 1./(MLocal[0]+mu_i);
  double an = 1. / (MLocal[0]);
  int incx = 1, incy = 1;
  double worktmp[3];
  double normUT;
  DCOPY(nLocal , qLocal, incx , worktmp , incy);

  DGEMV(LA_NOTRANS, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, worktmp, incy);
  normUT = sqrt(worktmp[1] * worktmp[1] + worktmp[2] * worktmp[2]);
  reaction[0] -= an * (worktmp[0] + mu_i * normUT);
  reaction[1] -= an * worktmp[1];
  reaction[2] -= an * worktmp[2];


  projectionOnCone(reaction, mu_i);
  return 0;

}



void frictionContact3D_projection_free(FrictionContactProblem* localproblem)
{
  free(localproblem->M->matrix0);
  localproblem->M->matrix0 = NULL;
}

void frictionContact3D_projection_with_regularization_free(FrictionContactProblem* localproblem)
{
  free(localproblem->M->matrix0);
  localproblem->M->matrix0 = NULL;
}



int frictionContact3D_projectionOnCone_velocity_solve(FrictionContactProblem* localproblem, double* velocity,  SolverOptions* options)
{
  /* int and double parameters */
  /*     int* iparam = options->iparam; */
  /*     double* dparam = options->dparam; */
  /* Current block position */

  /* Builds local problem for the current contact */
  /*   frictionContact3D_projection_update(contact, velocity); */

  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  double mu_i = localproblem->mu[0];
  int nLocal = 3;


  /*double an = 1./(MLocal[0]);*/
  /*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; */
  /*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] + MLocal[1*nLocal+2]; */
  /*   double beta = alpha*alpha - 4*det; */
  /*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */

  double an = 1. / (MLocal[0]);
  int incx = 1, incy = 1;
  double worktmp[3];
  double normUT;

  DCOPY(nLocal , qLocal, incx , worktmp , incy);
  DGEMV(LA_NOTRANS, nLocal, nLocal, 1.0, MLocal, 3, velocity, incx, 1.0, worktmp, incy);
  normUT = sqrt(velocity[1] * velocity[1] + velocity[2] * velocity[2]);
  velocity[0] -=  - mu_i * normUT + an * (worktmp[0]);
  velocity[1] -= an * worktmp[1];
  velocity[2] -= an * worktmp[2];
  double invmui = 1.0 / mu_i;
  projectionOnCone(velocity, invmui);

  normUT = sqrt(velocity[1] * velocity[1] + velocity[2] * velocity[2]);
  velocity[0] -= mu_i * normUT;
  return 0;
}



int frictionContact3D_projectionOnCylinder_solve(FrictionContactProblem *localproblem , double* reaction, SolverOptions* options)
{
  /* int and double parameters */
  /*   int* iparam = options->iparam; */
  /*   double* dparam = options->dparam; */
  double * MLocal = localproblem->M->matrix0;
  double * qLocal = localproblem->q;
  int nLocal = 3;


  /* Builds local problem for the current contact */
  /*   frictionContact3D_projection_update(contact, reaction); */


  /*double an = 1./(MLocal[0]);*/
  /*   double alpha = MLocal[nLocal+1] + MLocal[2*nLocal+2]; */
  /*   double det = MLocal[1*nLocal+1]*MLocal[2*nLocal+2] - MLocal[2*nLocal+1] + MLocal[1*nLocal+2]; */
  /*   double beta = alpha*alpha - 4*det; */
  /*   double at = 2*(alpha - beta)/((alpha + beta)*(alpha + beta)); */

  double an = 1. / (MLocal[0]);

  int incx = 1, incy = 1;
  double worktmp[3];

  double R  = localproblem->mu[0];
  DCOPY(nLocal , qLocal, incx , worktmp , incy);
  DGEMV(LA_NOTRANS, nLocal, nLocal, 1.0, MLocal, 3, reaction, incx, 1.0, worktmp, incy);
  reaction[0]   -= an * worktmp[0];
  reaction[1] -= an * worktmp[1];
  reaction[2] -= an * worktmp[2];

  projectionOnCylinder(reaction, R);
  return 0;

}

