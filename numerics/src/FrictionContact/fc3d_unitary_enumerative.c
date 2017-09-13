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

#define _XOPEN_SOURCE 700


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "fc3d_Solvers.h"
#include "Friction_cst.h"
#include "op3x3.h"
#include "fc3d_unitary_enumerative.h"
//#define FC3D_UE_DEBUG

#include "quartic.h"
#include "projectionOnCone.h"
#include "fc3d_compute_error.h"

#define FC3D_UE_TEST_NULL(EXPR)  (fabs(EXPR)<1e-15)

#pragma GCC diagnostic ignored "-Wmissing-prototypes"


static void solve2x2(double *a, double *b, double *c, double *a1, double *b1, double *c1, double *x, double *y);


void compute_racines(double * Poly, int *nbRealRacines, double *Racines)
{
  double r[3][5];
  //Roots of poly p[0] x^4 + p[1] x^3...+p[4]=0
  //x=r[1][k] + i r[2][k]  k=1,...,4
#ifdef FC3D_UE_DEBUG
  double Psav[5];
  printf("compute_racines: polynome(x)=%e.x4+%e.x3+%e.x2+%e.x+%e\n", Poly[0], Poly[1], Poly[2], Poly[3], Poly[4]);
  for (int k = 0; k < 5; k++)
    Psav[k] = Poly[k];
#endif
  if (fabs(Poly[1] / Poly[0]) > 1e7)
  {
    Poly[0] = 0.0;
#ifdef FC3D_UE_DEBUG
    printf("compute_racines: WARNING, Poly[1]/Poly[0] =%e. set Poly[0]=0\n", fabs(Poly[1] / Poly[0]));
#endif
  }
  int degp1 = 5;
  if (Poly[0] != 0.0)
    BIQUADROOTS(Poly, r);
  else if (Poly[1] != 0.0)
  {
    CUBICROOTS(Poly + 1, r);
    degp1 = 4;
  }
  else if (Poly[2] != 0.0)
  {
    QUADROOTS(Poly + 2, r);
    degp1 = 3;
  }
  else if (Poly[3] != 0.0)
  {
    r[1][1] = -Poly[0] / Poly[3];
    r[2][1] = 0;
    degp1 = 2;
  }
  else
  {
    if (verbose) printf("FC3D_unitary_enumerative: degre of polynom is 0.");
    degp1 = 1;
  }
  (*nbRealRacines) = 0;
  for (int k = 1; k < degp1; k++)
  {
    if (fabs(r[2][k]) < 1e-10)
    {
      Racines[*nbRealRacines] = r[1][k];
      (*nbRealRacines)++;
    }
  }
#ifdef FC3D_UE_DEBUG
  for (int k = 0; k < *nbRealRacines; k++)
  {
    printf("compute_racines debug : Psav(Racines[%d]=%e)=%e\n", k, Racines[k],
           Psav[0]*Racines[k]*Racines[k]*Racines[k]*Racines[k] +
           Psav[1]*Racines[k]*Racines[k]*Racines[k] +
           Psav[2]*Racines[k]*Racines[k] +
           Psav[3]*Racines[k] +
           Psav[4]);
  }
#endif
}

/*

M=|ac|
  |cb|
  V is such that V*M*V'=|l1 0|
                        |0 l2|
  V'=|V1_x V2_x|
     |V1_y V2_y|
     so
  V =|V1_x V1_y|
     |V2_x V2_y|

 */
void FC3D_unitary_enum_factorize2x2(double *a, double *b, double *c, double *l1, double *l2, double *V)
{
#ifdef FC3D_UE_DEBUG
  printf("FC3D_unitary_enum_factorize2x2 matrix:\n %e %e \n %e %e \n", *a, *c, *c, *b);
#endif
  if (FC3D_UE_TEST_NULL(*c))
  {
    V[0] = 1;
    V[1] = 0;
    V[2] = 0;
    V[3] = 1;
    *l1 = *a;
    *l2 = *b;
    return;
  }
  double sqrtD = sqrt((*a - *b) * (*a - *b) + 4 * (*c) * (*c));
  *l1 = 0.5 * (*a + *b + sqrtD);
  *l2 = 0.5 * (*a + *b - sqrtD);
  double buff = 0.5 * (*b - *a + sqrtD) / (*c);
  double InvNormV = 1.0 / sqrt(1 + buff * buff);

  V[0] = InvNormV;
  V[2] = InvNormV * (*l1 - *a) / (*c);

  V[1] = InvNormV * (*l2 - *b) / (*c);
  V[3] = InvNormV;


#ifdef FC3D_UE_DEBUG
  double VM[4];
  double M[4];
  M[0] = *a;
  M[1] = *c;
  M[2] = *c;
  M[3] = *b;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
    {
      VM[i + 2 * j] = V[i] * M[2 * j] + V[i + 2] * M[2 * j + 1];
    }
  double Vt[4];
  Vt[0] = V[0];
  Vt[1] = V[2];
  Vt[3] = V[3];
  Vt[2] = V[1];
  double VM_Vt[4];
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
    {
      VM_Vt[i + 2 * j] = VM[i] * Vt[2 * j] + VM[i + 2] * Vt[2 * j + 1];
    }
  printf("FC3D_unitary_enum_factorize2x2 Debug: Following matrix must be diag : d11=l1=%e=%e; d21=%e ;d12=%e ;d22=l2=%e=%e\n", *l1, VM_Vt[0], VM_Vt[1], VM_Vt[2], *l2, VM_Vt[3]);
#endif
}
void fc3d_unitary_enumerative_free(FrictionContactProblem* localproblem)
{
  free(localproblem->M->matrix0);
  localproblem->M->matrix0 = NULL;
}
void fc3d_unitary_enumerative_initialize(FrictionContactProblem* localproblem)
{
  if (!localproblem->M->matrix0)
    localproblem->M->matrix0 = (double*)malloc(9 * sizeof(double));

}


int fc3d_unitary_enumerative_test_non_sliding(FrictionContactProblem* problem, double * reaction, double * velocity, SolverOptions* options)
{
  double * M = problem->M->matrix0;
  double * Q = problem->q;
  double * mu = problem->mu;
  double tol = options->dparam[0];
  SET3X3(M);
  M = M00;
  SET3(reaction);
  reaction = reaction0;
  SET3(velocity);
  velocity = velocity0;
  SET3(Q);
  Q = Q0;

  //printf("fc3d_unitary_enumerative M:\n");
  //print3x3(M);M=M00;
  //printf("fc3d_unitary_enumerative Q:\n");
  //print3(Q);Q=Q0;

  /*take off? R=0 ?*/
  if (*Q0 + tol > 0)
  {
    *reaction0 = 0.;
    *reaction1 = 0.;
    *reaction2 = 0.;
    *velocity0 = *Q0;
    *velocity1 = *Q1;
    *velocity2 = *Q2;
#ifdef FC3D_UE_DEBUG
    printf("FC3D: take off with q=v=(%e,%e,%e)\n", *Q0, *Q1, *Q2);
#endif
    return 0;
  }

  if (*mu == 0.)
  {
    *reaction0 = -(*Q0) / (*M00);
    if (*reaction0 < 0)
    {
      *reaction0 = 0;
      *reaction1 = 0;
      *reaction2 = 0;
      *velocity0 = 0;
      *velocity1 = *Q1;
      *velocity2 = *Q2;
      return -1;
    }
    *reaction1 = 0;
    *reaction2 = 0;
    *velocity0 = 0;
    *velocity1 = *reaction0 * *M10 + *Q1;
    *velocity2 = *reaction0 * *M20 + *Q2;

    return 0;
  }

  /*sticking ? 0=MR+q*/
  //int info = solv3x3(M, reaction, Q);
  reaction[0] = Q[0];
  reaction[1] = Q[1];
  reaction[2] = Q[2];
  int info = solve_3x3_gepp(M, reaction);
  if(info && (verbose > 0))
    numerics_warning("fc3d_unitary_enumerative_test_non_sliding", "NaN output in solve_3x3_gepp");
  M = M00;
  reaction = reaction0;
  Q = Q0;
  if (!isnan(*reaction0))
  {
    if (- *reaction0 + tol > 0.0)
      if ((*reaction0 * *mu) * (*reaction0 * *mu) + tol * tol > *reaction1 * *reaction1 + *reaction2 * *reaction2)
      {
        *reaction0 = - *reaction0;
        *reaction1 = - *reaction1;
        *reaction2 = - *reaction2;
        *velocity0 = 0;
        *velocity1 = 0;
        *velocity2 = 0;
#ifdef FC3D_UE_DEBUG
        printf("FC3D: (Rn*mu)^2+tol=%e. (Rt)^2=%e\n", (*reaction0 * *mu) * (*reaction0 * *mu) + tol, *reaction1 * *reaction1 + *reaction2 * *reaction2);
        printf("FC3D:  0=MR+q with R in cone r=(%e,%e,%e), mu=%e,tol=%e\n", *reaction0, *reaction1, *reaction2, *mu, tol);
#endif
        return 0;
      }
  }
  return 1;

}
/*API for the nsgs*/
int fc3d_unitary_enumerative_solve(FrictionContactProblem* problem, double * reaction,   SolverOptions* options)
{
  int info;
  double velocity[3];
  return fc3d_unitary_enumerative(problem, reaction, velocity, &info, options);
}
int fc3d_unitary_enumerative(FrictionContactProblem* problem, double * reaction, double * velocity, int *info, SolverOptions* options)
{
  *info = fc3d_unitary_enumerative_test_non_sliding(problem, reaction, velocity, options);
  if (!(*info))
    return *info ;
#ifdef FC3D_UE_DEBUG
  printf("FC3D: Not a trivial case, sliding?\n");
#endif
  if (options->solverId == SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU)
  {
    *info = fc3d_unitary_enumerative_solve_poly_nu_sliding(problem, reaction, options);
  }
  else
  {
    *info = fc3d_unitary_enumerative_solve_sliding(problem, reaction, options);
  }
  if (!(*info))
  {
    double * M = problem->M->matrix0;
    double * Q = problem->q;
    SET3(Q);
    Q = Q0;
    SET3(reaction);
    reaction = reaction0;
    SET3(velocity);
    velocity = velocity0;
    mv3x3(M, reaction, velocity);
    reaction = reaction0;
    velocity = velocity0;
    *velocity0 += *Q0;
    *velocity1 += *Q1;
    *velocity2 += *Q2;
  }
#ifdef FC3D_UE_DEBUG
  if (*info)
    printf("fc3d_unitary_enumerative FAILED\n");

  double err = 0.0;
  fc3d_unitary_compute_and_add_error(reaction, velocity, *(problem->mu), &err);
  printf("error is %e.", err);

  printf("fc3d_unitary_enumerative_solve: end\n");
#endif
  return *info ;
}
int fc3d_unitary_enumerative_solve_sliding(FrictionContactProblem* problem, double * reaction,   SolverOptions* options)
{


  double * M = problem->M->matrix0;
  double * Q = problem->q;
  double * mu = problem->mu;
  double D1, D2;
  double alpha;
  double tol = options->dparam[0];
  SET3X3(M);
  M = M00;
  SET3(reaction);
  reaction = reaction0;
  SET3(Q);
  Q = Q0;
#ifdef FC3D_UE_DEBUG
  printf("fc3d_unitary_enumerative (begin) M:\n");
  print3x3(M);
  M = M00;
  printf("fc3d_unitary_enumerative (begin) Q:\n");
  print3(Q);
  Q = Q0;

#endif

  //printf("fc3d_unitary_enumerative M:\n");
  //print3x3(M);M=M00;
  //printf("fc3d_unitary_enumerative Q:\n");
  //print3(Q);Q=Q0;

  /*Sliding? */
#ifdef FC3D_UE_DEBUG
  if (*mu < 10e-20)
  {
    printf("FC3D_UE_DEBUG : wrong value of mu\n");
    return -1;
  }

#endif
  double *Q_2 = Q + 1;
  double V[4];
  double * V00 = V, * V10 = V00 + 1, *V01 = V10 + 1, *V11 = V01 + 1;
  /*D is the projection of the origine on the directrice of the conic (R_T1,R_T2)*/
  double OD[2];
  /*D2 is the projection of the origine on the directrice of the conic (V*R_T1,V*R_T2)*/
  double OD2[2];
  double Q2b[2];
  double M1_b[2];
  double NormD2 = (*M01) * (*M01) + (*M02) * (*M02);
  double NormD = sqrt(NormD2);
  double e = (*mu) * NormD / (*M00);
  double d, p;
  double phi = NAN;
  double cosphi = NAN;
  double sinphi = NAN;
  FC3D_unitary_enum_factorize2x2(M11, M22, M12, &D1, &D2, V);

#ifdef FC3D_UE_DEBUG
  double D_Dir[2];
#endif

  if (!FC3D_UE_TEST_NULL(e))
  {
    d = fabs(*Q0) / NormD;
    p = e * d;
    OD[0] = -((*M01) * (*Q)) / NormD2;


    if (!FC3D_UE_TEST_NULL(*M02))
    {
      OD[1] = -(((-*M01) * (*M01) * (*Q0)) / NormD2 + *Q0) / (*M02);

#ifdef FC3D_UE_DEBUG
      D_Dir[0] = 1;
      D_Dir[1] = -(*M01) / (*M02);
#endif
    }
    else
    {
      OD[1] = 0;
#ifdef FC3D_UE_DEBUG
      D_Dir[0] = 0;
      if (*M01 > 0)
        D_Dir[1] = 1;
      else
        D_Dir[1] = -1;
#endif
    }
    OD2[0] = (*V00) * OD[0] + (*V01) * OD[1];
    OD2[1] = (*V10) * OD[0] + (*V11) * OD[1];
    if (!FC3D_UE_TEST_NULL(OD2[0]))
    {
      if (FC3D_UE_TEST_NULL(OD2[1]))
      {
        phi = 0;
      }
      else
      {
        phi = atan(OD2[1] / OD2[0]);
      }
    }
    else
    {
      phi = M_PI / 2.0;
      if (OD2[1] < 0)
        phi = -M_PI / 2.0;
    }
    if (OD2[0] < 0)
      phi += M_PI;
    cosphi = cos(phi);
    sinphi = sin(phi);
  }
  else  /*e is null*/
  {
    d = 0;
    p = (-(*Q0) * (*mu)) / (*M00);
    if (p < 0)
    {
      projectionOnCone(reaction, *mu);
      return -1;
    }
  }
  double RTb[2];

  /*Q2b is V* \tilde q in the siconos dev note.*/
  Q2b[0] = (*V00) * Q_2[0] + (*V01) * Q_2[1];
  Q2b[1] = (*V10) * Q_2[0] + (*V11) * Q_2[1];
  /*M1_b is V* \tilde M_{1.} in the siconos dev note.*/
  M1_b[0] = (*V00) * (*M10) + (*V01) * (*M20);
  M1_b[1] = (*V10) * (*M10) + (*V11) * (*M20);
  double a1 = -M1_b[0] / (*mu);
  double a2 = -M1_b[1] / (*mu);

  double AA = -e * Q2b[1] * cosphi;
  double BB = e * Q2b[0] * sinphi;
  double CC = (D1 - D2) * p + e * cosphi * Q2b[0] - e * sinphi * Q2b[1];
  double DD = Q2b[0] - p * a1;
  double EE = -Q2b[1] + p * a2;
  double Poly4[5];
  Poly4[0] = AA - EE;
  Poly4[1] = -2 * CC + 2 * DD;
  Poly4[2] = 4 * BB - 2 * AA;
  Poly4[3] = 2 * CC + 2 * DD;
  Poly4[4] = AA + EE;
  int nbRealRacines;
  double Racines[4];
  double PossiblesTheta[8];

  compute_racines(Poly4, &nbRealRacines, Racines);
  for (int numR = 0; numR < nbRealRacines; numR++)
    PossiblesTheta[numR] = 2 * atan(Racines[numR]);
  /*Case RTb[x]==0*/
  PossiblesTheta[nbRealRacines] = -M_PI / 2;
  PossiblesTheta[nbRealRacines + 1] = M_PI / 2;
  PossiblesTheta[nbRealRacines + 2] = M_PI;
  PossiblesTheta[nbRealRacines + 3] = 0;
  for (int numR = 0; numR < nbRealRacines + 4; numR++)
  {
#ifdef FC3D_UE_DEBUG
    if (numR >= nbRealRacines)
      printf("FC3D_UE_DEBUG: The last attempt is R_T1 or R_T2 equal to 0?\n");
#endif
    //    double R=Racines[numR];
    double theta = PossiblesTheta[numR]; //2*atan(R);
    double costheta = cos(theta);
    double sintheta = sin(theta);
    double radius = p / (1 + e * cos(theta - phi));
    double fabsradius = fabs(radius);
    RTb[0] = radius * costheta;
    RTb[1] = radius * sintheta;
    *reaction1 = (*V00) * RTb[0] + (*V10) * RTb[1];
    *reaction2 = (*V01) * RTb[0] + (*V11) * RTb[1];
    *reaction0 = sqrt((*reaction1) * (*reaction1) + (*reaction2) * (*reaction2)) / (*mu);

    //In particular case RTb[x]==0, then check :
    if (numR >= nbRealRacines)
    {
      //RTb[0]==0
      if (numR < nbRealRacines + 2)
      {
        alpha = (-Q2b[1] + a2 * fabsradius) / RTb[1] - D2;
        if (fabs(Q2b[0] - a1 * fabsradius) > tol)
          continue;
      }
      else  //RTb[1]==0
      {
        alpha = (-Q2b[0] + a1 * fabsradius) / RTb[0] - D1;
        if (fabs(Q2b[1] - a2 * fabsradius) > tol)
          continue;
      }
    }
    else
    {
      alpha = (-Q2b[0] + a1 * fabsradius) / RTb[0] - D1;
    }

    if (alpha <= 0)
      continue;
#ifdef FC3D_UE_DEBUG
    double alpha1 = 0.0;
    double alpha2 = 0.0;
    if (numR >= nbRealRacines)
    {
      //RTb[0]==0
      if (numR < nbRealRacines + 2)
        alpha1 = (-Q2b[1] + a2 * fabsradius) / RTb[1] - D2;
      else
        alpha2 = (-Q2b[0] + a1 * fabsradius) / RTb[0] - D1;
    }
    else
    {
      if (RTb[1] != 0.0)
        alpha1 = (-Q2b[1] + a2 * fabsradius) / RTb[1] - D2;
      if (RTb[0] != 0.0)
        alpha2 = (-Q2b[0] + a1 * fabsradius) / RTb[0] - D1;
    }
    printf("FC3D_UE_DEBUG :: alpha1 = %e = %e =alpha2.(must be equal except if RTb[x]==0)\n", alpha1, alpha2);
    if (!FC3D_UE_TEST_NULL(e))
    {
      double zero1 = (D1 - D2) * radius * costheta * sintheta + Q2b[0] * sintheta - Q2b[1] * costheta - fabsradius * (a1 * sintheta - a2 * costheta);
      double zero2 = (D1 - D2) * RTb[0] * RTb[1] + Q2b[0] * RTb[1] - Q2b[1] * RTb[0] - (a1 * RTb[1] - a2 * RTb[0]) * sqrt(RTb[0] * RTb[0] + RTb[1] * RTb[1]);
      printf("FC3D_UE_DEBUG :: zero1 = %e = %e =zero2\n", zero1, zero2);
      double D_Dir2[2];
      D_Dir2[0] = (*V00) * D_Dir[0] + (*V01) * D_Dir[1];
      D_Dir2[1] = (*V10) * D_Dir[0] + (*V11) * D_Dir[1];
      double ONELIPSE2 = fabsradius / sqrt(
                           (D_Dir2[1] * RTb[0] - D_Dir2[0] * RTb[1] - OD2[0] * D_Dir2[1] + D_Dir2[0] * OD2[1]) *
                           (D_Dir2[1] * RTb[0] - D_Dir2[0] * RTb[1] - OD2[0] * D_Dir2[1] + D_Dir2[0] * OD2[1]) / (D_Dir2[0] * D_Dir2[0] + D_Dir2[1] * D_Dir2[1])) - e;
      double ONELIPSE = ((*M00) * (*M00) / ((*mu) * (*mu))) * radius * radius -
                        (Q[0] + (*M01) * (*reaction1) + (*M02) * (*reaction2)) *
                        (Q[0] + (*M01) * (*reaction1) + (*M02) * (*reaction2));
      printf("FC3D_UE_DEBUG :: OnElipse2 = %e; OnElipse = %e.(must be null)\n", ONELIPSE2, ONELIPSE);
    }
#endif
    double s1 = (*M00) * fabsradius / (*mu);
    double s2 = -Q[0] - (*M01) * (*reaction1) - (*M02) * (*reaction2);
#ifdef FC3D_UE_DEBUG
    printf("FC3D_UE_DEBUG: fabs(s1)=fabs(s2)=%e=%e\n", fabs(s1), fabs(s2));
#endif
    if ((s1 >= 0. && s2 >= 0.) || (s1 <= 0. && s2 <= 0.) || FC3D_UE_TEST_NULL(e))
    {
#ifdef FC3D_UE_DEBUG
      printf("R:\n");
      printf("%e\n%e\n%e\n", *reaction0, *reaction1, *reaction2);
      printf("-alphaRT:\n");
      printf("%e\n%e\n", -alpha* *reaction1, -alpha* *reaction2);
      printf("MR+q:\n");
      printf("%e\n%e\n%e\n",
             *M00 * (*reaction0) + *M01 * (*reaction1) + *M02 * (*reaction2) + *Q0,
             *M10 * (*reaction0) + *M11 * (*reaction1) + *M12 * (*reaction2) + *Q1,
             *M20 * (*reaction0) + *M21 * (*reaction1) + *M22 * (*reaction2) + *Q2);
      double err = 0.0;
      double velocity[3];

      velocity[0] = 0;
      velocity[1] = -alpha * (*reaction1);
      velocity[2] = -alpha * (*reaction2);
      fc3d_unitary_compute_and_add_error(reaction, velocity, *(problem->mu), &err);
      printf("Compute v with alpha=%e\n", alpha);
      printf("error is %e.", err);
#endif
      return 0;
    }
  }

#ifdef FC3D_UE_DEBUG
  printf("FC3D_UE_DEBUG: Solver failed\n");
#endif
  printf("fc3d_unitary_enumerative (failed) M:\n");
  print3x3(M);
  M = M00;
  printf("fc3d_unitary_enumerative (failed) Q:\n");
  print3(Q);
  Q = Q0;
  projectionOnCone(reaction, *mu);
  return -1;
}
int fc3d_unitary_enumerative_setDefaultSolverOptions(
  SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Setting default options for unitary enumerative solver.\n");
  }

  options->solverId =  SICONOS_FRICTION_3D_ONECONTACT_QUARTIC;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *)calloc(options->iSize, sizeof(int));
  options->dparam = (double *)calloc(options->dSize, sizeof(double));
  options->dWork = NULL;
  solver_options_nullify(options);
  options->dparam[0] = 1e-9;

  options->internalSolvers = NULL;

  return 0;
}
/*ax+by+c=0
  a1x+b1y+c1=0*/
void solve2x2(double *a, double *b, double *c, double *a1, double *b1, double *c1, double *x, double *y)
{
  double delta = *a * *b1 - *a1 * *b;
  if (delta == 0)
  {
    *x = NAN;
    *y = NAN;
  }
  else
  {
    double invd = 1.0 / delta;
    *x = (*b * *c1 - *c * *b1) * invd;
    *y = -(*a * *c1 - *a1 * *c) * invd;
  }

}
/*
 *Implementation from Gilles Davier : quatic formulation in respect of alpha.
 */
int fc3d_unitary_enumerative_solve_poly_nu_sliding(FrictionContactProblem* problem, double * reaction, SolverOptions* options)
{
  double * M = problem->M->matrix0;
  double * Q = problem->q;

  SET3X3(M);
  M = M00;
  SET3(reaction);
  reaction = reaction0;
  SET3(Q);
  Q = Q0;

  double cg[5];

  int nbRealRacines;
  double Racines[4];

  double W0 = *M00 ;
  double invW0 = 1 / W0;
  double alpha = *M01 ;
  double beta  = *M02 ;

  double q = *Q0 ;
  double q1 = -(- alpha * q * invW0 + *Q1);
  double q2 = -(- beta * q * invW0 + *Q2);


  double lambda1 = *M11 - alpha * alpha * invW0;
  double lambda2 = *M22 - beta * beta * invW0;
  double lambda12  = *M12 - alpha * beta * invW0;
  double mu = *(problem->mu);
#ifdef FC3D_UE_DEBUG
  printf("W0=%e\ninvW0=%e\n,alpha=%e\nbeta=%e\nlambda1=%e\nlambda2=%e\nlambda12=%e\nq=%e\nq1=%e\nq2=%e\nmu=%e", W0, invW0, alpha, beta, lambda1, lambda2, lambda12, q, q1, q2, mu);
#endif


  //MAPLE - 1contact.mws
  cg[0] = -mu * mu * q * q;
  cg[1] = -2 * mu * mu * q * alpha * q1 - 2 * mu * mu * q * beta * q2 - (2 * lambda1 + 2 * lambda2) * mu * mu * q * q;
  cg[2] = 2 * mu * mu * q * beta * q1 * lambda12 + W0 * W0 * q2 * q2 + 2 * mu * mu * q * q * lambda12 * lambda12 - 2 * mu * mu * beta * alpha * q1 * q2 + 2 * mu * mu * q * alpha * q2 * lambda12 - (lambda1 * lambda1 + 4 * lambda1 * lambda2 + lambda2 * lambda2) * mu * mu * q * q + W0 * W0 * q1 * q1 - beta * beta * mu * mu * q2 * q2 - 2 * (lambda1 + 2 * lambda2) * alpha * mu * mu * q * q1 - mu * mu * alpha * alpha * q1 * q1 - 2 * (2 * lambda1 + lambda2) * beta * mu * mu * q * q2;
  cg[3] = 2 * beta * mu * mu * q * q2 * lambda12 * lambda12 + 2 * (lambda1 + lambda2) * mu * mu * q * q * lambda12 * lambda12 - 4 * W0 * W0 * q1 * q2 * lambda12 - 2 * (lambda1 * lambda1 + 2 * lambda1 * lambda2) * beta * mu * mu * q * q2 + 2 * (lambda1 + lambda2) * alpha * mu * mu * q * q2 * lambda12 - (2 * lambda1 * lambda1 * lambda2 + 2 * lambda1 * lambda2 * lambda2) * mu * mu * q * q + 2 * (lambda1 + lambda2) * beta * mu * mu * q * q1 * lambda12 + 2 * alpha * beta * mu * mu * q2 * q2 * lambda12 + 2 * lambda1 * W0 * W0 * q2 * q2 + 2 * beta * beta * mu * mu * q1 * q2 * lambda12 - 2 * lambda2 * alpha * alpha * mu * mu * q1 * q1 - 2 * (lambda1 + lambda2) * alpha * beta * mu * mu * q1 * q2 - 2 * (2 * lambda1 * lambda2 + lambda2 * lambda2) * alpha * mu * mu * q * q1 + 2 * mu * mu * alpha * alpha * q1 * q2 * lambda12 + 2 * mu * mu * beta * alpha * q1 * q1 * lambda12 + 2 * mu * mu * q * alpha * q1 * lambda12 * lambda12 + 2 * lambda2 * W0 * W0 * q1 * q1 - 2 * lambda1 * beta * beta * mu * mu * q2 * q2;
  cg[4] = -beta * beta * mu * mu * q1 * q1 * lambda12 * lambda12 - alpha * alpha * mu * mu * q2 * q2 * lambda12 * lambda12 - lambda1 * lambda1 * beta * beta * mu * mu * q2 * q2 - lambda1 * lambda1 * lambda2 * lambda2 * mu * mu * q * q - lambda2 * lambda2 * alpha * alpha * mu * mu * q1 * q1 + lambda2 * lambda2 * W0 * W0 * q1 * q1 + lambda1 * lambda1 * W0 * W0 * q2 * q2 + W0 * W0 * q1 * q1 * lambda12 * lambda12 - 2 * beta * mu * mu * q * q1 *  pow((double) lambda12, (double) 3) - 2 * alpha * mu * mu * q * q2 *  pow((double) lambda12, (double) 3) + 2 * lambda1 * lambda2 * beta * mu * mu * q * q1 * lambda12 + 2 * lambda1 * lambda2 * alpha * mu * mu * q * q2 * lambda12 - 2 * lambda1 * lambda2 * alpha * beta * mu * mu * q1 * q2 - 2 * alpha * beta * mu * mu * q1 * q2 * lambda12 * lambda12 + 2 * lambda2 * alpha * beta * mu * mu * q1 * q1 * lambda12 + 2 * lambda2 * alpha * alpha * mu * mu * q1 * q2 * lambda12 + 2 * lambda1 * beta * mu * mu * q * q2 * lambda12 * lambda12 + 2 * lambda1 * beta * beta * mu * mu * q1 * q2 * lambda12 + 2 * lambda1 * alpha * beta * mu * mu * q2 * q2 * lambda12 - 2 * lambda1 * lambda2 * lambda2 * alpha * mu * mu * q * q1 - 2 * lambda1 * lambda1 * lambda2 * beta * mu * mu * q * q2 + 2 * lambda2 * alpha * mu * mu * q * q1 * lambda12 * lambda12 - mu * mu * q * q *  pow((double) lambda12, (double) 4) + W0 * W0 * q2 * q2 * lambda12 * lambda12 - 2 * lambda2 * W0 * W0 * q1 * q2 * lambda12 - 2 * lambda1 * W0 * W0 * q1 * q2 * lambda12 + 2 * lambda1 * lambda2 * mu * mu * q * q * lambda12 * lambda12;

  compute_racines(cg, &nbRealRacines, Racines);
  double nu = -1;
  for (int i = 0; i < nbRealRacines; i++)
  {
    if (nu < 0 && Racines[i] > 0)
      nu = Racines[i];
    else if (Racines[i] > 0 && Racines[i] < nu)
      nu = Racines[i];
  }
  if (nu > 0)
  {
#ifdef FC3D_UE_DEBUG
    printf("nb racines :%d\n", nbRealRacines);
    printf("nu=%e\n", nu);
#endif
    /*system Mt*Rt=qt*/
    double Mt00 = lambda1 + nu;
    double Mt01 = lambda12;
    double Mt10 = lambda12;
    double Mt11 = lambda2 + nu;
    double qt0 = -q1;
    double qt1 = -q2;
    solve2x2(&Mt00, &Mt01, &qt0, &Mt10, &Mt11, &qt1, reaction1, reaction2);
    *reaction0 = sqrt((*reaction1 * *reaction1) + (*reaction2 * *reaction2)) / mu;

    reaction = reaction0;

#ifdef FC3D_UE_DEBUG
    double velocity_[3];
    double * velocity = velocity_;
    SET3(velocity);
    velocity = velocity0;
    mv3x3(M, reaction, velocity);
    velocity = velocity0;
    reaction = reaction0;

    *velocity0 += *Q0;
    *velocity1 += *Q1;
    *velocity2 += *Q2;
    printf("reaction is:\n");
    print3(reaction);
    reaction = reaction0;
    printf("velocity is:\n");
    print3(velocity);
    velocity = velocity0;
    double err = 0.0;
    fc3d_unitary_compute_and_add_error(reaction, velocity, mu, &err);
    printf("error is %e.\n", err);
    printf("fc3d_unitary_enumerative_poly_nu  M:\n");
    print3x3(M);
    M = M00;
    printf("fc3d_unitary_enumerative_poly_nu  Q:\n");
    print3(Q);
    Q = Q0;
#endif
    return 0;
  }
#ifdef FC3D_UE_DEBUG

#endif
  if (verbose) printf("fc3d_unitary_enumerative_poly_nu (FAILED) M:\n");
  if (verbose) print3x3(M);
  M = M00;
  if (verbose) printf("fc3d_unitary_enumerative_poly_nu (FAILED) Q:\n");
  if (verbose) print3(Q);
  Q = Q0;
  return -1;
}
