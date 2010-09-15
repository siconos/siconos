#define _BSD_SOURCE

#include "LA.h"
#include "FrictionContact3D_Solvers.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "Friction_cst.h"
#include "op3x3.h"
#include "FrictionContact3D_unitary_enumerative.h"
//#define FC3D_UE_DEBUG

#include <stdio.h>
#include <math.h>
#include "quartic.h"


#define FC3D_UE_TEST_NULL(EXPR)  (EXPR==0)

void compute_racines(double * Poly, int *nbRealRacines, double *Racines)
{
  double r[3][5];
  //Roots of poly p[0] x^4 + p[1] x^3...+p[4]=0
  //x=r[1][k] + i r[2][k]  k=1,...,4
#ifdef FC3D_UE_DEBUG
  double Psav[5];
  for (int k = 0; k < 5; k++)
    Psav[k] = Poly[k];
#endif
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
    printf("FC3D_unitary_enumerative: degre of polynom is 0.");
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
  if (*c == 0)
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
void frictionContact3D_unitary_enumerative_free(FrictionContactProblem* localproblem)
{
  free(localproblem->M->matrix0);
  localproblem->M->matrix0 = NULL;
}
void frictionContact3D_unitary_enumerative_initialize(FrictionContactProblem* localproblem)
{
  if (!localproblem->M->matrix0)
    localproblem->M->matrix0 = (double*)malloc(9 * sizeof(double));

}
int frictionContact3D_unitary_enumerative_solve(FrictionContactProblem* problem, double * reaction, SolverOptions* options)
{
  double velocity[3];
  int info;
#ifdef FC3D_UE_DEBUG
  printf("frictionContact3D_unitary_enumerative_solve: begin\n");
#endif
  frictionContact3D_unitary_enumerative(problem, reaction, velocity, &info, options);
#ifdef FC3D_UE_DEBUG
  printf("frictionContact3D_unitary_enumerative_solve: end\n");
#endif
  return info;
}
int frictionContact3D_unitary_enumerative(FrictionContactProblem* problem, double * reaction, double * velocity, int *info, SolverOptions* options)
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
  SET3(velocity);
  velocity = velocity0;
  SET3(Q);
  Q = Q0;
  (*info) = -1;
  //printf("frictionContact3D_unitary_enumerative M:\n");
  //print3x3(M);M=M00;
  //printf("frictionContact3D_unitary_enumerative Q:\n");
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
    (*info) = 0;
    return 0;
  }

  if (*mu == 0.)
  {
    *reaction0 = -(*Q0) / (*M00);
    if (*reaction0 < 0)
      return -1;
    *reaction1 = 0;
    *reaction2 = 0;
    *velocity0 = 0;
    *velocity1 = *reaction0 * *M10 + *Q1;
    *velocity2 = *reaction0 * *M20 + *Q2;
    return 0;
  }

  /*sticking ? 0=MR+q*/
  solv3x3(M, reaction, Q);
  M = M00;
  reaction = reaction0;
  Q = Q0;
  (*info) = -1;
  if (!isnan(*reaction0))
  {
    if (- *reaction0 + tol > 0.0)
      if ((*reaction0 * *mu) * (*reaction0 * *mu) + tol > *reaction1 * *reaction1 + *reaction2 * *reaction2)
      {
        *reaction0 = - *reaction0;
        *reaction1 = - *reaction1;
        *reaction2 = - *reaction2;
        *velocity0 = 0;
        *velocity1 = 0;
        *velocity2 = 0;
        (*info) = 0;
        return 0;
      }
  }

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
  double OD[2];
  double OD2[2];
  double D_Dir[2];
  double Q2b[2];
  double M1_b[2];
  double NormD2 = (*M01) * (*M01) + (*M02) * (*M02);
  double NormD = sqrt(NormD2);
  double e = (*mu) * NormD / (*M00);
  double d, p;
  double phi;
  double cosphi;
  double sinphi;
  FC3D_unitary_enum_factorize2x2(M11, M22, M12, &D1, &D2, V);
  if (!FC3D_UE_TEST_NULL(e))
  {
    d = fabs(*Q0) / NormD;
    p = e * d;
    OD[0] = -((*M01) * (*Q)) / NormD2;
    OD[1] = -(((-*M01) * (*M01) * (*Q0)) / NormD2 + *Q0) / (*M02);
    D_Dir[0] = 1;
    D_Dir[1] = -(*M01) / (*M02);
    OD2[0] = (*V00) * OD[0] + (*V01) * OD[1];
    OD2[1] = (*V10) * OD[0] + (*V11) * OD[1];
    phi = atan(OD2[1] / OD2[0]);
    if (OD2[0] < 0)
      phi += M_PI;
    cosphi = cos(phi);
    sinphi = sin(phi);
  }
  else
  {
    d = 0;
    p = (-(*Q0) * (*mu)) / (*M00);
    if (p < 0)
    {
      return -1;
    }
  }
  double RTb[2];

  Q2b[0] = (*V00) * Q_2[0] + (*V01) * Q_2[1];
  Q2b[1] = (*V10) * Q_2[0] + (*V11) * Q_2[1];
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
  double PossiblesTheta[6];

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
      //RTb[0]==0
      if (numR < nbRealRacines + 2)
        alpha1 = (-Q2b[1] + a2 * fabsradius) / RTb[1] - D2;
      else
        alpha2 = (-Q2b[0] + a1 * fabsradius) / RTb[0] - D1;
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
#endif
      *velocity0 = 0;
      *velocity1 = -alpha * (*reaction1);
      *velocity2 = -alpha * (*reaction2);
      (*info) = 0;
      return 0;
    }
  }

#ifdef FC3D_UE_DEBUG
  printf("FC3D_UE_DEBUG: Solver failed\n");
#endif
  return -1;
}
int frictionContact3D_unitary_enumerative_setDefaultSolverOptions(
  SolverOptions* options)
{
  if (verbose > 0)
  {
    printf("Set the default solver options for the GLOBALAC Solver\n");
  }

  options->solverId =  SICONOS_FRICTION_3D_QUARTIC;
  options->numberOfInternalSolvers = 0;
  options->isSet = 1;
  options->filterOn = 1;
  options->iSize = 5;
  options->dSize = 5;
  options->iparam = (int *) malloc(options->iSize * sizeof(int));
  options->dparam = (double *) malloc(options->dSize * sizeof(double));
  options->dWork = NULL;
  options->iWork = NULL;
  for (unsigned int i = 0; i < 5; i++)
  {
    options->iparam[i] = 0;
    options->dparam[i] = 0.0;
  }
  options->dparam[0] = 1e-9;

  options->internalSolvers = NULL;

  return 0;
}
