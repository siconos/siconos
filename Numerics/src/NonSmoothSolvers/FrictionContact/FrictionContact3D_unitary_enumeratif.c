#include "LA.h"
#include "FrictionContact3D_Solvers.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "Friction_cst.h"
#include "op3x3.h"

#define FC3D_UE_DEBUG


#include <stdio.h>
#include <math.h>

//ROOTS OF LOW-ORDER POLYNOMIAL
//TERENCE R. F. NONWEILER (Recd. 14 Apr. 1967)
//James Watt Engineering Laboratories, The University,
//Glasgow W2, Scotland

int QUADROOTS(p, r)
double p[5], r[3][5];
{
  /*
  Array r[3][5]  p[5]
  Roots of poly p[0] x^2 + p[1] x+p[2]=0
  x=r[1][k] + i r[2][k]  k=1,2
  */
  double b, c, d;
  b = -p[1] / p[0] / 2;
  c = p[2] / p[0];
  d = b * b - c;
  if (d > 0)
  {
    if (b > 0) b = r[1][2] = sqrt(d) + b;
    else    b = r[1][2] = -sqrt(d) + b;
    r[1][1] = c / b;
    r[2][1] = r[2][2] = 0;
  }
  else
  {
    d = r[2][1] = sqrt(-d);
    r[2][2] = -d;
    r[1][1] = r[1][2] = b;
  }
  return(0);
}
int CUBICROOTS(p, r)
double p[5], r[3][5];
{
  /*
  Array r[3][5]  p[5]
  Roots of poly p[0] x^3 + p[1] x^2...+p[3]=0
  x=r[1][k] + i r[2][k]  k=1,...,3
  Assumes 0<arctan(x)<pi/2 for x>0
  */

  double s, t, b, c, d;
  int k;
  if (p[0] != 1)
    for (k = 1; k < 4; k++) p[k] = p[k] / p[0];
  p[0] = 1;
  s = p[1] / 3.0;
  t = s * p[1];
  b = 0.5 * (s * (t / 1.5 - p[2]) + p[3]);
  t = (t - p[2]) / 3.0;
  c = t * t * t;
  d = b * b - c;
  if (d >= 0)
  {
    d = pow((sqrt(d) + fabs(b)), 1.0 / 3.0);
    printf("d=%f\n", d);
    if (d != 0)
    {
      if (b > 0) b = -d;
      else b = d;
      c = t / b;
    }
    d = r[2][2] = sqrt(0.75) * (b - c);
    b = b + c;
    c = r[1][2] = -0.5 * b - s;
    if ((b > 0 && s <= 0) || (b < 0 && s > 0))
    {
      r[1][1] = c;
      r[2][1] = -d;
      r[1][3] = b - s;
      r[2][3] = 0;
    }
    else
    {
      r[1][1] = b - s;
      r[2][1] = 0;
      r[1][3] = c;
      r[2][3] = -d;
    }
  }  /* end 2 equal or complex roots */
  else
  {
    if (b == 0)
      d = atan(1.0) / 1.5;
    else
      d = atan(sqrt(-d) / fabs(b)) / 3.0;
    if (b < 0)
      b = sqrt(t) * 2.0;
    else
      b = -2.0 * sqrt(t);
    c = cos(d) * b;
    t = -sqrt(0.75) * sin(d) * b - 0.5 * c;
    d = -t - c - s;
    c = c - s;
    t = t - s;
    if (fabs(c) > fabs(t))
      r[1][3] = c;
    else
    {
      r[1][3] = t;
      t = c;
    }
    if (fabs(d) > fabs(t))
      r[1][2] = d;
    else
    {
      r[1][2] = t;
      t = d;
    }
    r[1][1] = t;
    for (k = 1; k < 4; k++) r[2][k] = 0;
  }
  return(0);
}
int BIQUADROOTS(p, r)
/* add _ if calling from fortran */
/*
Array r[3][5]  p[5]
Roots of poly p[0] x^4 + p[1] x^3...+p[4]=0
x=r[1][k] + i r[2][k]  k=1,...,4
*/

double p[5], r[3][5];
{
  double a, b, c, d, e;
  int k, j;
  if (p[0] != 1.0)
  {
    for (k = 1; k < 5; k++) p[k] = p[k] / p[0];
    p[0] = 1;
  }
  e = 0.25 * p[1];
  b = 2 * e;
  c = b * b;
  d = 0.75 * c;
  b = p[3] + b * (c - p[2]);
  a = p[2] - d;
  c = p[4] + e * (e * a - p[3]);
  a = a - d;
  p[1] = 0.5 * a;
  p[2] = (p[1] * p[1] - c) * 0.25;
  p[3] = b * b / (-64.0);
  if (p[3] < -1e-6)
  {
    CUBICROOTS(p, r);
    for (k = 1; k < 4; k++)
    {
      if (r[2][k] == 0 && r[1][k] > 0)
      {
        d = r[1][k] * 4;
        a = a + d;
        if (a >= 0 && b >= 0)
          p[1] = sqrt(d);
        else if (a <= 0 && b <= 0)
          p[1] = sqrt(d);
        else p[1] = -sqrt(d);
        b = 0.5 * (a + b / p[1]);
        goto QUAD;
      }
    }
  }
  if (p[2] < 0)
  {
    b = sqrt(c);
    d = b + b - a;
    p[1] = 0;
    if (d > 0) p[1] = sqrt(d);
  }
  else
  {
    if (p[1] > 0)
      b = sqrt(p[2]) * 2.0 + p[1];
    else
      b = -sqrt(p[2]) * 2.0 + p[1];
    if (b != 0)
      p[1] = 0;
    else
    {
      for (k = 1; k < 5; k++)
      {
        r[1][k] = -e;
        r[2][k] = 0;
      }
      goto END;
    }
  }
QUAD:
  p[2] = c / b;
  QUADROOTS(p, r);
  for (k = 1; k < 3; k++)
    for (j = 1; j < 3; j++) r[j][k + 2] = r[j][k];
  p[1] = -p[1];
  p[2] = b;
  QUADROOTS(p, r);
  for (k = 1; k < 5; k++) r[1][k] = r[1][k] - e;
END:
  ;
  return(0);
}


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
  BIQUADROOTS(Poly, r);
  (*nbRealRacines) = 0;
  for (int k = 1; k < 5; k++)
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
  |ca|
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
  if (c == 0)
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

int frictionContact3D_unitary_enumeratif(FrictionContactProblem* problem, double * reaction, double * velocity, int *info, SolverOptions* options)
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
  /*Decollage? R=0 ?*/
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

  /*Adherance? 0=MR+q*/
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

  /*Glissement? */
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
  double d = fabs(*Q0) / NormD;
  double p = e * d;
  double RTb[2];

  OD[0] = -((*M01) * (*Q)) / NormD2;
  OD[1] = -(((-*M01) * (*M01) * (*Q0)) / NormD2 + *Q0) / (*M02);
  D_Dir[0] = 1;
  D_Dir[1] = -(*M01) / (*M02);
  FC3D_unitary_enum_factorize2x2(M11, M22, M12, &D1, &D2, V);
  OD2[0] = (*V00) * OD[0] + (*V01) * OD[1];
  OD2[1] = (*V10) * OD[0] + (*V11) * OD[1];
  double phi = atan(OD2[1] / OD2[0]);
  if (OD2[0] < 0)
    phi += 3.14159265358979323846;
  double cosphi = cos(phi);
  double sinphi = sin(phi);
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

  compute_racines(Poly4, &nbRealRacines, Racines);

  for (int numR = 0; numR < nbRealRacines; numR++)
  {
    double R = Racines[numR];
    double theta = 2 * atan(R);
    double costheta = cos(theta);
    double sintheta = sin(theta);
    double radius = p / (1 + e * cos(theta - phi));
    RTb[0] = radius * costheta;
    RTb[1] = radius * sintheta;
    *reaction1 = (*V00) * RTb[0] + (*V10) * RTb[1];
    *reaction2 = (*V01) * RTb[0] + (*V11) * RTb[1];
    *reaction0 = sqrt((*reaction1) * (*reaction1) + (*reaction2) * (*reaction2));
    alpha = (-Q2b[1] + a2 * radius) / RTb[1] - D2;
    if (alpha <= 0)
      continue;
#ifdef FC3D_UE_DEBUG
    double alpha1 = (-Q2b[1] + a2 * radius) / RTb[1] - D2;
    double alpha2 = (-Q2b[0] + a1 * radius) / RTb[0] - D1;
    printf("FC3D_UE_DEBUG :: alpha1 = %e = %e =alpha2\n", alpha1, alpha2);
    double zero1 = (D1 - D2) * radius * costheta * sintheta + Q2b[0] * sintheta - Q2b[1] * costheta - radius * (a1 * sintheta - a2 * costheta);
    double zero2 = (D1 - D2) * RTb[0] * RTb[1] + Q2b[0] * RTb[1] - Q2b[1] * RTb[0] - (a1 * RTb[1] - a2 * RTb[0]) * sqrt(RTb[0] * RTb[0] + RTb[1] * RTb[1]);
    printf("FC3D_UE_DEBUG :: zero1 = %e = %e =zero2\n", zero1, zero2);
    double D_Dir2[2];
    D_Dir2[0] = (*V00) * D_Dir[0] + (*V01) * D_Dir[1];
    D_Dir2[1] = (*V10) * D_Dir[0] + (*V11) * D_Dir[1];
    double ONELIPSE2 = radius / sqrt(
                         (D_Dir2[1] * RTb[0] - D_Dir2[0] * RTb[1] - OD2[0] * D_Dir2[1] + D_Dir2[0] * OD2[1]) *
                         (D_Dir2[1] * RTb[0] - D_Dir2[0] * RTb[1] - OD2[0] * D_Dir2[1] + D_Dir2[0] * OD2[1]) / (D_Dir2[0] * D_Dir2[0] + D_Dir2[1] * D_Dir2[1])) - e;
    double ONELIPSE = ((*M00) * (*M00) / ((*mu) * (*mu))) * radius * radius -
                      (Q[0] + (*M01) * (*reaction1) + (*M02) * (*reaction2)) *
                      (Q[0] + (*M01) * (*reaction1) + (*M02) * (*reaction2));
    printf("FC3D_UE_DEBUG :: OnElipse2 = %e; OnElipse = %e.(must be null)\n", ONELIPSE2, ONELIPSE);
#endif
    double s1 = (*M00) * radius / (*mu);
    double s2 = -Q[0] - (*M01) * (*reaction1) - (*M02) * (*reaction2);
#ifdef FC3D_UE_DEBUG
    printf("fabs(s1)=fabs(s2)=%e=%e\n", fabs(s1), fabs(s2));
#endif
    if ((s1 >= 0. && s2 >= 0.) || (s1 <= 0. && s2 <= 0.))
    {
      *velocity0 = 0;
      *velocity1 = -alpha * (*reaction1);
      *velocity2 = -alpha * (*reaction2);
      (*info) = 0;
      return 0;
    }
  }
  return -1;
}
int frictionContact3D_unitary_enumeratif_setDefaultSolverOptions(
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
