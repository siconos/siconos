

#include "SiconosNumerics.h"

#ifdef __cplusplus

#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#define restrict __restrict

#else

#ifdef _WIN32
#define SICONOS_EXPORT __declspec(dllexport)
#else
#define SICONOS_EXPORT
#endif

#endif

#include <stdio.h>

typedef struct
{
  int id;
  double* xk;
  double h;
  double theta;
  double gamma;
  double g;
  double kappa;
  unsigned int f_eval;
  unsigned int nabla_eval;
} data;

SICONOS_EXPORT void compute_F(void* problem, int n, double* restrict l, double* restrict F)
{
  data* d = (data*) VI_get_env(problem);

  double xk0 = d->xk[0];
  double xk1 = d->xk[1];
  double l0 = l[0];
  double l1 = l[1];
  double invR = 1.0/(1.0 - d->kappa*l0*l1);
  F[1] = d->h*d->g*l0*invR;

  double v_gamma = xk1 + d->gamma*F[1];
  F[0] = -d->h*d->kappa*l0*l1*v_gamma;

  double v_theta = xk1 + d->theta*F[1];
  F[0] += xk0 + d->h*v_theta;
  F[1] += xk1;
  d->f_eval += 1;

}

SICONOS_EXPORT void compute_nabla_F(void* problem, int n, double* restrict l, NumericsMatrix* restrict nabla_F_mat)
{
  data* d = (data*) VI_get_env(problem);
  double* restrict nabla_F = nabla_F_mat->matrix0;

  double l0 = l[0];
  double l1 = l[1];
  double xk1 = d->xk[1];
  double invR = 1.0/(1.0 - d->kappa*l0*l1);
  double invR2 = invR*invR;
  double rr1 = d->g*l0*invR;

  double v_gamma = xk1 + d->gamma*d->h*rr1;

  nabla_F[1] = d->h*d->g*invR2;
  nabla_F[1 + 2] = d->h*(d->g*d->kappa*l0*l0)/invR2;
//  nabla_F[0] = d->h*(-d->kappa*l1*v_gamma + d->gamma*d->h*(1.0-d->kappa*l0*l1)*nabla_F[1]);
//  nabla_F[0 + 2] = d->h*(-d->kappa*l0*v_gamma + d->gamma*d->h*(1.0-d->kappa*l0*l1)*nabla_F[1 + 2]);
  nabla_F[0] = d->h*(-d->kappa*l1*v_gamma + d->gamma*d->h*(1.0-d->kappa*l0*l1)*nabla_F[1]);
  nabla_F[0 + 2] = d->h*(-d->kappa*l0*v_gamma + d->gamma*d->h*(1.0-d->kappa*l0*l1)*nabla_F[1 + 2]);
  d->nabla_eval += 1;
}

SICONOS_EXPORT void compute_Fmcp(void* env, int n1, int n2, double* restrict z, double* restrict F)
{
  data* d = (data*) env;
  double l0 = 2.0*z[0] - 1.0;
  double l1 = 2.0*z[2] - 1.0;
  double r1 = d->g*l0/(1.0 - d->kappa*l0*l1);
  double v_gamma = (d->xk[1] + d->gamma*(d->h*r1));
  double r0 = -d->kappa*l0*l1*(v_gamma);
  double v_theta = d->xk[1] + d->theta*(d->h*r1);
  F[0] = d->xk[0] + d->h*v_theta + d->h*r0 + z[1];
  F[2] = d->xk[1] + d->h*r1 + z[3];
  F[1] = 1.0 - z[0];
  F[3] = 1.0 - z[2];
  d->f_eval += 1;
}

SICONOS_EXPORT void compute_nabla_Fmcp(void* env, int n1, int n2, double* restrict z, double* restrict nabla_Fmcp)
{
  data* d = (data*) env;
  double l0 = 2.0*z[0] - 1.0;
  double l1 = 2.0*z[2] - 1.0;
  double invR = 1.0/(1.0 - d->kappa*l0*l1);
  double invR2 = invR*invR;
  double r1 = d->g*l0*invR;
  double v_gamma = d->xk[1] + d->gamma*(d->h*r1);

  nabla_Fmcp[2] = 2.0*d->h*d->g*invR2;
  nabla_Fmcp[2 + 4] = 0.0;
  nabla_Fmcp[2 + 2*4] = 2.0*d->h*(d->g*d->kappa*l0*l0)*invR2;
  nabla_Fmcp[2 + 3*4] = 1.0;

  nabla_Fmcp[0] = -2.0*d->h*d->kappa*l1*(v_gamma) + d->h*(d->theta - d->gamma*d->kappa*l0*l1)*nabla_Fmcp[2];
  nabla_Fmcp[0 + 1*4] = 1.0;
  nabla_Fmcp[0 + 2*4] = -2.0*d->h*d->kappa*l0*(v_gamma) + d->h*(d->theta - d->gamma*d->kappa*l0*l1)*nabla_Fmcp[2 + 2*4];
  nabla_Fmcp[0 + 3*4] = 0.0;

  nabla_Fmcp[1] = -1.0;
  nabla_Fmcp[1 + 4] = 0.0;
  nabla_Fmcp[1 + 8] = 0.0;
  nabla_Fmcp[1 + 12] = 0.0;

  nabla_Fmcp[3] = 0.0;
  nabla_Fmcp[3 + 4] = 0.0;
  nabla_Fmcp[3 + 2*4] = -1.0;
  nabla_Fmcp[3 + 3*4] = 0.0;
  d->nabla_eval += 1;
}

SICONOS_EXPORT void compute_nabla_Fmcp2(void* env, int n1, int n2, double* restrict z, double* restrict nabla_Fmcp)
{
  data* d = (data*) env;
  double l0 = 2.0*z[0] - 1.0;
  double l1 = 2.0*z[2] - 1.0;
  double R = (1.0 - d->kappa*l0*l1);
  double R2 = R*R;
  double r1 = d->g*l0/R;
  double v_gamma = d->xk[1] + d->gamma*(d->h*r1);

  nabla_Fmcp[2] = 2.0*d->h*d->g/R2;
  nabla_Fmcp[2 + 4] = 0.0;
  nabla_Fmcp[2 + 2*4] = 2.0*d->h*(d->g*d->kappa*l0*l0)/R2;
  nabla_Fmcp[2 + 3*4] = 1.0;

  nabla_Fmcp[0] = -2.0*d->h*d->kappa*l1*(v_gamma) + d->h*(d->theta - d->gamma*d->kappa*l0*l1)*nabla_Fmcp[2];
  nabla_Fmcp[0 + 1*4] = 1.0;
  nabla_Fmcp[0 + 2*4] = -2.0*d->h*d->kappa*l0*(v_gamma) + d->h*(d->theta - d->gamma*d->kappa*l0*l1)*nabla_Fmcp[2 + 2*4];
  nabla_Fmcp[0 + 3*4] = 0.0;

  nabla_Fmcp[1] = -1.0;
  nabla_Fmcp[1 + 4] = 0.0;
  nabla_Fmcp[1 + 8] = 0.0;
  nabla_Fmcp[1 + 12] = 0.0;

  nabla_Fmcp[3] = 0.0;
  nabla_Fmcp[3 + 4] = 0.0;
  nabla_Fmcp[3 + 2*4] = -1.0;
  nabla_Fmcp[3 + 3*4] = 0.0;
  d->nabla_eval += 1;
}
