#ifndef AlartCurnierGenerated_h
#define AlartCurnierGenerated_h
#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
extern "C"
{
#endif

#include "fc3d_AlartCurnierABGenerated.h"
#include "fc3d_AlartCurnierFABGenerated.h"
#include "fc3d_AlartCurnierFGenerated.h"
#include "fc3d_AlartCurnierJeanMoreauABGenerated.h"
#include "fc3d_AlartCurnierJeanMoreauFABGenerated.h"
#include "fc3d_AlartCurnierJeanMoreauFGenerated.h"


void fc3d_AlartCurnierFunctionGenerated(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *f,
  double *A,
  double *B);

void fc3d_AlartCurnierJeanMoreauFunctionGenerated(
  double *reaction,
  double *velocity,
  double mu,
  double *rho,
  double *f,
  double *A,
  double *B);

#if defined(__cplusplus) && !defined(BUILD_AS_CPP)
}
#endif
#endif
