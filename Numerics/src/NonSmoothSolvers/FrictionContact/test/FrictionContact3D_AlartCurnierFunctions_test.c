#undef NDEBUG
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "NumericsConfig.h"
#ifdef WITH_TIMERS
#define TIMER_FFTW_CYCLE
#endif
#include "timers_interf.h"
#include "op3x3.h"
#include "NonSmoothDrivers.h"
#include "FrictionContact3D_AlartCurnier.h"
#include "FrictionContact3D_globalAlartCurnier.h"

#define SIZE 1000

#define DO(X)                                               \
  do {for (unsigned int i=0;i<SIZE;++i) {X;};} while(0)

int main()
{
  DECL_TIMER(T0);
  DECL_TIMER(T1);


  int info = 0;
  int r = -1;

  FILE* file = fopen("./data/ACinputs.dat", "r");
  int dim = 0;
  double* reactions;
  double* velocities;
  double *mus;
  double *rhos;

  r = fscanf(file, "%d\n", &dim);
  assert(r > 0);
  if (r <= 0) return(r);

  reactions = (double *) malloc(3 * dim * sizeof(double));
  velocities = (double *) malloc(3 * dim * sizeof(double));
  mus = (double *) malloc(dim * sizeof(double));
  rhos = (double *) malloc(3 * dim * sizeof(double));

  for (unsigned int i = 0; i < dim * 3 ; ++i)
  {
    r = fscanf(file, "%lf\n", &reactions[i]);
    assert(r > 0);
  };

  for (unsigned int i = 0; i < dim * 3 ; ++i)
  {
    r = fscanf(file, "%lf\n", &velocities[i]);
    assert(r > 0);
  };

  for (unsigned int k = 0; k < dim ; ++k)
  {
    r = fscanf(file, "%lf\n", &mus[k]);
    assert(r > 0);
  };

  for (unsigned int i = 0; i < dim * 3 ; ++i)
  {
    r = fscanf(file, "%lf\n", &rhos[i]);
    assert(r > 0);
  };

  double F1[3], A1[9], B1[9],
         F2[3], A2[9], B2[9];
  for (unsigned int k = 0; k < dim; ++k)
  {

    double* p;

    p = F1;
    OP3(*p++ = NAN);

    p = F2;
    OP3(*p++ = NAN);

    p = A1;
    OP3X3(*p++ = NAN);

    p = B1;
    OP3X3(*p++ = NAN);

    p = B2;
    OP3X3(*p++ = NAN);

    p = A2;
    OP3X3(*p++ = NAN);

    START_TIMER(T0);
    DO(computeAlartCurnierSTDOld(&reactions[k * 3], &velocities[k * 3], mus[k], &rhos[k * 3], F1, A1, B1));
    STOP_TIMER(T0);

    START_TIMER(T1);
    DO(computeAlartCurnierSTD(&reactions[k * 3], &velocities[k * 3], mus[k], &rhos[k * 3], F1, A1, B1));
    STOP_TIMER(T1);

    PRINT_ELAPSED(T0);

    PRINT_ELAPSED(T1);

#ifdef WITH_TIMERS
    printf("T1/T0 = %g\n", ELAPSED(T1) / ELAPSED(T0));
#endif

    p = F1;
    OP3(info |= isnan(*p++));
    assert(!info);

    p = A1;
    OP3X3(info |= isnan(*p++));
    assert(!info);

    p = B1;
    OP3X3(info |= isnan(*p++));
    assert(!info);

    frictionContact3D_localAlartCurnierFunctionGenerated(&reactions[k * 3], &velocities[k * 3], mus[k], &rhos[k * 3], F2, A2, B2);

    p = F1;
    OP3(info |= isnan(*p++));
    assert(!info);

    p = A1;
    OP3X3(info |= isnan(*p++));
    assert(!info);

    p = B1;
    OP3X3(info |= isnan(*p++));
    assert(!info);

    sub3(F1, F2);
    sub3x3(A1, A2);
    sub3x3(B1, B2);

#define EPS 1e-6
    p = F2;
    OP3(info |= !(*p++ < EPS));
    assert(!info);

    p = A2;
    OP3X3(info |= !(*p++ < EPS));
    assert(!info);

    p = B2;
    OP3X3(info |= !(*p++ < EPS));
    assert(!info);

  }

  free(reactions);
  free(velocities);
  free(mus);
  free(rhos);

  return (info);
}

