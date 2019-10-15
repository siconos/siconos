#undef NDEBUG
#include "SiconosConfig.h" // for WITH_TIMERS // IWYU pragma: keep
#ifdef WITH_TIMERS
#define TIMER_FFTW_CYCLE
#endif
#include <assert.h>                       // for assert
#include <math.h>                         // for NAN, isnan
#include <stdio.h>                        // for fscanf, fclose, fopen, FILE
#include <stdlib.h>                       // for free, malloc
#include "AlartCurnierGenerated.h"        // for fc3d_AlartCurnierFunctionGe...
#include "fc3d_AlartCurnier_functions.h"  // for computeAlartCurnierSTD
#include "op3x3.h"                        // for OP3X3, sub3x3, OP3, sub3
#include "timers_interf.h"                // for DECL_TIMER, PRINT_ELAPSED


void computeAlartCurnierSTDOld(double R[3], double velocity[3], double mu, double rho[3], double F[3], double A[9], double B[9]);

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
  unsigned int dim = 0;
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

    fc3d_AlartCurnierFunctionGenerated(&reactions[k * 3], &velocities[k * 3], mus[k], &rhos[k * 3], F2, A2, B2);

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

  fclose(file);
  return (info);
}

