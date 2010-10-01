
/* to be compared with :
/* valgrind --tool=callgrind --dump-instr=yes --dump-line=yes --collect-jumps=yes  --simulate-cache=yes */

#define TIMER_FFTW_CYCLE

#include <stdlib.h>
#include <math.h>
#include "timers_interf.h"


#define SIZE 1000000

#define DO(X)                                               \
  do {for (unsigned int i=0;i<SIZE;++i) {X;};} while(0)

int main()
{
  double *a;
  double *b;
  double *c;

  a = (double *) malloc(SIZE * sizeof(double));
  b = (double *) malloc(SIZE * sizeof(double));
  c = (double *) malloc(SIZE * sizeof(double));

  DECL_TIMER(T0);
  DECL_TIMER(T1);
  DECL_TIMER(T2);
  DECL_TIMER(T3);
  DECL_TIMER(T4);
  DECL_TIMER(T5);

  START_TIMER(T0);

  START_TIMER(T1);
  DO(a[i] = 1.0);
  STOP_TIMER(T1);

  DO(b[i] = 2.0);
  DO(c[i] = 0.);

  START_TIMER(T2);
  DO(c[i] = a[i] + b[i]);
  STOP_TIMER(T2);

  START_TIMER(T3);

  DO(c[i] = a[i] * b[i]);

  STOP_TIMER(T3);

  START_TIMER(T4);
  DO(c[i] = a[i] * a[i]);
  STOP_TIMER(T4);

  START_TIMER(T5);
  DO(c[i] = pow(a[i], 2.));
  STOP_TIMER(T5);

  STOP_TIMER(T0);

  PRINT_ELAPSED(T0);

  PRINT_ELAPSED(T1);

  PRINT_ELAPSED(T2);

  PRINT_ELAPSED(T3);

  PRINT_ELAPSED(T4);

  PRINT_ELAPSED(T5);

#ifdef WITH_TIMERS
  printf("T1/T0 = %g\n", ELAPSED(T1) / ELAPSED(T0));
  printf("T2/T0 = %g\n", ELAPSED(T2) / ELAPSED(T0));
  printf("T3/T0 = %g\n", ELAPSED(T3) / ELAPSED(T0));
  printf("T4/T0 = %g\n", ELAPSED(T4) / ELAPSED(T0));
  printf("T5/T0 = %g\n", ELAPSED(T5) / ELAPSED(T0));
#endif

  exit(0);
}

